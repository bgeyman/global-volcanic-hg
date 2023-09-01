# system for flexibly defining and running environmental compartment models
# Author contact: Colin Thackray (thackray@seas.harvard.edu)

"""
Copyright (c) 2019 Colin Thackray.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
import yaml

import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

from scipy import interpolate
from numpy.linalg import eig
from copy import deepcopy

np.seterr(divide='ignore')


class Meta(object):
    """Object to handle meta-data input and application."""
    def __init__(self, metastring):
        """Create Meta object from yaml string."""
        self.data = yaml.safe_load(metastring)
        return

    def get(self, key):
        """Lookup data stored under 'key'."""
        return self.data.get(key, '')


class DummyFrag(object):
    def __init__(self, E, t):
        self.E = E
        self.t = t
        return


class EmissionFragment(object):
    """Object to handle operations involving specific emissions fragment.

    EmissionFragment will keep track of the magnitude of this unit of emissions
    over time, which compartment the emissions will enter, as well as optional
    weighting.
    """
    def __init__(self, E, t,
                 compartment=None, weight=None, numboxes=None,
                 weights=None):
        """Create EmissionFragment.
        
        Arguments:
            E: emissions magnitude (1-D array with dimension of t)
            t: corresponding time axis of emissions (1-D array)
            compartment: index of compartment emissions will enter
            weight: simple scaling for these emissions
        """
        assert ( (weights is not None) or
                 ( (compartment is not None)  and
                   (weight is not None) and
                   (numboxes is not None) )
        )
        if weights is not None:
            weights = np.array(weights)
            self.E = E[None,:]*weights[:,None]
        else:
            self.E = np.zeros((numboxes, len(E)))
            self.E[compartment,:] = E*weight

        self.t = t
        return

class ReservoirFragment(object):
    """Object to handle operations involving the reservoirs due to an emissions fragment."""
    def __init__(self,R,t):
        """Create ReservoirFragment object.
        
        Arguments:
            R: reservoir masses through time (array of shape Ncompartments by len(t))
            t: time axis of reservoir masses (1-D array)
        """
        self.R = R
        self.t = t
        return

class FragmentCollection(object):
    """Object to handle collection and combination of fragments.
    
    Collects/creates multiple emissions fragments based on Meta input, 
    can convert emissions fragments to reservoir fragments,
    output combinations of reservoir fragments. 
    """
    def __init__(self, meta):
        """Create FragmentCollection based on Meta object input.
        
        Arguments:
            meta: Meta object describing fragment collection
        """
        self.tags = meta.data['fragments'].keys()
        self.collection = {}
        self.meta = meta
        self.hasrun = False
        self.transfer_matrix = None
        return
    
    def _collect_emission_fragments(self, time_axis):
        """Bring together collection of emissions.

        :param time_axis: (array)
        """
        
        for i,tag in enumerate(self.tags):
            self.collection[tag] = {
                'emis':retrieve_emis_fragment(tag, self.meta, time_axis)
                }
            if i==0:
                Etot = deepcopy(self.collection[tag]['emis'].E)
            else:
                Etot += self.collection[tag]['emis'].E
        self.total = DummyFrag(Etot, time_axis)

    def _collect_reservoir_fragments(self,transfer_matrix):
        """Solve and gather solutions for emissions.

        :param transfer_matrix: (FlowMatrix)
        """
        self.transfer_matrix = transfer_matrix
        for tag in self.tags:
            self.collection[tag]['res'] = solve_fragment_euler(
                self.collection[tag]['emis'], self.meta,
                transfer_matrix )

    def _get_total_reservoir(self,transfer_matrix):
        """Solve sum of emissions rather than collection.

        :param transfer_matrix: (FlowMatrix)
        """
        self.transfer_matrix = transfer_matrix
        tot = solve_fragment_euler(
            self.total, self.meta, transfer_matrix)
        self.collection['total'] = {'res':tot}

    
    def _combine_reservoir_fragments(self):
        """Get total by summing collection. """
        first = True
        for tag in self.tags:
            frag = self.collection[tag]
            sign = frag.get('sign',1.0)
            res = frag.get('res')
            if not first:
                tot.R = tot.R + sign*res.R
            else:
                tot = deepcopy(res)
                tot.R = sign*tot.R
                first = False
        self.collection['total'] = {'res':tot}

    def run_total(self, time_axis, transfer_matrix):
        """Solve total emissions.

        :param time_axis: (array)
        :param transfer_matrix: (FlowMatrix)

        """
        self._collect_emission_fragments(time_axis)
        self._get_total_reservoir(transfer_matrix)
        self.hasrun = True
        
    def run_collection(self, time_axis, transfer_matrix):
        """Calculate reservoirs resulting from emissions.
        
        Arguments:
            time_axis: times onto which reservoirs will be mapped (1-D array)
            transfer_matrix: FlowMatrix object
        """
        self._collect_emission_fragments(time_axis)
        self._collect_reservoir_fragments(transfer_matrix)
        self._combine_reservoir_fragments()
        self.hasrun = True
        return

    def lookup_emisfrag(self, tag):
        """Get EmissionFragment named 'tag'."""
        if not self.hasrun:
            print('Collection has not been run yet...')
        else:
            return self.collection[tag]['emis']

    def lookup_matrix(self,):
        """Get FlowMatrix associated with this collection."""
        if self.transfer_matrix is None:
            print('no matrix defined yet')
        else:
            return self.transfer_matrix

    def lookup_resfrag(self, tag):
        """Get ReservoirFragment named 'tag'."""
        if not self.hasrun:
            print('Collection has not been run yet...')
        else:
            return self.collection[tag]['res']

    def combine_fragments(self, tags='all'):
        """Sum multiple fragments to get total reservoirs.
        Keyword Arguments:
            tags: list of tags to combine (default all)
        """
        if type(tags) in (list,tuple):
            pass
        elif tags in ('all','All','ALL'):
            tags = list(self.collection.keys()) # Add all by default
        resfrag1 = deepcopy(self.lookup_resfrag(tags[0]))
        reservoirs,tout = resfrag1.R, resfrag1.t
        for tag in tags[1:]:
            reservoirs += self.lookup_resfrag(tag).R
        return reservoirs,tout
 
    def get_total(self):
        """Return the totaled reservoir."""
        if not self.hasrun:
            print('Collection has not been run yet...')
        else:
            return self.lookup_resfrag('total')
    
class FragmentMeta(object):
    """Object to handle meta-data specific to a fragment."""
    def __init__(self,dicty,default_file):
        self.data = dicty
        self.default_file = default_file

    def get(self,key):
        """ Return values for fragment metadata here. Return the
        defaults if not specified. """
        if self.data.get(key,None) is None:
            if key == 'interpolation':
                return 'linear'
            elif key == 'compartment':
                return 0
            elif key == 'weight':
                return 1.
            elif key == 'csvheader':
                return 1
            elif key == 'time_column':
                return 'Year'
            elif key == 'filename':
                return self.default_file
            elif key == 'numboxes':
                return 25
            else:
                return None
        else:
            return self.data.get(key)

class FlowMatrix(object):
    """Object to define mass transfer between compartments.
    
    Attributes:
    matrix: matrix representation of mass transfer
    compartments: names of comparments
    
    Methods:
    get_matrix
    interpolate_matrix
    decompose
    get_timescales
    perturbation_analysis
    """
    def __init__(self,matrix,compartment_names,STATIC=True):
        """Create FlowMatrix from matrix, list of compartment names.
        
        Arguments:
            matrix: matrix of relative transfer rates (N x N array)
            comparment_names: list of comparment names associated with matrix (length N)
        
        Keyword Arguments:
            STATIC (Boolean): is the mass transfer static in time? (default True)
        """
        self.matrix = matrix
        self.compartments = compartment_names
        self.numcompartments = len(compartment_names)
        self.decompose()
        self.STATIC = STATIC
        return

    def get_matrix(self,t=None):
        """Get matrix of mass transfer at given time.
        
        If time specified, matrix will be time-interpolate to correspond. If no
        time is given, default static matrix will be used.
        
        Keyword Argument:
            t: time at which to calculate matrix (default None)
        """
        if self.STATIC:
            return self.matrix
        elif (t is None):
            return self.matrix
        else:
            return self._get_matrix_for_time(t)
    
    def _get_matrix_for_time(self,t):
        # self.Minterp is trained interpolator
        return self.Minterp(t)
    
    def interpolate_matrix(self,matrix_of_t,t_matrix):
        """Make matrix interpolation to requested time possible."""
        # matrix_of_t must be N x N x numt array
        self.Minterp = interpolate.interp1d(t_matrix,matrix_of_t)
        return
    
    def decompose(self):
        """ Decompose Matrix into eigenvalues/vectors. 
        
        Returns:
            eigenvalues (array): unordered array of eigenvalues
            eigenvectors (array): normalized right eigenvectors 
                   ( eigenvectors[:,i] corresponds to eigenvalues[i] ) 
        """
        eigenvalues,eigenvectors = eig(self.matrix)
        self.eigenvalues = eigenvalues
        self.eigenvectors = eigenvectors
        self.timescales = -1./np.float64(np.real(eigenvalues))
        self.residence_times = -1/np.float64(np.diagonal(self.matrix))
        return eigenvalues, eigenvectors

    def get_timescales(self):
        """Get timescales associated with matrix eigenvalues."""
        return self.timescales

    def perturbation_analysis(self, time_horizon, log=True,
                              compartment_index=0, numtimes=1000,
                              mintime=1e-3):
        """ Do a perturbation analysis over given time_horizon.
        
        Arguments:
            time_horizon (Scalar): number of years afterwhich to end
            log (Boolean): use a log time axis? (default=True)
            compartment_index (Int): compartment of original perturbation
                     (default=0)
            numtimes (Int): number of time points to calculate
                     (default=1000)
            mintime (float): number of years after perturbation for first
                    time point (default=1e-3)

        Returns:
            perturbation_output (array): compartment fractions of original perturbation.
                                    (numcompartments x numtimes)
            perturbation_times (array): time points (numtimes x 1)
        """

        perturbation_output = np.zeros((self.numcompartments,numtimes),dtype=complex)
        if log:
            times = np.logspace(np.log10(mintime), np.log10(time_horizon),
                               numtimes)
        else:
            times = np.linspace(mintime, time_horizon,
                               numtimes)
    
        initcond = np.zeros(self.numcompartments)
        initcond[compartment_index] = 1. # initial conditions
        inverse_eigenvectors = np.linalg.inv(self.eigenvectors)
        eiginitcond = np.dot(inverse_eigenvectors, initcond)
        for i,t in enumerate(times):
            in_eig_space = np.dot(  np.diag( np.exp(self.eigenvalues*t) ),
                                    eiginitcond  )
            perturbation_output[:,i] = np.dot(self.eigenvectors, 
                                              in_eig_space)
        return np.real(perturbation_output), times    


class Rate(object):
    """Stores rate and its metadata to be applied in model matrix."""

    def __init__(self, name, k=0., compartment_to=None, compartment_from=None,
                 compartment_to_ind=None, compartment_from_ind=None,
                 notes=" ", ref=" ", k_of_t=None, delta_k_func = None):
        self.k = k # magnitude
        self.compartment_to = compartment_to # flow to this compartment
        self.compartment_from = compartment_from # flow from this comp.
        self.compartment_to_ind = compartment_to_ind # index of comp. to
        self.compartment_from_ind = compartment_from_ind # index of comp. from
        self.notes = notes # e.g. text description for humans
        self.ref = ref # reference to source 
        self.delta_k_func = delta_k_func # function to scale for uncertainty
        self.dk = 1. # value to scale magnitude
        
        if self.compartment_to in ['','None',None]:
            self.TRANSFER = False
            self.LOSS = True # A loss is a transfer out of system
        else:
            self.TRANSFER = True
            self.LOSS = False
        self.name = name
        if k_of_t is None:
            self.k_of_t = lambda x:k # so constant k can be called with time
        else:
            self.k_of_t = k_of_t
            self.k = k_of_t(0.) # should be needless if using .get_k()
        return

    def assign_compartment_indices(self, compartments):
        """ Assigns indices to "to" and "from" compartments for matrix building."""
        if self.TRANSFER:
            if self.compartment_to_ind is None:
                self.compartment_to_ind = compartments.index(self.compartment_to)
        if self.compartment_from_ind is None:
            self.compartment_from_ind = compartments.index(self.compartment_from)
        return

    def get_k(self,t=None):
        """ Return current rate magnitude (at time t if given)."""
        dk = self.dk
        if t is None:
            return self.k * dk
        else:
            return self.k_of_t(t) * dk
    
    def set_dk(self,mode='Deterministic',whitelist=None):
        """ Set the factor by which to perturb rate magnitude.
        When mode=='Deterministic' this makes no change (dk=1)
        mode=='MonteCarlo' will sample the scaling factor from 
        this rate's uncertainty distribution.
        mode=='Sensitivity' will increase the rate by 10%
        """
        WHITE = False
        if ( (whitelist is not None) and (self.name in whitelist) ):
            WHITE = True

        if mode in ['Deterministic']:
            dk = 1.
            
        if mode in ['MonteCarlo']:
            if (self.delta_k_func is None) or (not WHITE):
                dk = 1.
            else:
                dk = self.delta_k_func(1)

        if mode in ['Sensitivity']:
            if WHITE:
                dk = 1.1
            else:
                dk = 1.
            
        if mode not in ['Deterministic','MonteCarlo','Sensitivity']:
            print("Unplanned condition in Rate.set_dk(): %s"%mode)

        self.dk = dk
        

    def get_notes(self):
        return self.notes

class Model(object):
    """Handles set-up with transfers and losses, applies solving, etc."""
    def __init__(self, compartments, mode='Deterministic', seed=1):
        self.rates = [] # list of Rate objects
        self.rates_dict = {} # dictionary with lookup names for Rate objects
        self.compartments = compartments # names of compartments
        ncomp = len(compartments)
        self.ncomp = ncomp
        self.transfer_matrix_raw = np.zeros((ncomp,ncomp)) # numerical matrix
        self.transfer_matrix = None # will be filled with FlowMatrix obj.
        self.loss_matrix_raw = np.zeros((ncomp,ncomp))
        self.loss_matrix = None
        self.initial_conditions = np.zeros(ncomp)
        self.matrix = None
        self.time_axis = None # on which solution will be given
        self.fragment_collection = None # FragmentCollection obj.
        self.time_units = 'time units'
        self.mass_units = 'mass units'
        self.mode = mode # Deterministic, MonteCarlo, Sensitivity
        self.seed = seed # numpy random seed for consistent (random) results
        self.dk_ratelist = None # whitelist for rates to alter in MC or Sens mode
        self.aux_outputs = {} # side outputs of interest such as implied POC flux.
        self.maximum_timestep = None
        return

    def set_mode(self,mode):
        """Set the mode of model run (Deterministic, MC, Sens.)"""
        self.mode = mode
        
    def set_seed(self,seed):
        """Set the numpy random seed for model run."""
        self.seed = seed

    def set_dk_ratelist(self,ratenames):
        """Set whitelist for rates to perturb."""
        self.dk_ratelist = ratenames

    def add_rate(self,rate):
        """Add new Rate object to registry for forming matrix."""
        assert type(rate) == type(Rate('')), "Rate object required"
        assert (rate.name not in self.rates_dict), "Rate name already in use"
        rate.assign_compartment_indices(self.compartments)
        self.rates.append(rate)
        self.rates_dict[rate.name] = rate
        return

    def gen_matrices(self, t=0):
        """Generate transfer, loss matrices using registered Rate objs."""
        self.loss_matrix_raw *= 0.
        self.transfer_matrix_raw *= 0.
        
        for rate in self.rates:
            if rate.LOSS:
                i = rate.compartment_from_ind
                self.loss_matrix_raw[i,i] -= rate.get_k(t)
            elif rate.TRANSFER:
                i = rate.compartment_from_ind
                j = rate.compartment_to_ind
                self.transfer_matrix_raw[i,i] -= rate.get_k(t)
                self.transfer_matrix_raw[j,i] += rate.get_k(t)
            else:
                raise "Rate must be transfer or loss"

        return

    def get_matrix(self,):
        """Get matrix obj."""
        return self.matrix

    def closest_time_ind(self, t):
        """Get the closest time on time_axis to given t."""
        return np.argmin(np.abs(t-self.time_axis))
    
    def get_flux(self, name, t=0, reservoirs=None):
        """Get the flux associated with a given rate (at time t if given)."""
        rate = self.rates_dict[name]
        if reservoirs is not None:
            x = reservoirs[rate.compartment_from_ind]
        else:
            t_ind = self.closest_time_ind(t)
            x = self.reservoirs[rate.compartment_from_ind,t_ind]

        return x*rate.get_k(t)

    def build(self,time_axis=None):
        """Build the model (generate matrix) given the parameters current set."""
        np.random.seed(self.seed)
        for rate in self.rates:
            rate.set_dk(mode=self.mode,whitelist=self.dk_ratelist)

        if time_axis is None:
            self.gen_matrices()
            self.matrix = FlowMatrix(self.transfer_matrix_raw + self.loss_matrix_raw,
                                     compartment_names=self.compartments,
                                     STATIC=True)
            self.transfer_matrix = FlowMatrix(self.transfer_matrix_raw,
                                              compartment_names=self.compartments,
                                              STATIC=True)
            self.loss_matrix = FlowMatrix(self.loss_matrix_raw,
                                          compartment_names=self.compartments,
                                          STATIC=True)
        else:
            flows = np.zeros((self.ncomp,self.ncomp,len(time_axis)))
            losss = np.zeros((self.ncomp,self.ncomp,len(time_axis)))
            tots  = np.zeros((self.ncomp,self.ncomp,len(time_axis)))
            for i,t in enumerate(time_axis):
                self.gen_matrices(t)
                flows[:,:,i] = self.transfer_matrix_raw
                losss[:,:,i] = self.loss_matrix_raw
                tots[:,:,i] = flows[:,:,i] + losss[:,:,i]
            self.matrix = FlowMatrix(self.transfer_matrix_raw + self.loss_matrix_raw,
                                     compartment_names=self.compartments,
                                     STATIC=False)
            self.transfer_matrix = FlowMatrix(self.transfer_matrix_raw,
                                              compartment_names=self.compartments,
                                              STATIC=False)
            self.loss_matrix = FlowMatrix(self.loss_matrix_raw,
                                          compartment_names=self.compartments,
                                          STATIC=False)
            self.matrix.interpolate_matrix(tots,time_axis)
            self.transfer_matrix.interpolate_matrix(flows,time_axis)
            self.loss_matrix.interpolate_matrix(losss,time_axis)
            
                
        self.residence_times = self.matrix.residence_times
        self.set_auto_timestep()
        self.removal_timescales = self.loss_matrix.residence_times
        return

    def set_auto_timestep(self):
        min_time = np.min(self.residence_times)
        self.maximum_timestep = 0.5*min_time
    
    def get_residence_times(self,t=0):
        """Get residence times in each compartment."""
        return self.residence_times

    def get_removal_timescales(self,t=0):
        """Get removal timescales (from the system) in each compartment."""
        return self.removal_timescales
    
    def get_steady_state(self, E):
        """ Calculate steady-state reservoirs for given emissions array.

        Arguments:
        matrix (array): (NxN) rates
        E (array): (Nx1) emissions

        Returns:
        (Nx1) array of reservoir masses
        """
        return np.linalg.solve(self.matrix.matrix,-E)
    
    def set_time_axis(self, time):
        """Set the time axis on which solution will be calculated."""
        if len(time) == 2:
            start_time = time[0]
            end_time = time[1]
            dt = self.maximum_timestep
            #if dt < 0:
            #    print(self.residence_times)
            steps = int((end_time-start_time)/dt)+1
            self.time_axis = np.linspace(start_time, end_time, steps)
        else:
            self.time_axis = time
        return self.time_axis

    def setup_emissions(self, meta):
        """Generate emissions based on given meta object."""
        self.emissions_collection = FragmentCollection(meta)
        return

    def run(self, FAST_MODE=False):
        """Run the model using currently defined rates, emissions."""
        if FAST_MODE:
            self.emissions_collection.run_total(self.time_axis,self.matrix)
#            self.reservoirs = self.emissions_collection.run_total(self.time_axis,self.matrix) \
#                              + self.initial_conditions[:,None]
        else:
            self.emissions_collection.run_collection(self.time_axis,self.matrix)
        self.reservoirs = self.emissions_collection.get_total().R \
                          + self.initial_conditions[:,None]
        return

    def get_reservoir(self,name,t):
        """Get given reservoir magnitude at given time."""
        compi = self.compartments.index(name)
        ti = np.argmin(np.abs(t-self.time_axis))
        return self.reservoirs[compi,ti]

    def set_initial_conditions(self,IC):
        """Set initial conditions for model run."""
        self.initial_conditions = IC
        return

    def perturbation_plot(self,tmax,rev=False,colors=None,**kwargs):
        """Plot perturbation diagram to time horizon tmax.

        :param tmax: (float) maximum time to show
        :param rev: (Boolean) whether to flip plot
        :param colors: (List) colors for each compartment

        """
        if rev:
            inc = 1
        else:
            inc = -1

        sizex,sizey = max(self.ncomp,6),max(self.ncomp/2,3)
        plt.figure(figsize=(sizex,sizey))
        pert,tpert = self.matrix.perturbation_analysis(tmax,**kwargs)
        if colors:
            plt.stackplot(tpert,pert[::inc,:],
                          labels=self.compartments[::inc],edgecolor='k',
                          colors=colors)
        else:
            plt.stackplot(tpert,pert[::inc,:],
                          labels=self.compartments[::inc],edgecolor='k')
            
        plt.semilogx()
        plt.legend(loc='upper right')
        plt.xlabel('Time (%s)'%self.time_units,fontsize=15)
        plt.ylabel('Fraction of perturbation',fontsize=15)
        return pert,tpert
        
    def timescale_analysis(self):
        """Print residence and removal timescales for each compartment. """
        res = self.get_residence_times()
        rem = self.get_removal_timescales()
        print('Residence times:')
        for i,comp in enumerate(self.compartments):
            print(f'{comp} {res[i]:.1f} {self.time_units}')
        print('\nRemoval timescales:')
        for i,comp in enumerate(self.compartments):
            print(f'{comp} {np.abs(rem[i]):.1f} {self.time_units}')
        return

    def eigen_analysis_plot(self,no_losses=False):
        """Plot eigenvalues and vectors for current model setup.

        :param no_losses: (Boolean)

        """
        if no_losses:
            vecs = self.transfer_matrix.eigenvectors
            vals = self.transfer_matrix.eigenvalues
        else:
            vecs = self.matrix.eigenvectors
            vals = self.matrix.eigenvalues

        plt.figure(figsize=(self.ncomp*2,self.ncomp))
        inds = range(self.ncomp)
        for i in inds:
            vec = 0.5*vecs[:,i]
            plt.plot(i+np.real(vec),inds,'o-',color='k')
            plt.plot(i+np.imag(vec),inds,'o-',color='g')
            plt.axvline(i,color='gray',linestyle='--')
            plt.axhline(i,color='gray',linestyle='--')
        printvals = []
        for x in -1/vals:
            if x.imag==0.:
                printvals.append('%.1f'%x.real)
            else:
                printvals.append('%.1f%+.1fi'%(x.real,x.imag))
        plt.xticks(inds,printvals,fontsize=15)
        plt.yticks(inds,self.compartments,fontsize=15)
        return
    
    def set_units(self,time=None,mass=None):
        """Set time and/or mass units for model.

        :param time: (str), optional; time units
        :param mass: (str), optional; mass units

        """
        if time is not None:
            self.time_units = time
        if mass is not None:
            self.mass_units = mass
        return
    
def get_meta(metastring):
    """ Collect the metadata needed for model run. """
    meta = Meta(metastring)
    return meta

def retrieve_emis_fragment(tag, meta, t_out):
    """ Use meta object to look up emis fragment for given tag. """
    fragment_meta = FragmentMeta(meta.data['fragments'][tag],
                                 meta.get('emissions_filename'))
    emissions_table = pd.read_csv(fragment_meta.get('filename'),
                                  header=fragment_meta.get('csvheader'))
    emisfrag = get_emis_fragment(emissions_table,
                                 fragment_meta.get('column_name'),
                                 fragment_meta.get('time_column'),
                                 t_out,
                                 fragment_meta.get('interpolation'),
                                 fragment_meta.get('compartment'),
                                 fragment_meta.get('weight'),
                                 fragment_meta.get('numboxes'),
                                 fragment_meta.get('weights'))
    return emisfrag

def get_emis_fragment(emissions, fragment_id, time_id, t_out,
                      interpolation, compartment, weight, numboxes,
                      weights):
    """ Find and load the given emissions data over the requested time period.
  
    EmisFragment object should include time axis
    """
    
    e_in, t_in = emissions[fragment_id], emissions[time_id]
    finterp = interpolate.interp1d(t_in, e_in, kind=interpolation)
    e_out = finterp(t_out)
    emisfrag = EmissionFragment(e_out, t_out,
                                compartment=compartment,
                                weight=weight,
                                numboxes=numboxes,
                                weights=weights)
    return emisfrag

def solve_fragment_euler(emissionfragment,meta,transfer_matrix):
    """ Solve EmissionFragment to get Reservoir on output time axis. """

    matrix = transfer_matrix.get_matrix()
    E = emissionfragment.E
    R = np.zeros_like(E)
    t = emissionfragment.t
    for i in range(len(t)-1):
        dRdt = np.dot(transfer_matrix.get_matrix(t[i]),R[:,i])+E[:,i]
        dt = t[i+1]-t[i]
        dR = dRdt*dt
        R[:,i+1] = R[:,i] + dR
    reservoirfragment = ReservoirFragment(R,t)
    return reservoirfragment

def assemble_fragments(meta):
    """ Run EmissionFragments and assemble them and their resulting
    ReservoirFragments by addition or substraction 
    to get model results. """
    return fragmentcollection

def plot_results(meta,fragmentcollection):
    """ Make requested plots of model output. """
    return

