import requests
import os.path

print('STARTING SCRIPT: 0-Fetch_Dataverse_Files.py')

# dictionary containing emission filenames (keys) and corresponding Dataverse IDs (values)
file_dict = {
    'Hg_eff_2005.nc':7127988,
    'Hg_eff_2006.nc':7127979,
    'Hg_eff_2007.nc':7127996,
    'Hg_eff_2008.nc':7128012,
    'Hg_eff_2009.nc':7127975,
    'Hg_eff_2010.nc':7127985,
    'Hg_eff_2011.nc':7128009,
    'Hg_eff_2012.nc':7128002,
    'Hg_eff_2013.nc':7127984,
    'Hg_eff_2014.nc':7128004,
    'Hg_eff_2015.nc':7127989,
    'Hg_eff_2016.nc':7128021,
    'Hg_eff_2017.nc':7128008,
    'Hg_eff_2018.nc':7128013,
    'Hg_eff_2019.nc':7128017,
    'Hg_eff_2020.nc':7128007,
    'Hg_eff_2021.nc':7127991,

    'Hg_exp_2005.nc':7128018,
    'Hg_exp_2006.nc':7128000,
    'Hg_exp_2007.nc':7127995,
    'Hg_exp_2008.nc':7127973,
    'Hg_exp_2009.nc':7128010,
    'Hg_exp_2010.nc':7127990,
    'Hg_exp_2011.nc':7128011,
    'Hg_exp_2012.nc':7127993,
    'Hg_exp_2013.nc':7127972,
    'Hg_exp_2014.nc':7128003,
    'Hg_exp_2015.nc':7127989,
    'Hg_exp_2016.nc':7128016,
    'Hg_exp_2017.nc':7127974,
    'Hg_exp_2018.nc':7128020,
    'Hg_exp_2019.nc':7128001,
    'Hg_exp_2020.nc':7127987,
    'Hg_exp_2021.nc':7127976,

    'Hg_pvf_2005_2D.nc':7127998,
    'Hg_pvf_2006_2D.nc':7127992,
    'Hg_pvf_2007_2D.nc':7127981,
    'Hg_pvf_2008_2D.nc':7127977,
    'Hg_pvf_2009_2D.nc':7128014,
    'Hg_pvf_2010_2D.nc':7127971,
    'Hg_pvf_2011_2D.nc':7127982,
    'Hg_pvf_2012_2D.nc':7128006,
    'Hg_pvf_2013_2D.nc':7127986,
    'Hg_pvf_2014_2D.nc':7128015,
    'Hg_pvf_2015_2D.nc':7127978,
    'Hg_pvf_2016_2D.nc':7127997,
    'Hg_pvf_2017_2D.nc':7127980,
    'Hg_pvf_2018_2D.nc':7127983,
    'Hg_pvf_2019_2D.nc':7128005,
    'Hg_pvf_2020_2D.nc':7128019,
    'Hg_pvf_2021_2D.nc':7127999,
    }

# Download model output first
out_path = f'../Output/Publication_Data/Model_Data.nc'
# delete file if it exists:
check_file = os.path.exists(out_path)
if check_file == True:
    os.remove(out_path)
# download file
r = requests.get(f'https://dataverse.harvard.edu/api/access/datafile/7346480', stream=True)
with open(out_path, 'wb') as outfile:
    outfile.write(r.content)

# now loop through and download emission files
for file, id in file_dict.items():
    out_path = f'../Output/Publication_Data/Emission_Files/{file}'
    # delete file if it exists:
    check_file = os.path.exists(out_path)
    if check_file == True:
        os.remove(out_path)
    # download file
    r = requests.get(f'https://dataverse.harvard.edu/api/access/datafile/{id}', stream=True)
    with open(out_path, 'wb') as outfile:
        outfile.write(r.content)

print('COMPLETED SCRIPT: 0-Fetch_Dataverse_Files.py')