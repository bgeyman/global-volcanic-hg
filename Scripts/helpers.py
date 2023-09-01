import pandas as pd

def subset_df(df, kv_dict:dict):
    ''' Subsets a dataframe based on an arbitrary number of column condition matches.
        Example:
            select_df(df, kv_dict={'Type':'pvf', 'year':2008})
        Usage Note: 
            for now, equals ("==") is the only conditional supported. '''
    
    for key in kv_dict.keys():
        df = df[df[key]==kv_dict[key]]
    
    return df

def make_rafaj_custom_cmap():
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from Additional_Colormaps import make_hex_list, make_custom_cmap_from_hex_list, load_all_cmaps

    new_cmaps = load_all_cmaps()
    n_colors = [4, 4, 6]
    total_n = np.sum(np.array(n_colors))

    rgb_values = (np.array(new_cmaps['rafaj_AQ'].colors[0:2])[:,:3]*255).astype(int)
    test_hex_list = make_hex_list(rgb_values)
    colors1       = make_custom_cmap_from_hex_list(hex_list=test_hex_list, n_colors=n_colors[0], cmap_name='c1')
    colors1       = colors1(np.linspace(0, 1, n_colors[0]))

    rgb_values    = (np.array(new_cmaps['rafaj_AQ'].colors[2:5])[:,:3]*255).astype(int)
    test_hex_list = make_hex_list(rgb_values)
    colors2       = make_custom_cmap_from_hex_list(hex_list=test_hex_list, n_colors=n_colors[1], cmap_name='c2')
    colors2       = colors2(np.linspace(0, 1, n_colors[1]))

    rgb_values    = (np.array(new_cmaps['rafaj_AQ'].colors[5:])[:,:3]*255).astype(int)
    test_hex_list = make_hex_list(rgb_values)
    colors3       = make_custom_cmap_from_hex_list(hex_list=test_hex_list, n_colors=n_colors[2], cmap_name='c3')
    colors3       = colors3(np.linspace(0, 1, n_colors[2]))

    all_colors = np.vstack((colors1, colors2, colors3))
    newcmap = ListedColormap(all_colors)

    return newcmap
