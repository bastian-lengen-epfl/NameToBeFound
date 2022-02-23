"""
This module contains a function that fits the distance ladder to the different dataframes given.
"""
import numpy as np
import pandas as pd
import Fit_parameters as Fp

def fit_distance_ladder(DF_dict, work_dir='./'):
    '''
    Return the raw result of the fit parameters

    :type   DF_dict: dictionary of pandas DataFrame
    :param  DF_dict: Dictionary that contains the DataFrame that will be fitted.
    :type   work_dir: string
    :param  work_dir: Working directory that contains the Data/ folder.
    '''
    Data_dir = work_dir + 'Data/'

    # We first create an empty signal vector y, design matrix L and parameters vector q and Sigma
    [y, L, q, Sigma] = 4*[np.empty(0)]
    # For the Cepheids
    if Fp.include_Cepheids == True:
        ### Load the DF
        # Use the Cepheids DF and check for its consistency with the Fit_parameters.py file
        Cepheids = DF_dict['Cepheids']
        galaxies_Cep = Cepheids['Gal'].drop_duplicates().reset_index(drop=True) # List of Cepheids-host galaxy
        if len(galaxies_Cep)!=Fp.N_galaxies_Cep:
            print('ERROR: make sure that the value for N_galaxies_Cep in the Fit_parameters.py corresponds to \
                   the number of the galaxy in the Cepheids.csv file!')
            return
        # The anchors Cepheids
        Cepheids_anchors = DF_dict['Cepheids_anchors']
        anchors_Cep =  Cepheids_anchors['Gal'].drop_duplicates().reset_index(drop=True) # List of Cepheids-host galaxy
        if len(anchors_Cep)!=Fp.N_anchors_Cep:
            print('ERROR: make sure that the value for N_anchors_Cep in the Fit_parameters.py corresponds to \
                   the number of the galaxy in the Cepheids_anchors.csv file!')
            return
        # The MW Cepheids
        if Fp.include_MW == True:
            Cepheids_MW = DF_dict['Cepheids_MW']

        ### fill y
        y = np.append(y, Cepheids['mW'])
        y = np.append(y, Cepheids_anchors['mW'])
        if Fp.include_MW == True:
            if Fp.fixed_zp == True:
                y = np.append(y, Cepheids_MW['mW'] \
                              - 10 + 5 * np.log10(Cepheids_MW['pi'])\
                              + 5/np.log(10)*Fp.zp/Cepheids_MW['pi'])
            else:
                y = np.append(y, Cepheids_MW['mW'] \
                              - 10 + 5 * np.log10(Cepheids_MW['pi']))
        if Fp.fixed_Zw == True:
            y = y-Fp.Zw*Cepheids['M/H']

        ### Create q (q_string for future dict so it's easier to access)
        # The parameters will be in that ordre:
        # mu, MW, b, (bs), (zw), (zp)
        q_string = []
        for gal in galaxies_Cep:
            q_string.append(f'mu_{gal}')
        if Fp.PLR_break == True:
            q_string.append('b_s')
            q_string.append('b_l')
        else:
            q_string.append('b')
        if Fp.fixed_Zw == False:
            q_string.append('Zw')
        if ((Fp.include_MW == True) and (Fp.fixed_zp==False)):
            q_string.append('zp')
        q = np.zeros(len(q_string))


        ### fill L
        # First need to check and split the period for a PL break

    # Solve

    return y