"""
This module contains functions that load the previously pre-processed data and return the corresponding dataset
that is used by the other functions.
"""
import pandas as pd
from Fit_parameters import include_Cepheids, include_MW, include_TRGB, fit_aB

def load_data(work_dir='./'):
    '''
    :param work_dir:
    :return:
    '''
    Data_dir = work_dir + 'Data/'

    # Create an empty dictionnary of dataset
    DF_dict={}

    # Loads the Cepheids related DataFrame
    if include_Cepheids==True:
        print('Loading Cepheids & Cepheids anchors...')
        DF_dict['Cepheids'] = pd.read_csv(Data_dir+'Cepheids.csv', sep=',')
        DF_dict['Cepheids'].columns = ['Gal', 'logP', 'mW', 'sig_mW', '[O/H]', 'z'] # Give a consistent name to the columns
        DF_dict['Cepheids_anchors'] = pd.read_csv(Data_dir + 'Cepheids_anchors.csv', sep=',')
        DF_dict['Cepheids_anchors'].columns = ['Gal', 'logP', 'mW', 'sig_mW', '[O/H]', 'z', 'mu', 'sig_mu']  # Give a consistent name to the columns
        if include_MW==True:
            print('Loading MW Cepheids...')
            DF_dict['Cepheids_MW'] = pd.read_csv(Data_dir+'Cepheids_MW.csv', sep=',')
            DF_dict['Cepheids_MW'].columns = ['Gal', 'logP', 'mW', 'sig_mW', '[O/H]', 'z', 'pi', 'sig_pi']  # Give a consistent name to the columns
        print('Loading SNe for the Cepheids-host galaxies...')
        DF_dict['SNe_Cepheids'] = pd.read_csv(Data_dir + 'SNe_Cepheids.csv', sep=',')
        DF_dict['SNe_Cepheids'].columns = ['Gal', 'mB', 'sig_mB'] # Give a consistent name to the columns

    # Loads the TRGB related DataFrame
    if include_TRGB==True:
        print('Loading TRGB & TRGB anchors...')
        DF_dict['TRGB'] = pd.read_csv(Data_dir+'TRGB.csv', sep=',')
        DF_dict['TRGB'].columns = ['Gal', 'm', 'sig_m', 'A', 'V-I', 'z'] # Give a consistent name to the columns
        DF_dict['TRGB_anchors'] = pd.read_csv(Data_dir + 'TRGB_anchors.csv', sep=',')
        DF_dict['TRGB_anchors'].columns = ['Gal', 'm', 'sig_m', 'A', 'V-I', 'z', 'mu', 'sig_mu'] # Give a consistent name to the columns
        print('Loading SNe for the TRGB-host galaxies...')
        DF_dict['SNe_TRGB'] = pd.read_csv(Data_dir + 'SNe_TRGB.csv', sep=',')
        DF_dict['SNe_TRGB'].columns = ['Gal', 'mB', 'sig_mB'] # Give a consistent name to the columns

    # Loads the SNe for the aB fit
    if fit_aB==True:
        print('Loading SNe for the aB fit...')
        DF_dict['SNe_Hubble'] = pd.read_csv(Data_dir + 'SNe_Hubble.csv', sep=',')
        DF_dict['SNe_Hubble'].columns = ['name', 'mB', 'sig_mB', 'z', 'sig_z'] # Give a consistent name to the columns

    return DF_dict