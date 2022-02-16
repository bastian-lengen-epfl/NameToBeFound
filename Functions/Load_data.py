"""
This module contains functions that load the previously pre-processed data and return the correspondig dataset
that is used by the other functions.
"""
import pandas as pd
from Fit_parameters import include_Cepheids, include_MW, include_TRGB, fit_aB

def load_data(work_dir='./'):
    Data_dir = work_dir + 'Data/'

    # Create an empty dictionnary of dataset
    DF_dict={}

    # Loads the Cepheids related DataFrame
    if include_Cepheids==True:
        print('Loading Cepheids & Cepheids anchors...')
        DF_dict['Cepheids'] = pd.read_csv(Data_dir+'Cepheids.csv', sep=',')
        DF_dict['Cepheids_anchors'] = pd.read_csv(Data_dir + 'Cepheids_anchors.csv', sep=',')
        if include_MW==True:
            print('Loading MW Cepheids...')
            DF_dict['Cepheids_MW'] = pd.read_csv(Data_dir+'Cepheids_MW.csv', sep=',')
        print('Loading the SNe for the Cepheids')
        DF_dict['SNe_Cepheids'] = pd.read_csv(Data_dir + 'SNe_Cepheids.csv', sep=',')

    # Loads the TRGB related DataFrame
    if include_TRGB==True:
        print('Loading TRGB & TRGB anchors...')
        DF_dict['TRGB'] = pd.read_csv(Data_dir+'TRGB.csv', sep=',')
        DF_dict['TRGB_anchors'] = pd.read_csv(Data_dir + 'TRGB.csv', sep=',')
        print('Loading the SNe for the TRGB')
        DF_dict['SNe_TRGB'] = pd.read_csv(Data_dir + 'SNe_Cepheids.csv', sep=',')

    # Loads the SNe for the aB fit
    if fit_aB==False:
        print('Loading SNe for the aB fit...')
        DF_dict['SNe_Hubble'] = pd.read_csv(Data_dir + 'SNe_Hubble.csv', sep=',')

    return DF_dict