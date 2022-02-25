import sys
import Fit_parameters as Fp
from Fitting import fit_distance_ladder
from Load_data import load_data
from Relativistic_corrections import RLB_correction, K_corr_Cep


def main() -> int:
    # Load the data
    DF_dict = load_data()

    # Cepheids relativistic corrections
    if Fp.include_Cepheids == True:
        if Fp.RLB_correction == True:
            print('Correcting Cepheids for the RLB...')
            DF_dict = RLB_correction(DF_dict)
        if Fp.Kcorr_Cep == True:
            print('K-correcting the Cepheids...')
            DF_dict = K_corr_Cep(DF_dict)
    # TRGB relativistic corrections
    if Fp.include_TRGB == True:
        if Fp.Kcorr_TRGB == True:
            pass

    # Fit
    y, q_dict, L, chi2_reduced = fit_distance_ladder(DF_dict)
    # Outliers rejection (Cepheids and SNe only)
    if Fp.outlier_rejection == True:
        pass
    print(10**(q_dict['5logH0']/5))
    print(chi2_reduced)
    return 0

if __name__ == '__main__':
    sys.exit(main())  # next section explains the use of sys.exit