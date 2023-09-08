import pandas as pd
import scipy.stats
import numpy as np

rename_dict ={
   'NAME': 'snp',
   # NAME: 'pos_name',
   'CHR': 'chromosome',
   'POS': 'position',
   'REF_ALLELE' : 'reference_allele',
   'EFFECT_ALLELE' : 'effect_allele',
   'BETA' : 'beta',
   'SE' : 'se',
    #None: 'z',
    #None: 'pval',
   'MAF': 'eaf',
   'N_OBS': 'n_iids'
}

necessary_colnames = [
    # 'snp',
    'pos_name',
    'chromosome',
    'position',
    'effect_allele',
    'reference_allele',
    'beta',
    'se',
    'z',
    'pval',
    # 'eaf',
    'n_iids', ]

def convert_to_mr_link2_format(filename):
    this_df = pd.read_csv(filename, sep='\t')

    this_df = this_df.rename(rename_dict, axis=1)
    this_df['pos_name'] = this_df.snp
    this_df['z'] = this_df.beta / this_df.se
    this_df['pval'] = 2*scipy.stats.norm.sf(np.abs(this_df.z))
    this_df['beta'] = this_df.z / np.sqrt( this_df.n_iids  + this_df.z**2)
    this_df['se'] = 1 / np.sqrt(this_df.n_iids + this_df.z ** 2)
    this_df = this_df[necessary_colnames]
    return this_df


if __name__ == '__main__':
    prepend = 'simulated_sumstats/mr_link2_yes_causalex1-n-100_ex2-n-100_ex1-b-0.0_ex2-b-0.4_overl-0_dir_pleio-1_run_0'
    exposure_df = convert_to_mr_link2_format(prepend + '_exposure_sumstats.txt')
    exposure_df.to_csv('testing_files/non_causal_exposure.txt', sep='\t', index=False)
    outcome_df  = convert_to_mr_link2_format(prepend + '_outcome_sumstats.txt')
    outcome_df.to_csv('testing_files/non_causal_outcome.txt', sep='\t', index=False)

    prepend = 'simulated_sumstats/mr_link2_yes_causalex1-n-100_ex2-n-100_ex1-b-0.4_ex2-b-0.4_overl-0_dir_pleio-1_run_0'
    exposure_df = convert_to_mr_link2_format(prepend + '_exposure_sumstats.txt')
    exposure_df.to_csv('testing_files/yes_causal_exposure.txt', sep='\t',index=False)
    outcome_df = convert_to_mr_link2_format(prepend + '_outcome_sumstats.txt')
    outcome_df.to_csv('testing_files/yes_causal_outcome.txt', sep='\t', index=False)
