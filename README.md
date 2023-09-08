# MR-link 2

MR-link 2 is a _cis_ MR method that is pleiotropy robust. 
Inference with MR-link 2 requires **pre-harmonized** summary statistics of an exposure and an outcome, and a 
genotype reference file. 

Please keep an eye on this space for our paper where we validate MR-link 2 on an extensive set of validation datasets.
The main benefit of MR-link 2 is that it has lower false positive rates than other _cis_ MR methods, while retaining
similar or better discriminative ability.

Disclaimer: MR-link 2 is still under development, so it can be subject to change.

### Requirements
MR-link 2 has been tested on macOS X and Linux combined with Python 3.9, 3.10 and 3.11. 
Although not tested, every Python version from 3.6 onwards should work.    
We require some (standard) python packages to be installed, these are: `numpy`, `scipy` and `pandas`.
If they haven't been installed, please install these using pip.
In the command line (shell, terminal), type: 

```{bash}
pip3 install numpy scipy pandas
```
On top of this, we require plink1.9 to be present in your PATH variable. 
Check this by typing `which plink` in your  shell. 
If this prints a path, you have plink installed in your path.


### Example
If you want to test MR-link 2 we have two examples:

This command tests for a causal effect in a region of synthetic data  
```{bash}
python3 mr_link_2_standalone.py \
            --reference_bed example_files/reference_cohort \
            --sumstats_exposure example_files/yes_causal_exposure.txt \
            --sumstats_outcome example_files/yes_causal_outcome.txt \
            --out example_of_a_causal_effect.txt
```

This command tests for a non-causal effect in a region of synthetic data:
```{bash}
python3 mr_link_2_standalone.py \
            --reference_bed example_files/reference_cohort \
            --sumstats_exposure example_files/non_causal_exposure.txt \
            --sumstats_outcome example_files/non_causal_outcome.txt \
            --out example_of_a_non_causal_effect.txt
```

After running these two commands (takes about 2 seconds each), they will output two tab separated files with results: 
`example_of_a_causal_effect.txt` and `example_of_a_non_causal_effect.txt`.
```
# causal effect
region                  var_explained   alpha                   se(alpha)               p(alpha)                sigma_y                 se(sigma_y)             p(sigma_y)              sigma_x                 function_time
2:101532661-103480976   0.99            0.5283473895075494      0.05920907116104074     4.521077372667832e-19   0.0001648650686890427   6.636626468861713e-06   3.179330642654169e-136  0.0005997100141796689   0.10385298728942871
```
In the above line we see that the causal effect `alpha` is 0.52, with a _P_ value of 4.5x10^-19. The `sigma_y` 
estimate is small (0.00016), but very significant (P: 6.6x10^-136). Indicating a causal effect, as well as a pleiotropic effect.  
```
# non causal effect
region                  var_explained   alpha                   se(alpha)               p(alpha)                sigma_y                 se(sigma_y)             p(sigma_y)              sigma_x                 function_time
2:101515908-103411057   0.99            -0.007902101622960967   0.05244604177702067     0.880235189572735       0.00014733217344441418  5.936191467934854e-06   5.5482348469826166e-136 0.0005383690303972176   0.07675600051879883
```
In the following example, line we see that the causal effect `alpha` is close to zero, with a _P_ value of 0.88. The `sigma_y` 
estimate again is small (0.00016), but very significant (P: 5.5x10^-136). This indicates that the locus is very pleiotropic.

Nb. results may be slightly different in your version, which may be due to the stochastic nature of the methods' inference, and 
or differences in software versions.

## Usage

MR-link 2 accepts full summary statistics files from which it will do the following:
1. Identify all the associated regions from the exposure files
2. For each associated region, load an LD correlation matrix, and make an MR-link 2 estimate

So it is not necessary to pre-select regions. MR-link 2 also performs some rudimentary allele harmonization, but please 
do your own checks beforehand as well.

The `mr_link_2_standalone.py` script uses plink like syntax to as a command. To see all the options, type 
`python3 mr_link_2_standalone.py --help`, which will output the following:
```
usage: mr_link_2_standalone.py [-h] --reference_bed REFERENCE_BED --sumstats_exposure SUMSTATS_EXPOSURE --sumstats_outcome SUMSTATS_OUTCOME --out OUT [--tmp TMP] [--p_threshold P_THRESHOLD] [--region_padding REGION_PADDING] [--maf_threshold MAF_THRESHOLD] [--max_correlation MAX_CORRELATION]
                               [--max_missingness MAX_MISSINGNESS] [--var_explained_grid VAR_EXPLAINED_GRID [VAR_EXPLAINED_GRID ...]] [--continue_analysis] [--no_normalize_sumstats] [--verbose VERBOSE]

MR-link 2: Pleiotropy robust cis Mendelian randomization

options:
  -h, --help            show this help message and exit
  --reference_bed REFERENCE_BED
                        The plink bed file prepend of the genotype file that can be used as an LD reference. Usage is the same as in the plink --bed command
  --sumstats_exposure SUMSTATS_EXPOSURE
                        The summary statistics file of the exposure file. Please see the README file or the example_files folder for examples on how to make these files.
  --sumstats_outcome SUMSTATS_OUTCOME
                        The summary statistics file of the outcome file. Please see the README file or the example_files folder for examples on how to make these files.
  --out OUT             The path where to output results
  --tmp TMP             a prepend on where to save temporary files
  --p_threshold P_THRESHOLD
                        The P value threshold for which select regions. This is the same as the clumping P value threshold
  --region_padding REGION_PADDING
                        The base pair padding (on one side) on each clumped SNP which defines the region in which MR-link 2 will perform its inference
  --maf_threshold MAF_THRESHOLD
                        The minor allele frequency threshold used for clumping, and for calculating the LD matrix.This will not be applied to the summary statistics files
  --max_correlation MAX_CORRELATION
                        The maximum correlation allowed in the LD matrix, if the correlation is higher than this value between a pair of SNPs, will only keep one of them. This value is used to reduce eigenvalue decomposition failures
  --max_missingness MAX_MISSINGNESS
                        This is the maximum amount of individual missingness that is allowed in a summary statistic MR-link 2 can be sensitive to large differences in summary statistic missingness, so by default each SNP association should have at least 0.95 of observations available.
  --var_explained_grid VAR_EXPLAINED_GRID [VAR_EXPLAINED_GRID ...]
                        This field specifies the amount of variance explained of the LD matrix that is used by MR-link 2. You can add onto this field, and all variances explained will be added: --var_explained_grid 0.99 0.999 0.2 0.96 0.1 will perform an MR-link 2 estimate for all these values.
  --continue_analysis   Flag to continue an already started analysis, if specified this will look for a temporary file, and if it present, reuse its results. This can be handy if you have hundreds of associated regions, which can sometimes take a long time to run.
  --no_normalize_sumstats
                        flag to _not_ normalize summary statistics
  --verbose VERBOSE     Set to 1 if you want to read more output, for debugging purposes
```



## The input files for MR-link 2
MR-link 2 requires 2 summary statistics files that are formatted the same and a genotype file in the plink 
.bed/.bim/.fam format

The summary statistics files should at have the following columns, which are mostly self explanatory

```
'pos_name', # SNP names of the variants (must match those of the plink bed file) 
'chromosome', # chromosome names (must match those in the plink bed file) 
'position', # base pair positions (must match those in the plink bed file)
'effect_allele', # The effect allele of the summary statistic
'reference_allele', # The reference allele of the summary statistics
'beta',  # The effect size of the SNP onto the trait
'se', # The standard error of the effect size of the SNP onto the trait
'z', # The Z score of the effect of the SNP on the trait
'pval', # The p values of the SNP 
'n_iids' # the number of individuals available for the SNP effect estimate
```

### Understanding the output of MR-link 2
MR-link 2 is a likelihood function that estimates 3 parameters, and tests 2 of these.
The main parameter of interest is alpha, which represents the causal effect. 
Then, MR-link 2 also simulates a parameter that measures vertical pleiotropy we call sigma_y, this parameter can be 
related to the amount of vertical pleiotropy present in a locus.
Finally, MR-link 2 estimates a parameter sigma_x that represents the amount of heritability there is for the exposure.

An MR-link 2 output file will contain the following columns and an explanation:
##### region
The region that was used for the MR-link 2 inference
##### var_explained
The amount of variance that was kept from the correlation matrix used to correct for LD

#### alpha
The point estimate of alpha, the causal effect
#### se(alpha)       
The standard error of the point estimate of alpha, the standard error of the causal effect
#### p(alpha)        
The p value of the point estimate of alpha
#### sigma_y 
The point estimate of sigma_y, the amount of pleiotropic variance there is in the region (per SNP in the region)
#### se(sigma_y)     
The standard error of the point estimate of sigma_y, the amount of pleiotropic variance there is in the region
#### p(sigma_y)      
The p value of the presence of pleiotropic effect.
#### sigma_x 
The point estimate of the exposure heritability (per SNP in the region)
#### function_time
The time it took to run the MR-link 2 estimate. This does not include preprocessing time, which can be substantial

