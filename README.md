# MR-link 2

MR-link 2 is a _cis_ MR method that is pleiotropy robust. 
Inference with MR-link 2 requires **pre-harmonized** summary statistics of an exposure and an outcome, and a 
genotype refererence file. 

Please keep an eye on this space for our paper where we validate MR-link 2 on an extensive set of validation datasets.
The main benefit of MR-link 2 is that it has lower false positive rates than other _cis_ MR methods, while retaining
discriminative ability.

Disclaimer: MR-link 2 is still under development, so it can be subject to change, although we find this unlikely.

### Requirements
MR-link 2 has been tested on macOS X and Linux combined with Python 3.9, 3.10 and 3.11. 
Although not tested, every python3 version from 3.6 onwards should work.    
We require some (standard) packages to be installed, these are: `numpy`, `scipy` and `pandas`.

On top of this we require plink1.9 to be present in your bash path, you can check this by typing `which plink` in your
terminal. If this prints a path, you have plink in your path.


### Example
If you want to test MR-link 2 we have two examples:

This command tests for a causal effect in a region of synthetic data 
```{bash}
python3 mr_link_2_standalone.py \
            --reference_bed genotype_files/reference_cohort \
            --sumstats_exposure testing_files/yes_causal_exposure.txt \
            --sumstats_outcome testing_files/yes_causal_outcome.txt \
            --out example_of_a_causal_effect.txt
```

This command tests for a non-causal effect in a region of synthetic data:
```{bash}
python3 mr_link_2_standalone.py \
            --reference_bed genotype_files/reference_cohort \
            --sumstats_exposure testing_files/non_causal_exposure.txt \
            --sumstats_outcome testing_files/non_causal_outcome.txt \
            --out example_of_a_non_causal_effect.txt
```

After running these two commands (takes about 2 seconds each), they will output two tab separated files with results: 
`example_of_a_causal_effect.txt` and `example_of_a_non_causal_effect.txt`.
```
# causal effect
region  var_explained   alpha   se(alpha)       p(alpha)        sigma_y se(sigma_y)     p(sigma_y)      sigma_x function_time
2:101532661-103480976   0.99    0.5283473895075494      0.05920907116104074     4.521077372667832e-19   0.0001648650686890427   6.636626468861713e-06   3.179330642654169e-136  0.0005997100141796689   0.10385298728942871
```
```
# non causal effect
region  var_explained   alpha   se(alpha)       p(alpha)        sigma_y se(sigma_y)     p(sigma_y)      sigma_x function_time
2:101515908-103411057   0.99    -0.007902101622960967   0.05244604177702067     0.880235189572735       0.00014733217344441418  5.936191467934854e-06   5.5482348469826166e-136 0.0005383690303972176   0.07675600051879883
````
Nb. results may be different in your version, which may be due to the stochastic nature of the methods' inference, and 
or differences in software versions.

### The input files for MR-link 2
MR-link 2 requires 2 summary statistics files that are formatted the same and a genotype file in the plink 
.bed/.bim/.fam format

The summary statistics files should at have the following columns, which are mostly self explanatory
.
```
sumstats_necessary_colnames = {'pos_name',  'chromosome', 'position', 'effect_allele', 'reference_allele',
                                   'beta',  'se', 'z', 'pval',  'n_iids'}
```
Where the `pos_name` identifier should fully match the corresponding variants in the plink bed file.  


### Understanding the output of MR-link 2
MR-link 2 is a likelihood function that estimates 3 parameters, and tests 2 of these.
The main parameter of interest is alpha, which represents the causal effect. 
Then, MR-link 2 also simulates a parameter that measures vertical pleiotropy we call sigma_y, this parameter can be 
related to the amount of vertical pleiotropy present in a locus.
Finally, MR-link 2 estimates a parameter sigma_x that represents the amount of heritability there is for the exposure.

An MR-link 2 output file will contain the following columns and an explanation:
#### region
The region that was used for the MR-link 2 inference
#### var_explained
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

