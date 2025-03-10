[![DOI](https://zenodo.org/badge/688924886.svg)](https://doi.org/10.5281/zenodo.14961110)
# MR-link-2

MR-link-2 is a _cis_ MR method that is pleiotropy robust. 
Inference with MR-link-2 requires **pre-harmonized** summary statistics of an exposure and an outcome, and a 
genotype reference file.  We have validated this method in 3 different real-world datasets of causality. 

Please find details of our validations and more information on the method in our [preprint](https://www.medrxiv.org/content/10.1101/2024.01.22.24301400v1).

If you have any questions or suggestions, feel free to open an issue.
We appreciate everybody trying to use our software, so we try to come back to you as soon as possible!  

### Requirements
MR-link-2 has been tested on MacOS X and Linux combined with Python 3.9, 3.10 and 3.11. 
Although not tested, every Python version from 3.6 onwards should work.    
We require some (standard) python packages to be installed, these are: 
`numpy`, `scipy`, `pandas`, `pyarrow`, `bitarray` and `duckdb`.
If you want to ensure all the tests run, `pytest` is also necessary.
If they haven't been installed, please install these using pip.
In the command line (shell, terminal), type: 

```{bash}
pip3 install numpy scipy pandas bitarray pytest duckdb pyarrow
```
On top of this, we require plink1.9 to be present in your PATH variable. 
Check this by typing `which plink` in your  shell. 
If this prints a path, you have plink installed in your path.

##### Testing if everything works as expected
This repository uses pytest to analyze results
If you want to make sure that everything works as expected, please ensure you have pytest installed and run the 
following command.
```{bash}
pytest tests/*
```
For this you need to have everything installed from the requirements, including `pytest`
If everything passes, you are ready to go! 
If not all the tests work, please open a github issue, and we'll get back to you ASAP.

### Example
If you want to test MR-link-2 we have two examples:

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
region                  var_explained   m_snps_overlap   alpha                   se(alpha)               p(alpha)                sigma_y                 se(sigma_y)             p(sigma_y)              sigma_x                 function_time
2:101532661-103480976	0.99	        1131	         0.4193872544798938	     0.0565878953774011	     1.25110890260147e-13	 0.1549295595449322	     0.00624199019634178	 5.3812527742659674e-136 0.564994403645872	     0.09744000434875488
```
In the above line we see that the causal effect `alpha` is 0.42 , with a _P_ value of 1.3x10^-13. The `sigma_y` 
estimate is large (0.155), and very significant (P: 5.3x10^-136). Indicating a causal effect, as well as a pleiotropic effect.  
```
# non causal effect
region                  var_explained   m_snps_overlap   alpha                   se(alpha)               p(alpha)                sigma_y                 se(sigma_y)             p(sigma_y)              sigma_x                 function_time
2:101515908-103411057	0.99	        1131	         -0.03309665745471561	 0.052906181109434916	 0.5315953131580522	     0.16421519190135686	 0.006618939960069835	 7.011323342257353e-136	 0.5641943699065054	     0.0805368423461914
```
In the following example, line we see that the causal effect `alpha` is close to zero, with a _P_ value of 0.52. The `sigma_y` 
estimate again is large (0.16) and very significant (P: 7.0x10^-136). This indicates that the locus is very pleiotropic.

Alternatively, you can use the `test_commands.sh` script to run a few tests, and have a look at their results.

Nb. results may be slightly different in your version, which may be due to the stochastic nature of the methods' inference, and 
or differences in software versions.

## Usage

MR-link-2 accepts full summary statistics files from which it will do the following:
1. Identify all the associated regions from the exposure files.
2. For each associated region, generate an LD correlation matrix, and make an MR-link-2 estimate.

So it is not necessary to pre-select regions. MR-link-2 also performs some rudimentary allele harmonization, but please 
do your own checks beforehand as well.

The `mr_link_2_standalone.py` script uses plink like syntax to as a command. To see all the options, type 
`python3 mr_link_2_standalone.py --help`, which will output the following:
```
usage: mr_link_2_standalone.py [-h] --reference_bed REFERENCE_BED --sumstats_exposure SUMSTATS_EXPOSURE --sumstats_outcome SUMSTATS_OUTCOME --out OUT [--tmp TMP] [--p_threshold P_THRESHOLD] [--region_padding REGION_PADDING] [--maf_threshold MAF_THRESHOLD] [--max_correlation MAX_CORRELATION]
                               [--max_missingness MAX_MISSINGNESS] [--var_explained_grid VAR_EXPLAINED_GRID [VAR_EXPLAINED_GRID ...]] [--continue_analysis] [--no_normalize_sumstats] [--verbose VERBOSE]

MR-link-2: Pleiotropy robust cis Mendelian randomization

options:
  -h, --help            show this help message and exit
  --reference_bed REFERENCE_BED
                        The plink bed file prepend of the genotype file that can be used as an LD reference. Usage is the same as in the plink --bed command
  --sumstats_exposure SUMSTATS_EXPOSURE
                        The summary statistics file of the exposure file. Please see the README file or the example_files folder for examples on how to make these files.
  
  --sumstats_outcome SUMSTATS_OUTCOME [SUMSTATS_OUTCOME]
                        The summary statistics file of the outcome file. Please see the README file or the example_files folder for examples on how to make these files.
                        We allow multiple outcomes to be analyzed at the same time, include the files separated by spaces
  --out OUT             The path where to output results
  --tmp TMP             Not necessary anymore: a prepend on where to save temporary files DEPRECATED
  --p_threshold P_THRESHOLD
                        The P value threshold for which select regions. This is the same as the clumping P value threshold
  --region_padding REGION_PADDING
                        The base pair padding (on one side) on each clumped SNP which defines the region in which MR-link-2 will perform its inference
  --maf_threshold MAF_THRESHOLD
                        The minor allele frequency threshold used for clumping, and for calculating the LD matrix.This will not be applied to the summary statistics files
  --max_correlation MAX_CORRELATION
                        The maximum correlation allowed in the LD matrix, if the correlation is higher than this value between a pair of SNPs, will only keep one of them. This value is used to reduce eigenvalue decomposition failures
  --max_missingness MAX_MISSINGNESS
                        This is the maximum amount of individual missingness that is allowed in a summary statistic MR-link-2 can be sensitive to large differences in summary statistic missingness, so by default each SNP association should have at least 0.95 of observations available.
  --var_explained_grid VAR_EXPLAINED_GRID [VAR_EXPLAINED_GRID ...]
                        This field specifies the amount of variance explained of the LD matrix that is used by MR-link-2. You can add onto this field, and all variances explained will be added: --var_explained_grid 0.99 0.999 0.2 0.96 0.1 will perform an MR-link-2 estimate for all these values.
  --continue_analysis   Flag to continue an already started analysis, if specified this will look for a temporary file, and if it present, reuse its results. This can be handy if you have hundreds of associated regions, which can sometimes take a long time to run.
  --regions_to_read_at_the_same_time NUMBER OF REGIONS
                        Number of regions to read at the same time from a file. Will reduce I/O times, but can increase memory usage
                        Only change this when you run into memory errors, or if you want to make a lot of comparisons.
  --prespecified_regions PRESPECIFIED_REGIONS
                        Specify which regions to do. format is the following: 
                        `{chr_1}:{start_1}-{end_1},{chr_2}:{start_2}-{end_2}`. 
                        i.e. these are regions that are separated with the comma character.
                        This step will skip the clumping step to identify regions, and will blindly perform MR on 
                        regions that may or may not be associated to the exposure. 
                        This feature has been included to allow investigators to identify associated regions
                        in one discovery cohort, and assess the MR effect in a different cohort. 
                        This is a scenario for instance to remove winners curse.

  --max_snps_in_correlation_matrix
                        How many SNPs are allowed in the correlation matrix, as it can take a lot of time to 
                        do the eigendecompositions and Rscripts otherwise.
                        default: 5250
  
  --run_other_cis_mr_and_coloc
                        This is a flag that you can add to run other MR analyses: MR-IVW, MR-IVW LD and MR-PCA,
                        as well as coloc and SuSIE coloc. Requires R to be installed with the coloc and data.table 
                        packages. Will append columns of these methods' results to the output file.
                        
  --no_normalize_sumstats
                        flag to _not_ normalize summary statistics
  --no_exclude_hla      
                        flag to _not_ exclude the HLA.
                        Please don't use this flag unless you absolutely know what you're doing.
                        The HLA has very long range LD that can interfere with inference of MR-link 2.                  
  --verbose VERBOSE     
                        Set to 1 if you want to read more output, for debugging purposes
```

## The input files for MR-link-2
MR-link-2 requires 2 summary statistics files that are formatted the same and a genotype file in the plink 
.bed/.bim/.fam format. 


**NEW**: MR-link-2 now accepts 2 types of file formats, the 
[GWAS summary statistics file format](https://www.ebi.ac.uk/gwas/docs/summary-statistics-format) and for legacy purposes
a file format that was specifically generated for this program.

If you don't have your summary statistics in the GWAS summary statistics format, you can format it still in the following way:
The summary statistics files should be a tab ('\t') separated file at have the following columns, which are mostly self 
explanatory.
```
'pos_name', # SNP names of the variants (must match those of the plink bed file) 
'chromosome', # chromosome names (must match those in the plink bed file) 
'position', # base pair positions (must match those in the plink bed file)
'effect_allele', # The effect allele of the summary statistic
'reference_allele', # The reference allele of the summary statistics
'beta',  # The effect size of the SNP onto the trait
'se', # The standard error of the effect size of the SNP onto the trait
'z', # The Z score of the effect of the SNP on the trait, if you don't have this, you can calculate it using the formula "z = beta / se" 
'pval', # The p values of the SNP 
'n_iids' # the number of individuals available for the SNP effect estimate
```

## Understanding the output of MR-link-2
MR-link-2 is a likelihood function that estimates 3 parameters, and tests 2 of these.
The main parameter of interest is alpha, which represents the causal effect. 
Then, MR-link-2 also simulates a parameter that measures horizontal pleiotropy we call sigma_y, this parameter can be 
related to the amount of vertical pleiotropy present in a locus.
Finally, MR-link-2 estimates a parameter sigma_x that represents the amount of heritability there is for the exposure.

### results files
After a succesful estimate, MR-link-2 output files contain the following columns and an explanation:
##### region
The region that was used for the MR-link-2 inference
##### var_explained
The amount of variance that was kept from the correlation matrix used to correct for LD
##### m_snps_overlap
The number of SNPs in the region on which the estimate was based.
##### alpha
The point estimate of alpha, the causal effect
##### se(alpha)       
The standard error of the point estimate of alpha, the standard error of the causal effect
##### p(alpha)        
The p value of the point estimate of alpha
##### sigma_y 
The point estimate of sigma_y, the amount of pleiotropic variance there is in the region 
##### se(sigma_y)     
The standard error of the point estimate of sigma_y, the amount of pleiotropic variance there is in the region
##### p(sigma_y)      
The p value of the presence of pleiotropic effect.
##### sigma_x 
The point estimate of the exposure heritability in the region
##### function_time
The time it took to run the MR-link-2 estimate. This does not include preprocessing time, which can be substantial


### _no_estimate files
The no_estimate files will be written when something went wrong in the MR-link 2 function. This is usually caused by incorrectly harmonized summary statistics.
Please don't hesitate to contact us if you have questions about a specific analysis. 

### _no_region files
The `_no_region` files will be written when, after matching the SNPs within the two summary statistics files, there is no exposure genetic association available at the  `--p_threshold` instrument selection threshold.

## Frequently asked questions 

### How do I improve performance of my analysis
Two tips on how to improve the performance of a large analysis run:
1. Use parquet files as input for MR-link-2
2. Use multiple phenotypes as outcomes for the same exposure.
This will make sure that there is the least amount of double work to do that is possible.

### I need to input harmonized summary statistics, what does this mean?
How to correctly harmonize a study depends on what you're studying. 
For the analysis in our original publication, we performed the following steps:

```
We processed summary statistics of all studies in the same way: 
First, if necessary, we lifted over summary statistics into human chromosome build 37 using UCSCs liftover tool (https://genome.ucsc.edu/cgi-bin/hgLiftOver) combined with their chain files (https://hgdownload.soe.ucsc.edu/downloads.html). 
Then, we include SNP variants that have LD information available, by overlapping the variants (based on chromosome, position and alleles) present in the summary statistics file with the variants in our LD reference (UK10k). Due to potential strand inconsistencies, palindromic SNPs were removed. 
Genetic associations are retained if they have at least a minor allele frequency of 0.5% in the UK10K LD reference and if the variant has been measured in at least 95% of the maximum number of measured individuals (if the information was available).
```

It is not necessary to remove palindromic SNPs from your cohort if you're confident that they are measured on the correct strand. 
As we were using publicly available data, we didn't want to make that assumption, as it can be detrimental to the MR-link 2 analysis.  

### Why are the effect sizes different from my other MR analysis?
To retain statistical stability, we perform a normalization step of the effect sizes. This normalization step changes the units of the genetic associations into 'variance explained' units. 
If summary statistics information is in another unit, this will result in changes of the eventual causal effect estimate compared to more classical MR analysis. 

If you _really_ do not want this normalization to happen, you can use the `--no_normalize_sumstats` option. Please be careful when doing this though, as the MR-link 2 results here have not been tested, and we provide no guarantees on the accuracy of the estimates.
