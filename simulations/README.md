# Simulations as they were performed in the MR-link-2 paper

This script simulates exposure and outcome summary statistics for Mendelian randomization (MR) analyses,
accounting for SNP linkage disequilibrium (LD), genetic heritability, and causal effects.

## Requirements
Install the required Python libraries:

```bash
pip install numpy scipy pandas argparse
```

### Usage

Run the script with the following example command:

```bash
python script_name.py --n_exp 100000 --n_out 300000 --exp_h2 0.01 --out_h2 0.005 --alpha -0.2 --m_snps 2068 --m_causal_exp 100 --m_causal_out 100 --n_sims 1000 --ref_size 500
```

#### Key Arguments:
    `--n_exp`: Number of samples in exposure GWAS (default: 100000)
    `--n_out`: Number of samples in outcome GWAS (default: 300000)
    `--exp_h2`: Heritability of exposure trait (default: 0.01)
    `--out_h2`: Heritability of outcome trait (default: 0.005)
    `--alpha`: Causal effect of exposure on outcome (default: -0.2)
    `--m_snps`: Number of SNPs (default: 2068, number of SNPs in the LD matrix)
    `--n_sims`: Number of simulations (default: 1000)


 Output:
    Simulated summary statistics in .parq format: `simulated_phenotypes/{parameters}.parq`
    To save space, we only save a limited amount of data, ['beta', 'pvalue', 'N', 'simulation','phenotype_name']
