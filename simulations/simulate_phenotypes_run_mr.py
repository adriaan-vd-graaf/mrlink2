import numpy as np
import scipy.stats
import argparse
import time
import scipy.io
import pandas as pd
import os

def select_causal_snps(exposure_1_n_causal, exposure_2_n_causal,
                       overlapping_causal_snps, ld_mat,
                       n_reps,
                       upper_ld_bound,
                       lower_ld_bound,
                       max_failures=5):
    """
    Selects causal SNPs for two exposures over multiple repetitions, ensuring linkage disequilibrium (LD)
    bounds are respected.

    This function performs SNP selection by calling a base SNP selection function (`select_causal_snps_base`)
    for a specified number of repetitions (`n_reps`), while enforcing upper and lower LD bounds.
    If the SNP selection process fails (e.g., due to LD bounds not being met), it retries up to
    `max_failures` before raising an error.

    :param exposure_1_n_causal:
        The number of causal SNPs to be selected for the first exposure.
    :param exposure_2_n_causal:
        The number of causal SNPs to be selected for the second exposure (usually the pleiotropic exposure).
    :param overlapping_causal_snps:
        The number of overlapping causal SNPs shared between the two exposures.
    :param ld_mat:
        A NumPy array representing the linkage disequilibrium (LD) matrix, which governs the correlation between SNPs.
    :param n_reps:
        The number of repetitions to run for selecting SNPs.
    :param upper_ld_bound:
        The upper bound on LD that is allowed when selecting SNPs.
    :param lower_ld_bound:
        The lower bound on LD that is allowed when selecting SNPs.
    :param max_failures:
        The maximum number of failures allowed before raising an exception (default is 5).

    :return:
        A tuple containing two NumPy arrays:
            - `exposure_snps`: A boolean array indicating the selected SNPs for exposure across all repetitions.
            - `outcome_snps`: A boolean array indicating the selected SNPs for outcome across all repetitions.

    :raises ValueError:
        If the SNP selection process fails more than `max_failures` times.
    """
    iteration = 0
    failures = 0

    exposure_snps = np.zeros((ld_mat[0].shape[0], n_reps), dtype=bool)
    outcome_snps = np.zeros((ld_mat[0].shape[0], n_reps), dtype=bool)

    while iteration < n_reps:

        if failures >= max_failures:
            raise ValueError('Unable to perform correct instrument selection')

        try:
            exposure_index, outcome_index = select_causal_snps_base(exposure_1_n_causal,
                                                                    exposure_2_n_causal,
                                                                    overlapping_causal_snps,
                                                                    ld_mat,
                                                                    upper_ld_bound,
                                                                    lower_ld_bound)

            exposure_snps[exposure_index, iteration] = True
            outcome_snps[outcome_index, iteration] = True

            iteration += 1
            failures = 0

        except:
            failures += 1

    return exposure_snps, outcome_snps


def select_causal_snps_base(exposure_1_n_causal, exposure_2_n_causal,
                            overlapping_causal_snps,
                            ld_mat,
                            upper_ld_bound,
                            lower_ld_bound):
    """
    This is a helper function for simulate_phenotypes to choose which indices will be considered for causal SNPs.
    It's in it's own function because it can fail to find SNPs from a selection, so we make the parent function more
    robust to failures.

    It implements steps 1 and 2 in the simulation function:

    1.
        exposure 1 SNPs are chosen, and they are chosen not to be in higher ld than 0.95 R ^ 2

    2.
        exposure 2 SNPs are chosen in two ways. `overlapping_causal_snps` are chosen without replacement from the
        exposure 1 SNPs, the rest are chosen from SNPs that are in LD with exposure 1.

    :param exposure_1_n_causal:
    :param exposure_2_n_causal:
    :param overlapping_causal_snps:
    :param exposure_geno:
    :param exposure_ld:
    :param upper_ld_bound:
    :param known_exposure_lower_ld_bound:
    :param lower_ld_bound:
    :return:
    """
    # make sure the exposure Causal SNPs are not in > upper_ld_bound R^2 with another.

    r_sq_mat = ld_mat ** 2

    exposure_1_causal_snps = np.zeros([exposure_1_n_causal], dtype=int)
    i = 0
    # choose the first causal SNP
    exposure_1_causal_snps[i] = np.random.choice(ld_mat.shape[0], 1, replace=False)
    i += 1

    # i indicates how many snps have been chosen
    while i < exposure_1_n_causal:
        # choosing SNPS that are not in full linkage (R**2 > 0.95) with one another, but do share some linkage.
        lower_than_upper = np.sum(r_sq_mat[exposure_1_causal_snps[:i], :] < upper_ld_bound, axis=0) == i

        if lower_ld_bound >= 0.0:  # this part tries to find causal exposure variants in high LD.
            higher_than_lower = np.sum(r_sq_mat[exposure_1_causal_snps[:i], :] >= lower_ld_bound,
                                       axis=0) >= 1
            if sum(higher_than_lower) == 0:
                raise ValueError("Could not find known exposure variants within the lower LD bound")
            if sum(higher_than_lower & lower_than_upper) == 0:
                raise ValueError("Could not find known exposure variants within the lower AND higher LD bound")
            choice = np.random.choice(np.where(higher_than_lower & lower_than_upper)[0], 1, replace=False)
        else:
            choice = np.random.choice(np.where(lower_than_upper)[0], 1, replace=False)

        exposure_1_causal_snps[i] = choice
        i += 1

    """
    Step 2.
        Select overlapping.
        If there are non overlapping SNPs left, select SNPs in LD with exposure 1.
    """
    # Here, the snps of exposure are overlapped, depending on the `overlapping_causal_snps` number, could be zero.
    exposure_2_overlapping_snps = np.random.choice(exposure_1_causal_snps, overlapping_causal_snps, replace=False)

    """
        Here we select SNPs for exposure 2, that are in LD with exposure 1.
        But we do not look for LD with already overlapping causal snps for exposure 1.

        Only if there are no overlapping snps left for exposure 2 will we just sat our exposure 2 causal snps are done.
    """
    if exposure_2_n_causal - overlapping_causal_snps > 0:

        exposure_1_causal_snps_non_overlapping = np.asarray(
            [x for x in exposure_1_causal_snps if x not in exposure_2_overlapping_snps],
            dtype=int)

        # permute this vector, so ordering is random.
        permuted_exposure_1_causal_non_overlapping = np.random.permutation(exposure_1_causal_snps_non_overlapping)
        exposure_2_ld_snps = np.zeros((exposure_2_n_causal - overlapping_causal_snps), dtype=int)

        # iterate over all causal snps, and choose one that is within the ld window
        i = 0
        while i < (exposure_2_n_causal - overlapping_causal_snps):

            ld_with_causal = r_sq_mat[permuted_exposure_1_causal_non_overlapping[i], :]

            try:
                chosen_indice = np.random.choice(
                    np.where((ld_with_causal > lower_ld_bound) & (ld_with_causal < upper_ld_bound))[0], 1)
                exposure_2_ld_snps[i] = chosen_indice
                i += 1
            except:
                raise ValueError("Could not find SNPs in an ld window around the exposure snp.")

        exposure_2_causal_snps = np.concatenate((exposure_2_overlapping_snps, exposure_2_ld_snps))

    else:
        exposure_2_causal_snps = exposure_2_overlapping_snps

    return exposure_1_causal_snps, exposure_2_causal_snps


def simulate_exp_outcome_sumstats(r_mat,
                                  causal_effect=-0.2,
                                  n_exp=10_000, n_out=300_000,
                                  exp_h2=0.01, out_h2=5e-3,
                                  m_causal_exp=5, m_causal_out=5,
                                  n_reps=1,
                                  min_ld_between_causal=0.0,
                                  max_ld_between_causal=1.0
                                  ):
    """
        Simulates exposure and outcome summary statistics for Mendelian randomization analyses.

        This function simulates causal SNPs for both exposure and outcome traits, ensuring LD constraints
        (if any) are met, and generates the corresponding summary statistics.
        The function accounts for genetic heritability (h2) and the causal effect of the exposure on the outcome.

        :param r_mat:
            A NumPy array representing the linkage disequilibrium (LD) matrix, which defines the correlation between SNPs.
        :param causal_effect:
            The true causal effect of the exposure on the outcome. Default is -0.2.
        :param n_exp:
            The number of samples in the exposure GWAS. Default is 10,000.
        :param n_out:
            The number of samples in the outcome GWAS. Default is 300,000.
        :param exp_h2:
            The heritability of the exposure trait. Default is 0.01.
        :param out_h2:
            The heritability of the outcome trait. Default is 5e-3.
        :param m_causal_exp:
            The number of causal SNPs for the exposure trait. Default is 5.
        :param m_causal_out:
            The number of causal SNPs for the outcome trait. Default is 5.
        :param n_reps:
            The number of repetitions for the simulation. Default is 1.
        :param min_ld_between_causal:
            The minimum linkage disequilibrium (LD) allowed between causal SNPs. Default is 0.0.
        :param max_ld_between_causal:
            The maximum linkage disequilibrium (LD) allowed between causal SNPs. Default is 1.0.

        :return:
            A tuple containing the following elements:
            - `bX`: Simulated exposure summary statistics (NumPy array of size m_snps x n_reps).
            - `bY`: Simulated outcome summary statistics (NumPy array of size m_snps x n_reps).
            - `exposure_causal_snps`: A boolean NumPy array indicating the selected causal SNPs for exposure.
            - `outcome_causal_snps`: A boolean NumPy array indicating the selected causal SNPs for outcome.

        :raises ValueError:
            If the number of causal SNPs is larger than the available number of SNPs in the LD matrix.
        """

    m = r_mat.shape[0]
    if m < m_causal_exp or m < m_causal_exp:
        raise ValueError('number of causal SNPs needs to be smaller than the amount of causal SNPs')

    sigma_exp = np.sqrt(exp_h2 / m_causal_exp)
    sigma_out = np.sqrt(out_h2 / m_causal_out)

    if max_ld_between_causal != 1.0 or min_ld_between_causal != 0.0:
        exposure_causal_snps, outcome_causal_snps = select_causal_snps(m_causal_exp,
                                                                       m_causal_out,
                                                                       overlapping_causal_snps=0,
                                                                       ld_mat=r_mat,
                                                                       n_reps=n_reps,
                                                                       upper_ld_bound=max_ld_between_causal,
                                                                       lower_ld_bound=min_ld_between_causal)
    else:
        exposure_causal_snps = np.zeros((m, n_reps), dtype=bool)
        outcome_causal_snps = np.zeros((m, n_reps), dtype=bool)

        for i_rep in range(n_reps):
            exposure_causal_snps[np.random.choice(m_snps, m_causal_exp, replace=False), i_rep] = True
            outcome_causal_snps[np.random.choice(m_snps, m_causal_out, replace=False), i_rep] = True

    gamX = exposure_causal_snps * np.random.normal(0, sigma_exp, (m, n_reps))
    gamY = outcome_causal_snps * np.random.normal(0, sigma_out, (m, n_reps))

    betX = c_r_mat @ gamX + np.random.multivariate_normal(np.zeros(m), c_r_mat * (1 - exp_h2) / n_exp, size=n_reps).T

    betY = c_r_mat @ (causal_effect * gamX + gamY) + \
           np.random.multivariate_normal(np.zeros(m_snps),
                                         c_r_mat * (1 - causal_effect ** 2 * exp_h2 - out_h2) / n_out, size=n_reps).T

    seX = np.ones((m_snps, n_reps)) * np.sqrt((1 - exp_h2) / n_exp)
    seY = np.ones((m_snps, n_reps)) * np.sqrt((1 - causal_effect ** 2 * exp_h2 - out_h2) / n_out)

    bX = (betX / seX) / np.sqrt(n_exp)
    bY = (betY / seY) / np.sqrt(n_out)

    return bX, bY, exposure_causal_snps, outcome_causal_snps


def wishart_rvs(df, sigma, n_reps=1):
    d = np.linalg.cholesky(sigma).T
    n = int(sigma.shape[0])

    if n_reps == 1:
        to_return = np.zeros(sigma.shape, dtype=float)
    else:
        to_return = np.zeros((sigma.shape[0], sigma.shape[1], n_reps), dtype=float)

    for j in range(n_reps):

        if df < 81 + n:
            x = np.random.normal(size=(df, n)) @ d
        else:
            a = np.diag(np.sqrt(np.random.chisquare(df - np.arange(0, n))))
            a[np.triu_indices(n, 1)] = np.random.normal(size=(n * (n - 1)) // 2)
            x = a[:, 0:n] @ d

        if n_reps == 1:
            to_return = x.T @ x
        else:
            to_return[:, :, j] = x.T @ x

    return to_return


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--n_exp', type=int, default=100_000)
    parser.add_argument('--n_out', type=int, default=300_000)

    parser.add_argument('--exp_h2', type=float, default=0.1)
    parser.add_argument('--out_h2', type=float, default=1e-5)

    parser.add_argument('--alpha', type=float, default=0.00)  # causal effect.
    parser.add_argument('--m_snps', type=int, default=2068)

    parser.add_argument('--m_causal_exp', type=int, default=100)
    parser.add_argument('--m_causal_out', type=int, default=100)

    parser.add_argument('--min_ld_between_causal', type=float, default=0.0)
    parser.add_argument('--max_ld_between_causal', type=float, default=1.0)

    parser.add_argument('--n_sims', type=int, default=1000)
    parser.add_argument('--ref_size', type=int, default=500)

    args = parser.parse_args()

    """
    Parameters of this simulation run
    """

    n_exp = args.n_exp
    n_out = args.n_out
    exp_h2 = args.exp_h2
    out_h2 = args.out_h2
    alpha = args.alpha
    m_snps = args.m_snps
    m_causal_exp = args.m_causal_exp
    m_causal_out = args.m_causal_out
    min_ld_between_causal = args.min_ld_between_causal
    max_ld_between_causal = args.max_ld_between_causal

    n_sims = args.n_sims
    ref_size = args.ref_size

    folder_path = 'simulated_phenotypes'
    parameter_name = f'{alpha=}_{exp_h2=}_{out_h2=}_{n_exp=}_{n_out=}_{m_snps=}' \
                     f'_{n_sims=}_{ref_size=}_{min_ld_between_causal=}' \
                     f'_{max_ld_between_causal}_{m_causal_exp=}_{m_causal_out=}'
    out_file = f'{folder_path}/{parameter_name}.parq'

    print(f'Starting: {parameter_name}')

    sigma_exp = np.sqrt(exp_h2 / (m_causal_exp))
    sigma_out = np.sqrt(out_h2 / (m_causal_out))

    ##process the correlation matrix to have more or less snps

    full_r_mat = np.asarray(scipy.io.loadmat('Cmatrix.mat')['Cr'])
    snps_in_full_r_mat = full_r_mat.shape[0]

    r_mat = full_r_mat[:m_snps, :][:, :m_snps]

    w_grid = [0.0, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 5e-1]

    # Small regularization parameter to make this matrix invertible
    for w in w_grid:
        c_r_mat = (1 - w) * r_mat + w * np.eye(m_snps)
        if all(np.linalg.eigvals(c_r_mat) > 0):
            print(f'chosen {w} regularization constant')
            w = w
            break

    all_optim_results = [0.0] * n_sims
    results_mat = np.zeros((n_sims, 7), float)
    chisq_stats = np.zeros(n_sims, float)

    phenotype_simulation_time = time.time()
    bXs, bYs, exposure_causal, outcome_causal = simulate_exp_outcome_sumstats(c_r_mat,
                                                                              causal_effect=alpha,
                                                                              n_exp=n_exp, n_out=n_out,
                                                                              exp_h2=exp_h2, out_h2=out_h2,
                                                                              m_causal_exp=m_causal_exp,
                                                                              m_causal_out=m_causal_out,
                                                                              n_reps=n_sims,
                                                                              min_ld_between_causal=min_ld_between_causal,
                                                                              max_ld_between_causal=max_ld_between_causal
                                                                              )
    phenotype_simulation_time = time.time() - phenotype_simulation_time

    pXs = 2 * scipy.stats.norm.sf(abs(bXs * np.sqrt(n_exp)))
    pYs = 2 * scipy.stats.norm.sf(abs(bYs * np.sqrt(n_out)))
    mafs = np.asarray([0.5] * m_snps)

    all_dfs = []


    for i_pheno in range(bXs.shape[1]):
        exp_tmp_df = pd.DataFrame()
        exp_tmp_df['beta'] = bXs[:, i_pheno] * mafs * (1 - mafs)
        exp_tmp_df['p_values'] = pXs[:, i_pheno]
        exp_tmp_df['N'] = n_exp
        exp_tmp_df['simulation'] = i_pheno
        exp_tmp_df['phenotype_name'] = f'exposure_{i_pheno}'

        out_tmp_df = pd.DataFrame()
        out_tmp_df['beta'] = bYs[:, i_pheno] * mafs * (1 - mafs)
        out_tmp_df['p_values'] = pYs[:, i_pheno]
        out_tmp_df['N'] = n_out
        out_tmp_df['simulation'] = i_pheno
        out_tmp_df['exp_or_outcome'] = 'outcome'

        all_dfs.append(exp_tmp_df)
        all_dfs.append(out_tmp_df)

    full_df = pd.concat(all_dfs)


    # Check if the folder exists, if not, create it
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

    full_df.to_parquet(out_file, index=False, compression='zstd')







