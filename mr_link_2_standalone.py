import numpy as np
import scipy.optimize
import scipy.stats
import re
import pandas as pd
import argparse
import os
import time
import subprocess
import copy
import scipy
import scipy.io
from typing import Any, Tuple, Dict
import tempfile
import bitarray

import pyarrow.dataset as ds
import duckdb

from collections import OrderedDict

class PlinkGenoReader:
    """
    TODO documentation
    """

    def __init__(self, plink_bed_prepend):

        self.prepend = plink_bed_prepend
        self.bed_loc = f'{self.prepend}.bed'
        self.bim_loc = f'{self.prepend}.bim'
        self.fam_loc = f'{self.prepend}.fam'

        for filename in [self.bed_loc, self.bim_loc, self.fam_loc]:
            if not os.path.exists(filename):
                raise ValueError(f'Error Could not find {filename=} when reading the genotype reference')

        self.n_individuals = 0
        with open(self.fam_loc, 'r') as f:
            for line in  f:
                if not line in ['\n', '\r\n']:
                    self.n_individuals+=1

        last_chromosome_position = ('0', 0)

        self._decoder = {0: bitarray.bitarray('00'),  #homoz A1 (usually minor)
                         1: bitarray.bitarray('01'),  #heteroz
                         2: bitarray.bitarray('11'),  #homoz A2 (usually major)
                         3: bitarray.bitarray('10'),  # missing
                         }

        ## Ensured that this is ordered, so we can find a range.
        self.bim_data = OrderedDict()
        self.snp_name_to_variant_order = {}
        with open(self.bim_loc, 'r') as f:
            for i, line in enumerate(f):
                split = line.split()
                self.bim_data[split[1]] = (split[0], split[3], split[4], split[5])
                self.snp_name_to_variant_order[split[1]] = i

                if split[0] != last_chromosome_position[0]:
                    last_chromosome_position = (split[0], int(split[3]))
                elif last_chromosome_position[1] > int(split[3]):
                    raise ValueError(f'Error: bim file is not ordered correctly at {split}')
                else:
                    last_chromosome_position = (split[0], int(split[3]))

        self.n_variants = len(self.bim_data)

    def read_list_of_snps_into_geno_array(self, snps_to_load: set, missing_encoding=3, dtype=float):
        """
        # TODO fix documentation.
        Reads a bed file into a numpy array.
        """
        self._missing_encoding = missing_encoding

        snps_to_load = list(snps_to_load)
        variant_idxes_to_load = [self.snp_name_to_variant_order[x] for x in snps_to_load]
        n_variants = len(snps_to_load)

        bytes_per_variant = int(np.ceil(self.n_individuals / 4))
        byte_array = bytearray(b"\0" * bytes_per_variant * len(snps_to_load))

        with open(self.bed_loc, "rb") as bed_file:
            magick = bed_file.read(3)
            if magick != b'l\x1b\x01':
                raise ValueError("Plink file magic string is not correct.")

            byte_array_position = 0
            for variant_idx in variant_idxes_to_load:
                bed_file.seek((variant_idx * bytes_per_variant) + 3, 0)
                byte_array[byte_array_position: byte_array_position+bytes_per_variant] = bed_file.read(bytes_per_variant)
                byte_array_position = byte_array_position + bytes_per_variant

        genotypes = np.zeros((self.n_individuals, n_variants), dtype=np.uint8)

        offset = 0
        for i in range(n_variants):
            bits = bitarray.bitarray(endian="little")
            bits.frombytes(byte_array[offset:(offset+bytes_per_variant)])
            array = list(bits.decode(self._decoder))[:self.n_individuals]
            genotypes[:,i] = array
            offset+=bytes_per_variant

        genotypes = np.array(genotypes, dtype=dtype)
        genotypes[genotypes == 3] = missing_encoding

        #convert genotypes to minor allele is coded as one.
        tmp_geno = copy.deepcopy(genotypes)
        tmp_geno[genotypes == 2] = 0
        tmp_geno[genotypes == 0] = 2
        genotypes = tmp_geno


        self.genotypes = genotypes
        return genotypes, snps_to_load


class StartEndRegion:
    """
    Class that implements a genomic region, which is really a stretch of base pairs.

    Attributes
    ----------

    chromosome: str
        string representing the chromosome

    start: int
        start base pair position of the genomic region.

    end: int
        end position of the genomic region.

    Methods
    -------

    position_in_region(self, chr, position):
        returns a bool on if the single position is within the region.

    snp_in_region(chr, position)
        returns a bool on if the single position is within the region.
        synonym method of position in region.

    snp_object_in_region(snp_object)
        returns a bool on if the object of class SNP is within the region.

    region_overlaps(other_region)
        returns a bool if another of region of the same class overlaps.

    """

    def __init__(self, *args, **kwargs):

        if type(args[0]) == list and len(args) == 1:
            self.chromosome = str(args[0][0])
            self.start = int(args[0][1])
            self.end = int(args[0][2])
        elif len(args) == 3 and type(args[0]) != list:
            self.chromosome = str(args[0])
            self.start = int(args[1])
            self.end = int(args[2])
        elif type(args[0]) == str and len(args) == 1:
            self.chromosome, _, rest = args[0].partition(":")
            self.start, self.end = [int(x) for x in rest.split('-')]
        elif isinstance(args[0], StartEndRegion) and len(args) == 1:
            self.chromosome = args[0].chromosome
            self.start = args[0].start
            self.end = args[0].end
        else:
            raise ValueError(
                "Constructor only accepts a list [<chr>, <start>, <end>], three arguments (<chr>, <start>, <end>) "
                "or a string formatted as 'chr:start-end' ")

        # Runtime checks.
        if self.start > self.end:
            raise RuntimeError("Region cannot have a start position smaller than an end position")

        if self.start < 0 or self.end < 0:
            raise RuntimeError("Region cannot have negative positions.")

    def position_in_region(self, chr, position):
        """
        return if the position is in the region.

        :param chr: chromosome str or castable to str
        :param position: position int or castable to int
        :return: bool
        """
        return (self.chromosome == str(chr)) and (self.start <= int(position) <= self.end)

    def snp_in_region(self, chr, position):
        """
        return if the position is in the region.
        synonym method of position_in_region method.

        :param chr: chromosome str or castable to str
        :param position: position int or castable to int
        :return: bool
        """
        return self.position_in_region(chr, position)

    def snp_object_in_region(self, snp_object):
        """

        :param snp_object: object of type SNP (this package)
        :return: True or false if snp in region
        """
        return self.snp_in_region(snp_object.chromosome, snp_object.position)

    def region_overlaps(self, other_region):
        if self.chromosome == other_region.chromosome:
            # this may contain an error, and could be done more efficiently.
            if self.start <= other_region.start <= self.end \
                    or self.start <= other_region.end <= self.end \
                    or other_region.start <= self.start <= other_region.end \
                    or other_region.start <= self.end <= other_region.end:
                return True

        return False

    def __str__(self):
        return '{}:{}-{}'.format(self.chromosome, self.start, self.end)

    def __lt__(self, other):

        if not other.__class__ is self.__class__:
            return NotImplemented

        if not self.chromosome == other.chromosome:
            try:
                return int(self.chromosome) < int(other.chromosome)
            except:
                return self.chromosome < other.chromosome
        return self.start < other.start

    def __gt__(self, other):

        if not other.__class__ is self.__class__:
            return NotImplemented

        if not self.chromosome == other.chromosome:
            try:
                return int(self.chromosome) > int(other.chromosome)
            except:
                return self.chromosome > other.chromosome

        return self.start > other.end

    def __contains__(self, item):
        if isinstance(item, StartEndRegion):
            return self.region_overlaps(item)
        else:
            raise ValueError("Only classes (or inheritance allowed:) SNP.variant or gene_regions.StartEndRegion")

    def __repr__(self):
        return f'{self.chromosome}:{self.start}-{self.end}'


class StartEndRegions:
    """
    This class contains multiple start end regions

    Attributes
    ----------

    gene_regions: list of StartEndRegion objects

    Methods
    -------

    in_gene_regions(self, chr, position)
        identifies if a position is in any of the regions.

    make_non_overlapping_regions(self)
        combines the regions into contiguous non-overlapping regions.

    """

    def __init__(self, list_of_regions):
        self.gene_regions = list([StartEndRegion] * len(list_of_regions))
        i = 0
        for region in list_of_regions:
            tmp = StartEndRegion(region)
            self.gene_regions[i] = tmp
            i += 1

    def in_gene_regions(self, chr, position):
        """
        Identify if a snp is in any of the gene regions.

        :param chr: chromosome castable to str
        :param position: position castable to int
        :return: boolean
        """
        for i in self.gene_regions:
            if i.snp_in_region(chr, position):
                return True
        return False

    def make_non_overlapping_regions(self):
        """
        Combines all overlapping regions, and turns them into one big region.

        :return: new instance of StartEndRegions containing contiguous, non overlapping regions
        """
        sorted_regions = sorted(self.gene_regions)

        combined = False
        non_overlapping = []
        index = 0
        tmp_region = sorted_regions[index]

        while index < len(sorted_regions) - 1:

            if tmp_region.region_overlaps(sorted_regions[index + 1]):
                tmp_region.end = sorted_regions[
                    index + 1].end  # Assumes this is sorted, so we don't look at the start region
            else:
                non_overlapping.append(tmp_region)
                tmp_region = sorted_regions[index + 1]
            index += 1

        # finally, add the last tmp region to the file.
        non_overlapping.append(tmp_region)

        return StartEndRegions([str(x) for x in non_overlapping])

    def __next__(self):
        self.i += 1
        if self.i > len(self.gene_regions):
            raise StopIteration
        else:
            return self.gene_regions[self.i - 1]

    def __iter__(self):
        self.i = 0
        return self

    def __repr__(self):
        return (f'StartEndRegions with {len(self.gene_regions)} regions across chromosome(s) '
                f'{sorted(set([x.chromosome for x in self.gene_regions]))}')

    def __contains__(self, item):
        if isinstance(item, StartEndRegion):
            for region in sorted(self.gene_regions):
                if region > item:
                    return False
                elif region in item:
                    return True
                else:
                    continue
        else:
            raise ValueError("Only classes (or inheritance allowed:) SNP.variant or gene_regions.StartEndRegion")


def identify_regions(sumstats_exposure: str,
                     bed_prepend: str,
                     plink_p_threshold: float,
                     plink_maf_threshold: float,
                     padding: int,
                     plink_tmp_prepend: str,
                     r_sq_threshold=0.01,
                     verbosity_level=0) -> StartEndRegions:
    """
    This is a function that takes a summary_statistics file, and the location of a bed file
    and clumps the most associated variants using plink.

    :param sumstats_exposure: pandas.DataFrame containing at least the columns rsid, pval, where .rsid should
                              match the variants in the bed file.
    :param bed_prepend: str containing the path prepend of plink genotype files in the .bed, .bim, .fam file
    :param plink_p_threshold: float the p value threshold used for clumping
    :param plink_maf_threshold: the minor allele frequency threshold used for further filtering
    :param padding: int, the number of base poirs of padding of the plink region, and the combination of clumps/
    :param plink_tmp_prepend: str a prepend of a file location where files are stored
    :param r_sq_threshold: float, the rˆ2 threshold for clumping, default = 0.01
    :param verbosity_level: int the verbosity level. Will return the output of plink as stdout if it is not Null.
    :return: StartEndRegions that were found using all the clumps that were present.
    """

    stderr = subprocess.DEVNULL if not verbosity_level else None
    stdout = subprocess.DEVNULL if not verbosity_level else None

    files_to_remove = []
    subprocess.run(['plink',
                    '--bfile', bed_prepend,
                    '--maf', f'{plink_maf_threshold}',
                    '--clump', sumstats_exposure,
                    '--clump-p1', f'{plink_p_threshold:.2e}',
                    '--clump-r2', f'{r_sq_threshold}',
                    '--clump-kb', f'{padding / 1000}',
                    '--clump-snp-field', 'rsid',
                    '--clump-field', 'p_value',
                    '--out', f'{plink_tmp_prepend}_clumping_results'
                    ], check=True, stderr=stderr, stdout=stdout)
    files_to_remove += [f'{plink_tmp_prepend}_clumping_results.{x}' for x in ['clumped', 'log', 'nosex']]

    if not os.path.exists(f'{plink_tmp_prepend}_clumping_results.clumped'):
        return None

    all_regions = []
    try:
        with open(f'{plink_tmp_prepend}_clumping_results.clumped') as f:
            f.readline()  # header
            for line in f:
                if line in {'', '\n'}:
                    continue
                split = line.split()
                chromosome = split[0]
                position = split[3]
                data = [str(chromosome), int(int(position) - padding if int(position) - padding > 0 else 0),
                        int(position) + padding]
                all_regions.append(
                    StartEndRegion(data)
                )
        all_regions = sorted(all_regions)

    except Exception as x:
        raise ValueError(f'Couldn\'t parse the clumped file with error {x}')

    combined_regions = StartEndRegions(all_regions).make_non_overlapping_regions()

    for filename in files_to_remove:
        if os.path.exists(filename):
            os.remove(filename)

    return combined_regions.gene_regions


def mr_link2_loglik_reference_v0(th: np.ndarray, lam: np.ndarray,
                                 c_x: np.ndarray, c_y: np.ndarray,
                                 n_x: float, n_y: float) -> float:
    """
    The MR-link2 log likelihood function. This function calculates -1 * likelihood of three parameters:
    alpha, sigma_x and sigma_y.
    Designed to be used in optimization algorithms like those in scipy.minimize

    :param th:
        List or numpy array of floats with the parameters to optimize first is alpha, second the
        1 /  exposure heritability and third the 1/ outcome heritability.
    :param lam:
        np.ndarray of selected eigenvalues of the cX and cY parameters.
    :param c_x:
        The dot product of the selected eigenvectors and summary statistics vector of the exposure
    :param c_y:
        The dot product of the selected eigenvectors and summary statistics vector of the outcome
    :param n_x:
        The number of individuals in the exposure dataset
    :param n_y:
        The number of individuals in the outcome dataset
    :return:
        a single float that contains the likelihood of the parameters theta.
    """

    n_x = np.longdouble(n_x)
    n_y = np.longdouble(n_y)
    a = np.longdouble(th[0])

    tX = abs(np.longdouble(th[1]))
    tY = abs(np.longdouble(th[2]))
    lam = np.asarray(lam, dtype=np.longdouble)
    c_x = np.asarray(c_x, dtype=np.longdouble)
    c_y = np.asarray(c_y, dtype=np.longdouble)

    Dxx = 1. / (((a ** 2 * n_y + n_x) * lam + tX) - a ** 2 * n_y ** 2 * (lam ** 2) / (n_y * lam + tY))
    Dxy = -Dxx * (a * n_y * lam) / (n_y * lam + tY)

    # This if statement is used to catch a float overflow warning.
    # sometimes elements of (nY*lam*tY) can be infinite, but this will reduce the second term to zero.
    # Analogously, the _a_ value can be set to zero, leading to an overflow error as well, but setting the second term
    # to zero as well.
    Dyy = (1. / (n_y * lam + tY))
    if a != 0:  # catching a overflow error
        selection = (n_y * lam + tY) < np.sqrt(np.finfo(np.float64).max)
        if np.any(selection):
            Dyy[selection] = Dyy[selection] + (Dxx * (a ** 2 * n_y ** 2 * lam ** 2)) / ((n_y * lam + tY) ** 2)

    dX = n_x * c_x + a * n_y * c_y
    dY = n_y * c_y
    m = len(c_x)

    loglik = -m * np.log(2 * np.pi) + \
             -(1 / 2) * sum(
        np.log((a ** 2 * n_y + n_x) * lam + tX - a ** 2 * n_y ** 2 * (lam ** 2) / (n_y * lam + tY))) + \
             -(1 / 2) * sum(np.log(n_y * lam + tY)) + \
             +(1 / 2) * (sum(dX ** 2 * Dxx) + 2 * sum(dX * dY * Dxy) + sum(dY ** 2 * Dyy)) + \
             -(n_x / 2) * sum((c_x ** 2) / lam) + \
             -(n_y / 2) * sum((c_y ** 2) / lam) + \
             +(m / 2) * (np.log(n_x) + np.log(n_y)) - sum(np.log(lam)) + (m / 2) * (np.log(tX) + np.log(tY))

    return -loglik


def mr_link2_loglik_reference_v2(th: np.ndarray, lam: np.ndarray,
                                 c_x: np.ndarray, c_y: np.ndarray,
                                 n_x: float, n_y: float) -> float:
    """
    The MR-link2 log likelihood function. This function calculates -1 * likelihood of three parameters:
    alpha, sigma_x and sigma_y.
    Designed to be used in optimization algorithms like those in scipy.minimize

    :param th:
        List or numpy array of floats with the parameters to optimize first is alpha, second the
        1 /  exposure heritability and third the 1/ outcome heritability.
    :param lam:
        np.ndarray of selected eigenvalues of the cX and cY parameters.
    :param c_x:
        The dot product of the selected eigenvectors and summary statistics vector of the exposure
    :param c_y:
        The dot product of the selected eigenvectors and summary statistics vector of the outcome
    :param n_x:
        The number of individuals in the exposure dataset
    :param n_y:
        The number of individuals in the outcome dataset
    :return:
        a single float that contains the likelihood of the parameters theta.
    """

    n_x = float(n_x)
    n_y = float(n_y)
    a = th[0]
    tX = abs(th[1])
    tY = abs(th[2])

    Dyy = (1. / (n_y * lam + tY))

    if a != 0.0:
        Dxx = 1. / (np.exp(np.log(a ** 2 * n_y + n_x) + np.log(lam)) + tX -
                    np.exp(np.log(a ** 2 * n_y ** 2 * (lam ** 2)) - np.log(n_y * lam + tY)))
        Dxy = -Dxx * a * np.exp(np.log((n_y * lam)) - np.log(n_y * lam + tY))
        Dyy = Dyy + np.exp(np.log(Dxx * (a ** 2 * n_y ** 2 * lam ** 2)) - (2 * np.log(n_y * lam + tY)))
        asq_ny_sq_lam_sq_div_ny_lam_ty = np.exp(np.log(a ** 2 * n_y ** 2 * (lam ** 2)) - np.log(n_y * lam + tY))
    else:
        Dxx = 1. / (np.exp(np.log(n_x) + np.log(lam)) + tX)
        Dxy = -Dxx * a * np.exp(np.log((n_y * lam)) - np.log(n_y * lam + tY))
        Dyy = Dyy
        asq_ny_sq_lam_sq_div_ny_lam_ty = 0.0 * lam

    dX = n_x * c_x + a * n_y * c_y
    dY = n_y * c_y
    m = len(c_x)

    loglik = -m * np.log(2 * np.pi) + \
             -(1 / 2) * sum(np.log((a ** 2 * n_y + n_x) * lam + tX - asq_ny_sq_lam_sq_div_ny_lam_ty)) + \
             -(1 / 2) * sum(np.log(n_y * lam + tY)) + \
             +(1 / 2) * (sum(dX ** 2 * Dxx) + 2 * sum(dX * dY * Dxy) + sum(dY ** 2 * Dyy)) + \
             -(n_x / 2) * sum((c_x ** 2) / lam) + \
             -(n_y / 2) * sum((c_y ** 2) / lam) + \
             +(m / 2) * (np.log(n_x) + np.log(n_y)) - sum(np.log(lam)) + (m / 2) * (np.log(tX) + np.log(tY))

    return -loglik


def mr_link2_loglik_alpha_h0(th, lam, cX, cY, nX, nY) -> float:  # fix alpha to zero
    """
    The MR-link2 log likelihood function when the causal effect is zero.
    This function calculates -1 * likelihood of two parameters:
    sigma_x and sigma_y. the causal effect alpha is set to zero
    Designed to be used in optimization algorithms like those in scipy.minimize

    :param th:
        List or numpy array of floats with the 2 parameters to optimize first is alpha, second the exposure heritability
        and third the outcome heritability.
    :param lam:
        np.ndarray of selected eigenvalues of the cX and cY parameters.
    :param cX:
        The dot product of the selected eigenvectors and summary statistics vector of the exposure
    :param cY:
        The dot product of the selected eigenvectors and summary statistics vector of the outcome
    :param nX:
        The number of individuals in the exposure dataset
    :param nY:
        The number of individuals in the outcome dataset
    :return:
        a single float that contains the likelihood of the parameters theta.
    """

    return mr_link2_loglik_reference_v2(np.asarray([0.0, th[0], th[1]]), lam, cX, cY, nX, nY)


def mr_link2_loglik_sigma_y_h0(th, lam, c_x, c_y, n_x, n_y) -> float:  # this is for fixing sigma_y to zero
    """
    The MR-link2 log likelihood function when the pleiotropic effect is zero.
    This function calculates -1 * likelihood of two parameters:
    alpha and sigma_x. the causal effect alpha is set to zero
    Designed to be used in optimization algorithms like those in scipy.minimize

    :param th:
        List or numpy array of floats with the 2 parameters to optimize first is alpha, second the exposure heritability

    :param lam:
        np.ndarray of selected eigenvalues of the cX and cY parameters.
    :param c_x:
        The dot product of the selected eigenvectors and summary statistics vector of the exposure
    :param c_y:
        The dot product of the selected eigenvectors and summary statistics vector of the outcome
    :param n_x:
        The number of individuals in the exposure dataset
    :param n_y:
        The number of individuals in the outcome dataset
    :return:
        a single float that contains the likelihood of the parameters theta.
    """

    n_x = float(n_x)
    n_y = float(n_y)
    a = th[0]
    tX = abs(th[1])

    Dyy = np.zeros_like(lam)

    if a != 0.0:
        Dxx = 1. / (np.exp(np.log(a ** 2 * n_y + n_x) + np.log(lam)) + tX)
        Dxy = np.zeros_like(lam)  # -Dxx * a * np.exp(np.log((n_y * lam)) - np.log(n_y * lam + tY))
        Dyy = Dyy + np.zeros_like(
            lam)  # np.exp(np.log(Dxx * (a ** 2 * n_y ** 2 * lam ** 2)) - (2 * np.log(n_y * lam + tY)))
        asq_ny_sq_lam_sq_div_ny_lam_ty = np.zeros_like(
            lam)  # np.exp(np.log(a ** 2 * n_y ** 2 * (lam ** 2)) - np.log(n_y * lam + tY))
    else:
        Dxx = 1. / (np.exp(np.log(n_x) + np.log(lam)) + tX)
        Dxy = np.zeros_like(lam)  # -Dxx * a * np.exp(np.log((n_y * lam)) - np.log(n_y * lam + tY))
        Dyy = Dyy
        asq_ny_sq_lam_sq_div_ny_lam_ty = np.zeros_like(lam)

    dX = n_x * c_x + a * n_y * c_y
    dY = n_y * c_y
    m = len(c_x)

    loglik = -m * np.log(2 * np.pi) + \
             -(1 / 2) * sum(np.log((a ** 2 * n_y + n_x) * lam + tX - asq_ny_sq_lam_sq_div_ny_lam_ty)) + \
             +(1 / 2) * (sum(dX ** 2 * Dxx) + 2 * sum(dX * dY * Dxy) + sum(dY ** 2 * Dyy)) + \
             -(n_x / 2) * sum((c_x ** 2) / lam) + \
             -(n_y / 2) * sum((c_y ** 2) / lam) + \
             +(m / 2) * (np.log(n_x) + np.log(n_y)) - sum(np.log(lam)) + (m / 2) * (np.log(tX))
    # + (m/2) * np.log(tY)) -(1 / 2) * sum(np.log(n_y * lam + tY)) ## These two terms should cancel out.

    return -loglik


def mr_link2(selected_eigenvalues: np.ndarray, selected_eigenvectors: np.ndarray,
             exposure_betas: np.ndarray, outcome_betas: np.ndarray, n_exp: float, n_out: float,
             sigma_exp_guess: float, sigma_out_guess: float) -> dict[str, Any]:
    """
    Run MR-link-2, perform two likelihood ratio tests: one for the pleiotropic effect sigma_y and one for the
    causal effect alpha

    :param selected_eigenvalues: eigenvalues that are selected for this MR-link2 run
    :param selected_eigenvectors: Eigenvectors that are selected for this MR-link2 run
    :param exposure_betas: a vector of exposure betas
    :param outcome_betas: a vector of outcome betas
    :param n_exp: number of individuals in the exposure cohort
    :param n_out: number of individuals in the outcome cohort
    :param sigma_exp_guess: guess for sigma exposure
    :param sigma_out_guess: guess for sigma outcoem

    :return: returns a dictionary containing the results of the MR-link-2 optimizations
             that can be used as a row in a pandas dataframe
    """
    start_time = time.time()
    method = 'Nelder-Mead'
    options = {'maxiter': 300, 'disp': False}

    c_x = selected_eigenvectors.T @ exposure_betas
    c_y = selected_eigenvectors.T @ outcome_betas

    max_sigma = np.sqrt(np.finfo(np.float64).max)

    """ 
    Alpha h0 estimation
    """
    alpha_h0_guesses = [
        [sigma_exp_guess, sigma_out_guess],
        [max_sigma, max_sigma],
        [1, max_sigma],
        [1e3, 1e3]
    ]

    alpha_h0_results = scipy.optimize.minimize(mr_link2_loglik_alpha_h0, args=(selected_eigenvalues,
                                                                               c_x,
                                                                               c_y,
                                                                               n_exp,
                                                                               n_out),
                                               x0=np.asarray(alpha_h0_guesses[0]),
                                               method=method, options=options)
    for alpha_h0_guess in alpha_h0_guesses[1:]:  # first one already done

        if alpha_h0_results.success:
            break

        new_alpha_h0_results = scipy.optimize.minimize(mr_link2_loglik_alpha_h0, args=(selected_eigenvalues,
                                                                                       c_x,
                                                                                       c_y,
                                                                                       n_exp,
                                                                                       n_out),
                                                       x0=np.asarray(alpha_h0_guess),
                                                       method=method, options=options)
        if alpha_h0_results.fun >= new_alpha_h0_results.fun:
            alpha_h0_results = new_alpha_h0_results

    """
    Sigma_y estimation
    """
    sigma_y_guesses = [[0.0, sigma_exp_guess],
                       [1.0, sigma_exp_guess],
                       [0.0, alpha_h0_results.x[0]],
                       [1.0, alpha_h0_results.x[0]],
                       [0.0, max_sigma],
                       [1e-10, max_sigma]
                       ]

    sigma_y_h0_results = scipy.optimize.minimize(mr_link2_loglik_sigma_y_h0, args=(selected_eigenvalues,
                                                                                   c_x,
                                                                                   c_y,
                                                                                   n_exp,
                                                                                   n_out),
                                                 x0=np.asarray(sigma_y_guesses[0]),
                                                 method=method, options=options)

    for sigma_y_guess in sigma_y_guesses[1:]:
        if sigma_y_h0_results.success:
            break

        new_sigma_y_h0_results = scipy.optimize.minimize(mr_link2_loglik_sigma_y_h0, args=(selected_eigenvalues,
                                                                                           c_x,
                                                                                           c_y,
                                                                                           n_exp,
                                                                                           n_out),
                                                         x0=np.asarray(sigma_y_guess),
                                                         method=method, options=options)
        if new_sigma_y_h0_results.fun < sigma_y_h0_results.fun:
            sigma_y_h0_results = new_sigma_y_h0_results

    """
    Ha estimation
    """
    ha_guesses = [
        [0.0, alpha_h0_results.x[0], alpha_h0_results.x[1]],
        [sigma_y_h0_results.x[0], sigma_y_h0_results.x[1], np.sqrt(np.finfo(np.float64).max)],
        [1.0, alpha_h0_results.x[0], alpha_h0_results.x[1]],
        [1e-10, max_sigma, max_sigma]
    ]
    ha_results = scipy.optimize.minimize(fun=mr_link2_loglik_reference_v2, args=(selected_eigenvalues,
                                                                                 c_x,
                                                                                 c_y,
                                                                                 n_exp,
                                                                                 n_out),
                                         x0=np.asarray(ha_guesses[0], dtype=float),
                                         method=method, options=options)

    for ha_guess in ha_guesses[1:]:
        if ha_results.success:
            break
        new_ha_result = scipy.optimize.minimize(fun=mr_link2_loglik_reference_v2,
                                                args=(selected_eigenvalues,
                                                      c_x,
                                                      c_y,
                                                      n_exp,
                                                      n_out),
                                                x0=np.asarray(ha_guess, dtype=float),
                                                method=method, options=options)
        if new_ha_result.fun < ha_results.fun:
            ha_results = new_ha_result

    #
    if True or not ha_results.success:
        a = mr_link2_loglik_reference_v0(ha_results.x, selected_eigenvalues, c_x, c_y, n_exp, n_out)
        b = mr_link2_loglik_reference_v2(ha_results.x, selected_eigenvalues, c_x, c_y, n_exp, n_out)

        if not np.isclose(a, b):
            print(f'Functions are not the same {a} compared to {b}, difference: {a - b:.4e}')
            pass  # for debugging purposes

    ## Now take the likelihoods and start doing the likelihood ratio test
    alpha = ha_results.x[0]
    alpha_chi_sq = 2 * (
            alpha_h0_results.fun - ha_results.fun)  ## This is the likelihood ratio, because they are in log scale
    alpha_p_val = scipy.stats.chi2.sf(alpha_chi_sq, 1)  ## This is
    z_alpha = 0.0 if alpha_chi_sq <= 0 else np.sign(alpha) * np.sqrt(alpha_chi_sq)
    se_alpha = alpha / z_alpha if z_alpha != 0 else np.nan

    sigma_y = 1 / abs(ha_results.x[2])
    sigma_y_chi_sq = 2 * (sigma_y_h0_results.fun - ha_results.fun)
    sigma_y_p_val = scipy.stats.chi2.sf(sigma_y_chi_sq, 1)
    z_sigma_y = 0.0 if sigma_y_chi_sq <= 0 else np.sqrt(sigma_y_chi_sq)
    se_sigma_y = sigma_y / z_sigma_y if z_sigma_y != 0 else np.nan

    to_return = {'alpha': alpha, 'se(alpha)': se_alpha, 'p(alpha)': alpha_p_val,
                 'sigma_y': sigma_y, 'se(sigma_y)': se_sigma_y, 'p(sigma_y)': sigma_y_p_val,
                 'sigma_x': 1 / abs(ha_results.x[1]),
                 'alpha_h0_sigma_x': 1 / abs(alpha_h0_results.x[0]),
                 'alpha_h0_sigma_y': 1 / abs(alpha_h0_results.x[1]),
                 'alpha_h0_loglik': alpha_h0_results.fun,
                 'sigma_y_h0_alpha': sigma_y_h0_results.x[0],
                 'sigma_y_h0_sigma_x': 1 / abs(sigma_y_h0_results.x[1]),
                 'sigma_y_h0_loglik': sigma_y_h0_results.fun,
                 'ha_loglik': ha_results.fun,
                 'optim_alpha_h0_success': alpha_h0_results.success,
                 'optim_alpha_h0_nit': alpha_h0_results.nit,
                 'optim_sigma_y_h0_success': sigma_y_h0_results.success,
                 'optim_sigma_y_h0_nit': sigma_y_h0_results.nit,
                 'optim_ha_success': ha_results.success,
                 'optim_ha_nit': ha_results.nit,
                 'function_time': time.time() - start_time
                 }

    return to_return


def select_instruments_by_clumping(p_values: np.ndarray,
                                   correlation_matrix: np.ndarray,
                                   clumping_p_threshold: float,
                                   r_threshold=0.1) -> list[int]:
    """
    This function performs clumping based on a pre-computed correlation matrix.

    :param p_values: np.ndarry of p values
    :param correlation_matrix:  a square correlation matrix **not squared correlation**
    :param clumping_p_threshold: the p value threshold to select the SNPs on
    :param r_threshold: the correlation threshold in pearson correlation. **not squared**
    :return: a list of indices as selected SNPS.
    """
    selection = []
    mask = np.ones(len(p_values), dtype=bool)

    while min(p_values[mask]) <= clumping_p_threshold:
        mask_idxs = sorted(np.where(mask)[0])
        selected = mask_idxs[np.argmin(p_values[mask])]
        selection.append(selected)

        to_keep = np.abs(correlation_matrix[selected, :]) <= r_threshold
        mask = np.logical_and(mask, to_keep).reshape((len(p_values)))
        if np.sum(mask) == 0:
            break

    return selection


def match_n_variants_n_individuals_from_plink_log(plink_logfile: str) -> Tuple[int, int]:
    """
    This is a utitlity function to
    :param plink_logfile: str
        location of the plink logfile which will be parsed for the number of individuals and variants
    :return: tuple
        Tuple containing the numbeer of variants and number of individuals
    """

    regex = re.compile('^([0-9]+) variants and ([0-9]+) people pass filters and QC\.')
    with open(plink_logfile) as f:
        for line in f:
            match_obj = regex.match(line)
            if match_obj:
                return int(match_obj[1]), int(match_obj[2])

    raise ValueError('Did not find n_individuals in the plink log.')



def run_colocalization_analysis(bXs, bYs, pXs, pYs, n_exp, n_out, mafs, ld_mat: np.ndarray,
                          out_file='tmp__colocalization__testing_out.txt.gz', in_file='tmp__colocalization__testing_in.txt.gz', ld_file='tmp__colocalization__testing_ld.bin'):
    stderr = subprocess.DEVNULL if not verbosity else None
    stdout = subprocess.DEVNULL if not verbosity else None

    all_dfs = []
    # ['MAF', 'pvalues', 'N', 'type', 'simulation'])
    if len(bXs.shape) == 2:
        n_sims = bXs.shape[1]
        for i_pheno in range(n_sims):
            exp_tmp_df = pd.DataFrame()
            exp_tmp_df['beta'] = bXs[:, i_pheno] * mafs * (1 - mafs)
            exp_tmp_df['MAF'] = mafs
            exp_tmp_df['p_values'] = pXs[:, i_pheno]
            exp_tmp_df['N'] = n_exp
            exp_tmp_df['type'] = 'quant'
            exp_tmp_df['simulation'] = i_pheno
            exp_tmp_df['exp_or_outcome'] = 'exposure'

            out_tmp_df = pd.DataFrame()
            out_tmp_df['beta'] = bYs[:, i_pheno] * mafs * (1 - mafs)
            out_tmp_df['MAF'] = mafs
            out_tmp_df['p_values'] = pYs[:, i_pheno]
            out_tmp_df['N'] = n_out
            out_tmp_df['type'] = 'quant'
            out_tmp_df['simulation'] = i_pheno
            out_tmp_df['exp_or_outcome'] = 'outcome'

            all_dfs.append(exp_tmp_df)
            all_dfs.append(out_tmp_df)

    else:
        n_sims = 1
        for i_pheno in range(n_sims):
            exp_tmp_df = pd.DataFrame()
            exp_tmp_df['beta'] = bXs * mafs * (1 - mafs)
            exp_tmp_df['MAF'] = mafs
            exp_tmp_df['p_values'] = pXs
            exp_tmp_df['N'] = n_exp
            exp_tmp_df['type'] = 'quant'
            exp_tmp_df['simulation'] = i_pheno
            exp_tmp_df['exp_or_outcome'] = 'exposure'

            out_tmp_df = pd.DataFrame()
            out_tmp_df['beta'] = bYs * mafs * (1 - mafs)
            out_tmp_df['MAF'] = mafs
            out_tmp_df['p_values'] = pYs
            out_tmp_df['N'] = n_out
            out_tmp_df['type'] = 'quant'
            out_tmp_df['simulation'] = i_pheno
            out_tmp_df['exp_or_outcome'] = 'outcome'

            all_dfs.append(exp_tmp_df)
            all_dfs.append(out_tmp_df)

    full_df = pd.concat(all_dfs)
    full_df.to_csv(out_file, sep='\t', index=False)
    ld_mat.tofile(ld_file)

    subprocess.run(['Rscript', f'{os.path.dirname(os.path.abspath(__file__))}/r_analyses/colocalization_analysis.R',
                    out_file, in_file, ld_file],
                   check=True,
                   stdout=stdout,
                   stderr=stderr,
                   )

    if os.path.exists(in_file):
        results = pd.read_csv(in_file, sep='\t')
    else:
        raise ValueError('Could not perform coloc')

    os.remove(out_file)
    os.remove(in_file)
    os.remove(ld_file)

    return results


def run_external_mr_analysis(bXs, bYs, seXs, seYs, pXs, pYs, n_exp, n_out, mafs, ld_mat: np.ndarray, instruments,
                             out_file='tmp__colocalization_mr_testing_out.txt.gz', in_file='mr_testing_in.txt.gz',
                             ld_file='mr_testing_ld.bin'):

    """
    This function is written to run the MR analyses

    :param bXs:
    :param bYs:
    :param seXs:
    :param seYs:
    :param pXs:
    :param pYs:
    :param n_exp:
    :param n_out:
    :param mafs:
    :param ld_mat:
    :param instruments:
    :param out_file:
    :param in_file:
    :param ld_file:
    :return: a dataframe of results.
    """
    all_dfs = []

    print(instruments.shape)
    if len(bXs.shape) == 2:
        for i_pheno in range(bXs.shape[1]):
            exp_tmp_df = pd.DataFrame()
            exp_tmp_df['beta'] = bXs[:, i_pheno]  # * mafs * (1 - mafs)
            exp_tmp_df['ses'] = seXs[:, i_pheno]  # * mafs * (1 - mafs)
            exp_tmp_df['MAF'] = mafs
            exp_tmp_df['p_values'] = pXs[:, i_pheno]
            exp_tmp_df['N'] = n_exp
            exp_tmp_df['type'] = 'quant'
            exp_tmp_df['simulation'] = i_pheno
            exp_tmp_df['exp_or_outcome'] = 'exposure'
            exp_tmp_df['selected_as_instrument'] = instruments[:, i_pheno]

            out_tmp_df = pd.DataFrame()
            out_tmp_df['beta'] = bYs[:, i_pheno]  # * mafs * (1 - mafs)
            out_tmp_df['ses'] = seYs[:, i_pheno]
            out_tmp_df['MAF'] = mafs
            out_tmp_df['p_values'] = pYs[:, i_pheno]
            out_tmp_df['N'] = n_out
            out_tmp_df['type'] = 'quant'
            out_tmp_df['simulation'] = i_pheno
            out_tmp_df['exp_or_outcome'] = 'outcome'
            out_tmp_df['selected_as_instrument'] = instruments[:, i_pheno]

            all_dfs.append(exp_tmp_df)
            all_dfs.append(out_tmp_df)

    else:
        for i_pheno in range(1):
            exp_tmp_df = pd.DataFrame()
            exp_tmp_df['beta'] = bXs
            exp_tmp_df['ses'] = seXs
            exp_tmp_df['MAF'] = mafs
            exp_tmp_df['p_values'] = pXs
            exp_tmp_df['N'] = n_exp
            exp_tmp_df['type'] = 'quant'
            exp_tmp_df['simulation'] = i_pheno
            exp_tmp_df['exp_or_outcome'] = 'exposure'
            exp_tmp_df['selected_as_instrument'] = instruments

            out_tmp_df = pd.DataFrame()
            out_tmp_df['beta'] = bYs
            out_tmp_df['ses'] = seYs
            out_tmp_df['MAF'] = mafs
            out_tmp_df['p_values'] = pYs
            out_tmp_df['N'] = n_out
            out_tmp_df['type'] = 'quant'
            out_tmp_df['simulation'] = i_pheno
            out_tmp_df['exp_or_outcome'] = 'outcome'
            out_tmp_df['selected_as_instrument'] = instruments

            all_dfs.append(exp_tmp_df)
            all_dfs.append(out_tmp_df)

    full_df = pd.concat(all_dfs)
    full_df.to_csv(out_file, sep='\t', index=False, float_format='%.30e')
    ld_mat.tofile(ld_file)

    subprocess.run(['Rscript', f'{os.path.dirname(os.path.abspath(__file__))}/r_analyses/mr_analyses.R',
                    out_file, in_file, ld_file],
                   check=True  # , stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
                   )

    if os.path.exists(in_file):
        results = pd.read_csv(in_file, sep='\t')
    else:
        raise ValueError('Could not perform external mr analyses')

    return results


def mr_link2_on_region(region: StartEndRegion,
                       exposure_df: pd.DataFrame,
                       outcome_df: pd.DataFrame,
                       bed_file: str,
                       maf_threshold: float,
                       max_correlation: float,
                       regional_ld_matrix: np.ndarray,
                       snps_alleles_in_ld_matrix: list[tuple],
                       tmp_prepend: str,
                       var_explained_grid: list,
                       verbosity=0,
                       max_snp_threshold=5250,
                       ) -> pd.DataFrame:
    """

    This function takes summary statistics files, and a bed file and some parameters, and outputs MR-link 2 estimates.

    :param region: StartEndRegion
        a StartEndRegion object that is used to prune the exposure, outcome and bed file.
    :param exposure_df:
        pandas DataFrame from which the exosure summary statistics are used
    :param outcome_df:  pd.DataFrame
        pandas DataFrame from which the outcome summary statistics are used
    :param bed_file: str
        Prepend for the bed file from which the correlation matrix is determined.
    :param maf_threshold: float
        minor allele frequency threshold for selecting SNPs for the correlation matrix
    :param max_correlation: float
        The maximum allowed correlation value of a pair of variants. We use this value to remove any snps that are too
        close together
    :param tmp_prepend: str
        A prepend of where to store temporary files, if the function does not error out, these files will be cleaned up
    :param var_explained_grid: list or iterable containing floats
        with variance explained values for which to perform MR-link 2
    :param verbosity: int,
        A debugging parameter that  can be used to print diagnostic information
    :param max_snp_threshold: int,
        The number of SNPs that are allowed to be present in the region before

    :return: pd.DataFrame
        pd.DataFrame with the results of a row the MR-link 2
        inference for each variance explained in the parameter space

    """
    ##setup
    stderr = subprocess.DEVNULL if not verbosity else None
    stdout = subprocess.DEVNULL if not verbosity else None

    ## Isolate regions of interest.
    regional_exp_df = exposure_df[
        (exposure_df.chromosome.astype(str) == region.chromosome) &
        (exposure_df.base_pair_location.astype(int) >= region.start) &
        (exposure_df.base_pair_location.astype(int) <= region.end)
        ]

    regional_out_df = outcome_df[
        (outcome_df.chromosome.astype(str) == region.chromosome) &
        (outcome_df.base_pair_location >= region.start) &
        (outcome_df.base_pair_location <= region.end)
        ]
    if verbosity:
        print(regional_exp_df.shape)
        print(regional_out_df.shape)

    ld_matrix_snp_names = {x[0] for x in snps_alleles_in_ld_matrix}
    ld_matrix_snp_to_allele_dict = {x[0]: (x[1], x[2]) for x in snps_alleles_in_ld_matrix}
    shared_snps = set(regional_exp_df.rsid.tolist()) & set(regional_out_df.rsid.tolist()) & ld_matrix_snp_names

    if verbosity:
        print(f'Found {len(shared_snps)} shared SNPS')
    if len(shared_snps) == 0:
        print('Could not make an MR estimate due to no shared SNPs, returning None')
        return None

    regional_exp_df = regional_exp_df[regional_exp_df.rsid.isin(shared_snps)]
    regional_out_df = regional_out_df[regional_out_df.rsid.isin(shared_snps)]

    ## make sure that the LD matrix contains only the shared SNPs.
    indexes_to_keep = np.asarray([i for i in range(len(snps_alleles_in_ld_matrix)) if snps_alleles_in_ld_matrix[i][0] in shared_snps], dtype=int)
    ordered_snps = [snps_alleles_in_ld_matrix[i][0] for i in range(len(snps_alleles_in_ld_matrix)) if snps_alleles_in_ld_matrix[i][0] in shared_snps]
    correlation_mat = copy.copy(regional_ld_matrix[indexes_to_keep, :][:, indexes_to_keep])

    ## harmonize the betas, exposure is reference
    snp_to_beta_dict = {}
    exp_dict = dict(zip(regional_exp_df['rsid'],
                        zip(regional_exp_df['beta'],
                            regional_exp_df['effect_allele'],
                            regional_exp_df['other_allele'],
                            regional_exp_df['p_value'],
                            regional_exp_df['standard_error'])))

    out_dict = dict(zip(regional_out_df['rsid'],
                        zip(regional_out_df['beta'],
                            regional_out_df['effect_allele'],
                            regional_out_df['other_allele'],
                            regional_out_df['p_value'],
                            regional_out_df['standard_error'])))

    """ 
    harmonize the summary statistics 
    to the LD reference alleles.  
    """
    for rsid in shared_snps:
        ## Use the LD matrix alleles as the underlying reference
        ref_a1, ref_a2 = ld_matrix_snp_to_allele_dict[rsid]
        exp_beta, exp_a1, exp_a2, _, _ = exp_dict[rsid]
        out_beta, out_a1, out_a2, _, _ = out_dict[rsid]

        # test if the alleles are the same, if not we should not continue
        if {exp_a1, exp_a2} != {ref_a1, ref_a2} or {out_a1, out_a2} != {ref_a1, ref_a2}:
            raise ValueError(f'Harmonization error in snp {rsid}, ld mat: {(ref_a1, ref_a2)}, exposure: {(exp_a1, out_a2)} outcome: {(out_a1, out_a2)}')

        # as we tested that the alelles are the same set, only have to look at alignment.
        if exp_a1 != ref_a1:
            exp_beta = -1 * exp_beta
        if out_a1 != ref_a1:
            out_beta = -1 * out_beta

        snp_to_beta_dict[rsid] = [exp_beta, out_beta]

    ## these should have been normalized. In my case they are.
    exp_betas = np.asarray([snp_to_beta_dict[x][0] for x in ordered_snps], dtype=float)
    exp_pvals = np.asarray([exp_dict[x][3] for x in ordered_snps], dtype=float)
    exp_ses = np.asarray([exp_dict[x][4] for x in ordered_snps], dtype=float)

    out_betas = np.asarray([snp_to_beta_dict[x][1] for x in ordered_snps], dtype=float)
    out_pvals = np.asarray([out_dict[x][3] for x in ordered_snps], dtype=float)
    out_ses = np.asarray([out_dict[x][4] for x in ordered_snps], dtype=float)


    """
    Now, the optimization
    """
    m_snps = exp_betas.shape[0]

    print(f'We have {m_snps} snps in the {region} region')
    if m_snps > max_snp_threshold:
        raise ValueError(
            f'Too many SNPs in {region} region, after filtering, contains: {m_snps} snps > {max_snp_threshold}'
            'will take too long. Consider reducing the region size or pruning the correlation matrix more.')

    print('    Starting external analyses')
    n_exp = int(np.median(regional_exp_df.n))
    n_out = int(np.median(regional_out_df.n))

    mafs = np.zeros_like(exp_betas)
    mafs[:] = 0.5


    coloc_results = run_colocalization_analysis(exp_betas, out_betas, exp_pvals, out_pvals, n_exp, n_out, mafs,
                                          correlation_mat, out_file = tmp_prepend + '_coloc_out',
                                          in_file = tmp_prepend + '_coloc_in', ld_file = tmp_prepend + 'coloc_ld.bin')
    instruments = np.zeros(exp_betas.shape[0], dtype=bool)
    tmp_instruments = np.asarray(select_instruments_by_clumping(exp_pvals, correlation_mat, 5e-8, r_threshold=0.1), dtype=int)
    for instrument in tmp_instruments:
        instruments[instrument] = True
    r_mr_results = run_external_mr_analysis(exp_betas, out_betas, exp_ses, out_ses,  exp_pvals, out_pvals, n_exp, n_out, mafs,
                                          correlation_mat,instruments= instruments, out_file = tmp_prepend + '_mr_out',
                                          in_file = tmp_prepend + '_mr_in', ld_file = tmp_prepend + 'mr_ld.bin')

    print('     Starting MR-link2')

    pseudo_exp_h2 = np.sum(exp_betas ** 2 - (1 / n_exp))
    pseudo_out_h2 = np.sum(out_betas ** 2 - (1 / n_out))

    exposure_instruments = select_instruments_by_clumping(exp_pvals,
                                                          correlation_mat,
                                                          clumping_p_threshold=5e-6,
                                                          r_threshold=0.1)  # r is pearson
    pseudo_instrument_exp_h2 = np.sum(exp_betas[exposure_instruments] ** 2 - (1 / n_exp))

    outcome_instruments = select_instruments_by_clumping(out_pvals,
                                                         correlation_mat,
                                                         clumping_p_threshold=5e-6,
                                                         r_threshold=0.1)  # r is pearson
    pseudo_instrument_out_h2 = np.sum(exp_betas[outcome_instruments] ** 2 - (1 / n_exp))

    if np.isclose(pseudo_instrument_out_h2, 0.0):
        out_h2_guess = pseudo_out_h2
    else:
        out_h2_guess = pseudo_instrument_out_h2

    if np.isclose(pseudo_instrument_exp_h2, 0.0):
        exp_h2_guess = pseudo_exp_h2
    else:
        exp_h2_guess = pseudo_instrument_exp_h2

    sigma_guesses = np.zeros((2), dtype=float)
    sigma_guesses[:] = float(m_snps / exp_h2_guess), float(m_snps / out_h2_guess)

    lam, u = np.linalg.eigh(correlation_mat)

    if verbosity:
        print(f'original lambda {lam.shape}')
        print(lam[:10:-1])
        print(f'original u  {u.shape}')
        print(u[:10:-1, :10:-1])

    # var_explained_grid = np.arange(0.7, 1.0, 0.01)
    results_list = []

    for i, var_explained in enumerate(var_explained_grid):
        mr_link2_estimates = []
        variances_explained = np.cumsum(sorted(lam[lam > 0.0], reverse=True)) / sum(lam[lam > 0.0])
        threshold = sorted(lam[lam > 0.0], reverse=True)[np.argmin(variances_explained <= var_explained)]

        correlation_variance_explained = sum(lam[lam > threshold]) / sum(lam)
        n_components_selected = np.sum(lam > threshold)

        if verbosity:
            print(f'Selected {n_components_selected} / {m_snps} components at threshold {threshold}, '
                  f'representing {correlation_variance_explained * 100}% of variance explained')
        selection = lam >= threshold

        mr_link2_point_estimate = mr_link2(lam[selection], u[:, selection],
                                           exp_betas, out_betas,
                                           n_exp, n_out, sigma_guesses[0], sigma_guesses[1])
        if verbosity:
            print(mr_link2_point_estimate)

        mr_link2_point_estimate = pd.DataFrame(mr_link2_point_estimate, index=[0])
        mr_link2_estimates.append(mr_link2_point_estimate)

        this_run = pd.concat(mr_link2_estimates)

        this_run['var_explained'] = var_explained

        this_run = this_run.join(r_mr_results, rsuffix='')
        this_run = this_run.join(coloc_results, rsuffix='')

        results_list.append(this_run)

    mr_results_df = pd.concat(results_list)

    mr_results_df['region'] = str(region)
    mr_results_df['m_snps_overlap'] = m_snps


    columns = ['region', 'var_explained', 'm_snps_overlap',
                                   'alpha', 'se(alpha)', 'p(alpha)',
                                   'sigma_y', 'se(sigma_y)', 'p(sigma_y)',
                                   'sigma_x', 'function_time',
                                   'beta_ivw', 'se_ivw', 'p_ivw',
                                   'beta_ivw_r', 'se_ivw_r', 'p_ivw_r',
                                   'beta_pca', 'se_pca', 'p_pca',
                                   'PP.H1.abf', 'PP.H2.abf', 'PP.H3.abf', 'PP.H4.abf']

    if 'susie_max_PP.H4.abf' in mr_results_df.columns:
        columns += ['susie_max_PP.H1.abf',
                    'susie_max_PP.H2.abf',
                    'susie_max_PP.H3.abf',
                    'susie_max_PP.H4.abf']


    mr_results_df = mr_results_df[columns]

    ## this is a normalization step that is done as the likelihood function computes them per SNP
    mr_results_df['sigma_x'] = mr_results_df['sigma_x'] * mr_results_df['m_snps_overlap']

    mr_results_df['sigma_y'] = mr_results_df['sigma_y'] * mr_results_df['m_snps_overlap']
    mr_results_df['se(sigma_y)'] = mr_results_df['se(sigma_y)'] * mr_results_df['m_snps_overlap']

    ## clean up
    for filename in files_to_remove:
        if os.path.exists(filename):
            os.remove(filename)

    return mr_results_df


def read_ld_matrix_local(geno_file_obj: PlinkGenoReader, list_of_snps, maf_threshold, max_correlation):
    """
    return a correlation matrix of a list of snps.

    This function has been tested to produce the same results as the plink read ld matrix function, very happy about this!

    :param geno_file_obj:
    :param list_of_snps:
    :param maf_threshold:
    :param max_correlation:
    :return:
    """

    def maf_func(x):
        return np.sum(x) / (x.shape[0] * 2)

    overlapping_snps = set(list_of_snps) & geno_file_obj.snp_name_to_variant_order.keys()
    if len(overlapping_snps) == 0:
        raise ValueError(f'No overlapping SNPs between the LD reference and the region {list_of_snps[:10]=}, '
                         f'{list(geno_file_obj.snp_name_to_variant_order)[:10]=}')

    """
    Sorting the overlapping SNPs based on their position, to ensure correctness between the plink implementation
    """

    overlapping_snps = [x[0] for x in sorted([(snp_name, position) for snp_name, position in
                                              zip(overlapping_snps, [geno_file_obj.bim_data[x] for x in overlapping_snps])],
                                             key=lambda x: int(x[1][1]))]

    genotypes, all_snps = geno_file_obj.read_list_of_snps_into_geno_array(overlapping_snps)
    #drop genotypes with missing values.
    missing_values = np.sum(genotypes == 3, axis=0) >= 1

    # Calculate minor allele frequencies (MAFs)
    mafs = np.apply_along_axis(maf_func, 0, genotypes)
    maf_filter = np.logical_or(mafs <= maf_threshold,
                                mafs >= 1-maf_threshold)

    genotypes = genotypes[:, np.logical_and(~missing_values, ~maf_filter)]
    genotype_matrix_snps = [x for x, bool_1, bool_2, in zip(all_snps, missing_values, maf_filter)
                            if not (bool_1 or bool_2)]
    ld_matrix = np.corrcoef(genotypes.T)

    # it can happen that there are no SNPs available after this filter , we return None.
    if genotypes.shape[1] == 0:
        return None, None

    try:
        ld_matrix, ordered_snps = remove_highly_correlated(ld_matrix, genotype_matrix_snps, max_correlation)
    except:
        return None, None

    return ld_matrix, ordered_snps

def read_ld_matrix_plink(bed_prepend, list_of_snps, maf_threshold, max_correlation, tmp_prepend):

    stderr = subprocess.DEVNULL if not verbosity else None
    stdout = subprocess.DEVNULL if not verbosity else None

    ## Make a correlation matrix of these SNPs
    files_to_remove = []
    snp_file = f'{tmp_prepend}_overlapping_snps.txt'
    plink_out = f'{tmp_prepend}_plink_correlation'

    with open(snp_file, 'w') as f:
        for snp in sorted(list_of_snps):
            f.write(f'{snp}\n')

    subprocess.run(['plink',
                    '--bfile', bed_prepend,
                    '--extract', snp_file,
                    '--maf', f'{maf_threshold}',
                    '--make-just-bim',
                    '--r', 'square', 'bin',
                    '--out', plink_out], check=True, stderr=stderr, stdout=stdout)

    m_snps_from_log, n_individuals = match_n_variants_n_individuals_from_plink_log(f'{plink_out}.log')
    if verbosity:
        print(m_snps_from_log, n_individuals)

    files_to_remove += [f'{plink_out}.{x}' for x in ['bim', 'log', 'ld.bin', 'nosex']]
    files_to_remove.append(snp_file)

    # read in the correlation matrix
    correlation_mat = np.fromfile(f'{plink_out}.ld.bin')
    correlation_mat = correlation_mat.reshape((int(np.floor(np.sqrt(correlation_mat.shape))),
                                               int(np.floor(np.sqrt(correlation_mat.shape)))))

    ordered_snps = []
    with open(f'{plink_out}.bim', 'r') as f:
        for i, line in enumerate(f):
            ordered_snps.append(line.split()[1])

    correlation_mat, ordered_snps = remove_highly_correlated(correlation_mat, ordered_snps, max_correlation)

    ## cleanup
    for filename in files_to_remove:
        os.remove(filename)

    return correlation_mat, ordered_snps


def remove_highly_correlated(ld_matrix, snp_ordering, max_correlation):

    # remove NA's from correlation matrix because turns out there are still some (didnt expect)
    idxs_to_remove = np.any(np.isnan(ld_matrix), axis=1)

    # remove highly correlated, as it could be a source of eigenvalue non-convergence
    to_remove = set()
    indexes = np.where(np.tril(np.abs(ld_matrix), k=-1) >= max_correlation)
    for a, b in zip(indexes[0], indexes[1]):  # I don't understand why this happens in two tuples.
        if a in to_remove:
            continue
        elif b in to_remove:
            continue
        else:
            to_remove.add(b)  # keep the first row
            idxs_to_remove[b] = True

    ld_matrix = ld_matrix[~idxs_to_remove, :][:, ~idxs_to_remove]

    #Ugly, but works.
    ordered_snps = pd.DataFrame({"snp_names": snp_ordering})
    ordered_snps = ordered_snps.loc[~idxs_to_remove, 'snp_names'].tolist()

    return ld_matrix, ordered_snps



def read_multiple_regions_from_parquet_duckdb(filename: str, regions, position_name='base_pair_location'):
    # Create the dataset object
    dataset = ds.dataset(filename, format='parquet', partitioning='hive')

    # Connect to DuckDB
    con = duckdb.connect()

    # Register the dataset in DuckDB
    con.register('sumstats', dataset)

    # Construct the SQL query to fetch data for all regions
    query_conditions = []
    for region in regions:
        chromosome = region.chromosome
        start_position = region.start
        end_position = region.end
        condition = f"(chromosome = {chromosome} AND {position_name} > {start_position} AND {position_name} <= {end_position})"
        query_conditions.append(condition)

    query = "SELECT * FROM sumstats WHERE " + " OR ".join(query_conditions) + ";"

    # Execute the query and fetch the result as a DataFrame
    df = con.execute(query).df()
    return df

def read_regions_from_file(filename, region_or_regions):

    original_input_to_gwas_catalog_dict = {'pos_name': 'rsid', 'position': 'base_pair_location', 'reference_allele': 'other_allele',
                                           'se':'standard_error', 'pval': 'p_value', 'n_iids': 'n'}

    if isinstance(region_or_regions, StartEndRegion):
        regions = [region_or_regions]
    elif isinstance(region_or_regions, list):
        regions = region_or_regions
    else:
        raise NotImplementedError(f'Only StartEndRegion or list allowed in the constructor found: {region_or_regions.__type__}')

    if filename[-5:] == '.parq':
        ## parquet file

        try:
            return_df = read_multiple_regions_from_parquet_duckdb(filename, regions)
        except:
             return_df = read_multiple_regions_from_parquet_duckdb(filename, regions, position_name='position')

        return_df = return_df.rename(columns=original_input_to_gwas_catalog_dict)
        return return_df

    else:
        regions_in_df = []
        full_df = pd.read_csv(filename, sep='\t')
        full_df = full_df.rename(columns=original_input_to_gwas_catalog_dict)
        for region in regions:
            regions_in_df.append(full_df[(full_df['chromosome'].astype(str) == region.chromosome) &
                                         (full_df['base_pair_location'].astype(int) > region.start) &
                                         (full_df['base_pair_location'].astype(int) <= region.end)
                                         ])
        return pd.concat(regions_in_df)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="""                                    
                MR-link 2: 

Pleiotropy robust cis Mendelian randomization

""",
    )
    parser.add_argument('--reference_bed',
                        required=True,
                        help='The plink bed file prepend of the genotype file that can be used as an LD reference. '
                             'Usage is the same as in the plink --bed command')
    parser.add_argument('--sumstats_exposure',
                        required=True,
                        help='The summary statistics file of the exposure file. Please see the README file or the\n'
                             'example_files folder for examples on how to make these files.')
    parser.add_argument('--sumstats_outcome',
                        required=True,
                        help='The summary statistics file of the outcome file. Please see the README file or the\n'
                             'example_files folder for examples on how to make these files.'
                             'You can specify more than one outcome per exposure',
                        action='extend',
                        nargs='+',
                        type=str
                        )
    parser.add_argument('--out',
                        required=True,
                        help='The path where to output results'
                        )

    parser.add_argument('--tmp',
                        default='tmp_',
                        help='a prepend on where to save temporary files DEPRECATED')

    parser.add_argument('--p_threshold',
                        default=5e-8,
                        help='The P value threshold for which select regions. This is the same as the clumping P value'
                             ' threshold'
                        )
    parser.add_argument('--region_padding',
                        default=5e5,
                        help='The base pair padding (on one side) on each clumped SNP which defines the region in which '
                             'MR-link 2 will perform its inference'
                        )
    parser.add_argument('--maf_threshold',
                        default=0.005,
                        help='The minor allele frequency threshold used for clumping, and for calculating the LD matrix.'
                             'This will not be applied to the summary statistics files'
                        )
    parser.add_argument('--max_correlation',
                        default=0.99,
                        help='The maximum correlation allowed in the LD matrix, if the correlation is higher than this '
                             'value between a pair of SNPs, will only keep one of them. This value is used to reduce '
                             'eigenvalue decomposition failures'
                        )

    parser.add_argument('--max_missingness',
                        default=0.95,
                        help='This is the maximum amount of individual missingness that is allowed in a summary statistic '
                             'MR-link 2 can be sensitive to large differences in summary statistic missingness, so '
                             'by default each SNP association should have at least 0.95 of observations available.')

    parser.add_argument('--var_explained_grid',
                        nargs='+',
                        default=[0.99],
                        type=float,
                        help='This field specifies the amount of variance explained of the LD matrix that is used by MR'
                             '-link 2. You can add onto this field, and all variances explained will be added: '
                             '--var_explained_grid 0.99 0.999 0.2 0.96 0.1 will perform an MR-link 2 estimate for '
                             'all these values.')

    parser.add_argument('--continue_analysis',
                        action='store_true',
                        help='Flag to continue an already started analysis, if specified this will look for a temporary '
                             'file, and if it present, reuse its results. This can be handy if you have hundreds of '
                             'associated regions, which can sometimes take a long time to run.')

    parser.add_argument('--no_normalize_sumstats',
                        action='store_true',
                        help='flag to _not_ normalize summary statistics')

    parser.add_argument('--regions_to_read_at_the_same_time',
                        type=int,
                        default=10,
                        help='Number of regions to read at the same time from a file. Will reduce I/O times, but can increase memory usage')

    parser.add_argument('--no_exclude_hla',
                        action='store_true',
                        help='flag to _not_ exclude the HLA region.')

    parser.add_argument('--prespecified_regions',
                        required=False,
                        help='Specify which regions to do. format is the following: '
                             '`{chr_1}:{start_1}-{end_1},{chr_2}:{start_2}-{end_2}`. i.e. these are regions that are '
                             'separated with the comma character.'
                             'This step will skip the clumping step to identify regions, '
                             'and will blindly perform MR on regions that may or may not be associated to the exposure.\n'
                             'This feature has been included to allow investigators to identify associated regions in '
                             'one discovery cohort, and assess the MR effect in a different cohort, this is done '
                             'to remove winners curse.'
                        )

    parser.add_argument('--max_snps_in_correlation_matrix',
                        required=False,
                        default=5250,
                        help='How many SNPs are allowed in the correlation matrix, as it can take a lot of time to '
                             'Do the eigendecompositions and Rscripts otherwise.'
                        )

    parser.add_argument('--verbose',
                        default=0,
                        help='Set to 1 if you want to read more output, for debugging purposes ')

    args = parser.parse_args()

    # necessary params
    reference_bed = args.reference_bed
    sumstats_exposure = args.sumstats_exposure
    sumstats_outcome = args.sumstats_outcome

    # optional params that could influence results
    p_threshold = float(args.p_threshold)
    region_padding = int(float(args.region_padding))
    maf_threshold = float(args.maf_threshold)
    max_correlation = float(args.max_correlation)

    max_correlation_matrix_snps = int(args.max_snps_in_correlation_matrix)

    read_chunk = int(args.regions_to_read_at_the_same_time)


    max_missingness = float(args.max_missingness)
    if 0.0 > max_missingness > 1.0:
        raise ValueError('Missingness needs to be between 0.0 and 1.0')

    var_explained_grid = args.var_explained_grid
    for var_explained_val in var_explained_grid:
        if 0.0 > var_explained_val >= 1.0:
            raise ValueError('All var explained values should be in the [0.0, 1.0] interval')

    # housekeeping params
    verbosity = args.verbose
    files_to_remove = set()

    plink_geno_obj = PlinkGenoReader(reference_bed)

    original_input_to_gwas_catalog_dict = {'pos_name': 'rsid', 'position': 'base_pair_location', 'reference_allele': 'other_allele',
                                           'se':'standard_error', 'pval': 'p_value', 'n_iids': 'n'}

    with tempfile.TemporaryDirectory() as tmp_prepend:
        tmp_dir = f'{tmp_prepend}/tmp_'
        # This is some preprocessing of the summary statistic files to make sure that we filter for at least 95% sample size
        sumstats_necessary_colnames = {'rsid', 'chromosome', 'base_pair_location','effect_allele', 'other_allele',
                                       'beta', 'standard_error', 'p_value', 'n'}

        print(f'Performing sumstats file preprocessing')
        """
        Exposure summary statistics loading and parsing
        """
        tmp_exposure_sumstats_file = f'{tmp_dir}_exposure_sumstats.txt'
        files_to_remove.add(tmp_exposure_sumstats_file)

        ## Can also be a parquet file.
        if sumstats_exposure[-5:] == '.parq':
            exposure_df = pd.read_parquet(sumstats_exposure)
        else:
            exposure_df = pd.read_csv(sumstats_exposure, sep='\t')

        exposure_df = exposure_df.rename(columns=original_input_to_gwas_catalog_dict)
        if len(set(exposure_df.columns) & sumstats_necessary_colnames) != len(sumstats_necessary_colnames):
            raise ValueError('Exposure summary statistics do not contain the necessary columns.\n'
                             f'The following columns should at least be present:\n{sumstats_necessary_colnames}\n'
                             f'The following columns are present: {exposure_df.columns}')

        if verbosity:
            print(f'before missingness filter: {exposure_df.shape=}')


        if not args.no_normalize_sumstats:
            exposure_df['z'] = exposure_df['beta'] / exposure_df['standard_error']
            exposure_df['beta'] = exposure_df.z / np.sqrt(exposure_df.n + exposure_df.z ** 2)
            exposure_df['standard_error'] = 1 / np.sqrt(exposure_df.n + exposure_df.z ** 2)

        exposure_df.to_csv(tmp_exposure_sumstats_file, sep='\t', index=False)
        if verbosity:
            print(f'after missingness filter: {exposure_df.shape=}')

        """
        Perform clumping or use regions that were previously prespecified.
        """

        if args.prespecified_regions is not None:
            exposure_regions = [StartEndRegion(x) for x in args.prespecified_regions.split(',')]
            print(f"The following regions have been prespecified: {exposure_regions}. Not doing clumping")
        else:
            print('Starting to identify regions through clumping')
            exposure_regions = identify_regions(tmp_exposure_sumstats_file,
                                                reference_bed,
                                                p_threshold,
                                                maf_threshold,
                                                region_padding,
                                                tmp_dir,
                                                verbosity_level=verbosity)

            if exposure_regions is None:
                print(f'Did not find any regions for {sumstats_exposure}, at p < {p_threshold:.2e}')
                # write to an empty file
                with open(args.out + "_noregions", 'w') as f:
                    f.write('\n')
                exit()

        regions_to_do = StartEndRegions(exposure_regions)

        if len(regions_to_do.gene_regions) == 0:
            print(f'Did not find any regions for {sumstats_exposure}, at p < {p_threshold:.2e}')
            # write to an empty file
            with open(args.out + "_no_overlap", 'w') as f:
                f.write('\n')
            exit()

        if not args.no_exclude_hla:
            hla_region = StartEndRegion('6:25000000-37000000')
            regions_to_do = [StartEndRegion(x) for x in regions_to_do if StartEndRegion(x) not in hla_region]
        else:
            regions_to_do = [StartEndRegion(x) for x in regions_to_do]

        print(f'Finished identifying {len(regions_to_do)} regions, now continueing with MR-link2 for each region')

        if verbosity:
            print(exposure_df.head())



        all_results = []
        outcome_region_combinations_already_done = []
        previously_done_df = pd.DataFrame()

        combined_df = None

        if args.continue_analysis:
            if os.path.exists(args.out + '_tmp'):
                previously_done_df = pd.read_csv(args.out + '_tmp', sep='\t')
                outcome_region_combinations_already_done = set(zip(previously_done_df.region, previously_done_df.outcome_file))

            elif os.path.exists(args.out):
                previously_done_df = pd.read_csv(args.out, sep='\t')
                outcome_region_combinations_already_done = set(zip(previously_done_df.region, previously_done_df.outcome_file))

            all_results.append(previously_done_df)
            combined_df = pd.concat(all_results)

        exceptions = []
        if os.path.exists(args.out + '_no_estimate'):
            with open(args.out + '_no_estimate', 'a') as f:
                for region, exception in exceptions:
                    f.write(f'{region}\t{exception}\n')


        """
        HERE, we perform the reading of chunks 
        """
        for i in range(0, len(regions_to_do), read_chunk):
            regions_chunk = regions_to_do[i:i + read_chunk]
            chunked_regions_per_outcome = {}
            """
            This is an optimization and I don't love it. 
            But can save like 5x of the time, so I'll take the time saving
            
            In essence what we do is load in up to <read_chunks> regions of the outcome files.
            After which we apply MR-link-2 on all the loaded files. 
            
            In profiling, I found that this is substantially faster than reading a region one at a time.
            
            This is especially valuable when we analyze parquet files. 
            
            """
            for outcome_location in sumstats_outcome:
                """
                It is necessary to also check if the region has been done before.
                makes it look ugly, but retains robustness and removes double work.
                """
                regions_for_this_outcome = []
                for region in regions_chunk:
                    if (str(region), outcome_location) in outcome_region_combinations_already_done:
                        print(f'Already done {outcome_location} with {region}, continueing with the next one')
                        continue
                    regions_for_this_outcome.append(region)

                if len(regions_for_this_outcome):
                    start = time.time()
                    chunked_regions_per_outcome[outcome_location] = read_regions_from_file(
                        outcome_location, regions_chunk)
                    end = time.time()

                    if True or verbosity:
                        print(f'Read {outcome_location} {len(regions_chunk)} in {end - start:.2f}s')

            """
            And now to continue with the MR-link-2
            """
            for region in regions_chunk:
                regional_exposure_df = exposure_df[
                    (exposure_df.chromosome.astype(str) == region.chromosome) &
                    (exposure_df.base_pair_location.astype(int) >= region.start) &
                    (exposure_df.base_pair_location.astype(int) <= region.end)
                    ].copy()


                regional_exposure_df = regional_exposure_df[regional_exposure_df.n >= (max_missingness * regional_exposure_df.n.max())]
                if verbosity:
                    print(f'after missingness filter for {region}: {regional_exposure_df.shape=}')

                regional_ld_matrix, snps_in_ld_matrix = read_ld_matrix_local(plink_geno_obj,
                                                                             regional_exposure_df.rsid,
                                                                             maf_threshold=maf_threshold,
                                                                             max_correlation=max_correlation,
                                                                             )
                if regional_ld_matrix is None:
                    exceptions.append((region, 'NO_SNPS_IN_LD_MATRIX'))
                    print(f'Unable to identify mr-link2 results in {region} region due to no SNPs in the LD matrix')
                    continue

                snps_and_alleles_in_ld_matrix = [(x, plink_geno_obj.bim_data[x][2], plink_geno_obj.bim_data[x][3]) for x
                                                 in
                                                 snps_in_ld_matrix]

                for outcome_location in sumstats_outcome:

                    regional_results = None
                    if (str(region), outcome_location) in outcome_region_combinations_already_done:
                        print(f'Already done {outcome_location} on {region}, continueing with the next one')
                        continue

                    """
                    Outcome summary statistics loading and parsing
                    """
                    outcome_df = chunked_regions_per_outcome[outcome_location]
                    regional_outcome_df = outcome_df[
                        (outcome_df.chromosome.astype(str) == region.chromosome) &
                        (outcome_df.base_pair_location.astype(int) >= region.start) &
                        (outcome_df.base_pair_location.astype(int) <= region.end)
                        ].copy()

                    regional_outcome_df = regional_outcome_df.rename(columns=original_input_to_gwas_catalog_dict)

                    if len(set(outcome_df.columns) & sumstats_necessary_colnames) != len(sumstats_necessary_colnames):
                        raise ValueError('Outcome summary statistics do not contain the necessary columns.\n'
                                         f'The following columns should at least be present:\n{sumstats_necessary_colnames}\n'
                                         f'The following columns are present: {outcome_df.columns}')
                    if verbosity:
                        print(f'before missingness filter: {regional_outcome_df.shape=}')
                    regional_outcome_df = regional_outcome_df[regional_outcome_df.n >= (max_missingness * regional_outcome_df.n.max())]

                    if not args.no_normalize_sumstats:
                        regional_outcome_df['z'] = regional_outcome_df['beta'] / regional_outcome_df['standard_error']
                        regional_outcome_df['beta'] = regional_outcome_df.z / np.sqrt(regional_outcome_df.n + regional_outcome_df.z ** 2)
                        regional_outcome_df['standard_error'] = 1 / np.sqrt(regional_outcome_df.n + regional_outcome_df.z ** 2)

                    if verbosity:
                        print(f'after missingness filter: {regional_outcome_df.shape=}')


                    try:
                        regional_results = mr_link2_on_region(region,
                                                              regional_exposure_df,
                                                              regional_outcome_df,
                                                              reference_bed,
                                                              maf_threshold,
                                                              max_correlation,
                                                              regional_ld_matrix,
                                                              snps_alleles_in_ld_matrix=snps_and_alleles_in_ld_matrix,
                                                              tmp_prepend=tmp_dir,
                                                              verbosity=verbosity,
                                                              var_explained_grid=var_explained_grid,
                                                              max_snp_threshold=max_correlation_matrix_snps)
                    except Exception as x:
                        exceptions.append((region, x))
                        print(f'Unable to make an MR-link2 estimate in {region} due to {x}')
                        continue

                    if regional_results is None:
                        exceptions.append((region, 'NONE_RESULT'))
                        print(f'Unable to identify mr-link2 results in {region} region')
                        continue

                    regional_results['exposure_file'] = sumstats_exposure
                    regional_results['outcome_file'] = outcome_location

                    all_results.append(regional_results)
                    combined_df = pd.concat(all_results)
                    combined_df.to_csv(args.out + '_tmp', sep='\t', index=False)

                # write exceptions to a file
                if len(exceptions) != 0:
                    with open(args.out + '_no_estimate', 'a') as f:
                        for region, exception in exceptions:
                            f.write(f'{region}\t{exception}\n')

        # finally, clean up
        for filename in files_to_remove:
            if os.path.exists(filename):
                os.remove(filename)

        # Finally, write results to txt.
        if len(all_results) != 0 and combined_df is not None:
            combined_df.to_csv(args.out, sep='\t', index=False)
            if os.path.exists(args.out + '_tmp'):
                os.remove(args.out + '_tmp')


