import pytest
import os
import sys
import tempfile

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


from mr_link_2_standalone import *


"""
StartEndRegion class tests
"""
class MockSNP:
    def __init__(self, chromosome, position):
        self.chromosome = chromosome
        self.position = position

def test_startendregion_init_with_list():
    region = StartEndRegion(['1', 100, 200])
    assert region.chromosome == '1'
    assert region.start == 100
    assert region.end == 200

def test_startendregion_init_with_three_arguments():
    region = StartEndRegion('1', 100, 200)
    assert region.chromosome == '1'
    assert region.start == 100
    assert region.end == 200

def test_startendregion_init_with_string():
    region = StartEndRegion('1:100-200')
    assert region.chromosome == '1'
    assert region.start == 100
    assert region.end == 200

def test_startendregion_init_with_another_region():
    region1 = StartEndRegion('1', 100, 200)
    region2 = StartEndRegion(region1)
    assert region2.chromosome == '1'
    assert region2.start == 100
    assert region2.end == 200

def test_startendregion_init_with_invalid_arguments():
    with pytest.raises(ValueError):
        StartEndRegion('invalid')

def test_startendregion_invalid_start_end_positions():
    with pytest.raises(RuntimeError):
        StartEndRegion('1', 200, 100)

def test_startendregion_negative_positions():
    with pytest.raises(RuntimeError):
        StartEndRegion('1', -100, 200)

def test_startendregion_position_in_region():
    region = StartEndRegion('1', 100, 200)
    assert region.position_in_region('1', 150) is True
    assert region.position_in_region('1', 250) is False

def test_startendregion_snp_in_region():
    region = StartEndRegion('1', 100, 200)
    assert region.snp_in_region('1', 150) is True
    assert region.snp_in_region('1', 250) is False

def test_startendregion_snp_object_in_region():
    region = StartEndRegion('1', 100, 200)
    snp = MockSNP('1', 150)
    assert region.snp_object_in_region(snp) is True
    snp_outside = MockSNP('1', 250)
    assert region.snp_object_in_region(snp_outside) is False

def test_startendregion_region_overlaps():
    region1 = StartEndRegion('1', 100, 200)
    region2 = StartEndRegion('1', 150, 250)
    region3 = StartEndRegion('1', 300, 400)
    assert region1.region_overlaps(region2) is True
    assert region1.region_overlaps(region3) is False

def test_startendregion_str():
    region = StartEndRegion('1', 100, 200)
    assert str(region) == '1:100-200'

def test_startendregion_lt():
    region1 = StartEndRegion('1', 100, 200)
    region2 = StartEndRegion('1', 150, 250)
    assert (region1 < region2) is True
    assert (region2 < region1) is False

def test_startendregion_gt():
    region1 = StartEndRegion('1', 100, 200)
    region2 = StartEndRegion('1', 150, 250)
    region3 = StartEndRegion('1', 201, 250)
    assert (region1 > region2) is False
    assert (region2 > region1) is False
    assert (region3 > region1) is True

def test_startendregion_contains():
    region1 = StartEndRegion('1', 100, 200)
    region2 = StartEndRegion('1', 150, 250)
    region3 = StartEndRegion('1', 300, 400)
    assert (region2 in region1) is True
    assert (region3 in region1) is False

def test_startendregion_repr():
    region = StartEndRegion('1', 100, 200)
    assert repr(region) == '1:100-200'


"""
StartEndRegions class tests
"""

def test_startendregions_init():
    regions = StartEndRegions([['1', 100, 200], ['1', 300, 400]])
    assert len(regions.gene_regions) == 2
    assert regions.gene_regions[0].chromosome == '1'
    assert regions.gene_regions[0].start == 100
    assert regions.gene_regions[0].end == 200

def test_startendregions_in_gene_regions():
    regions = StartEndRegions([['1', 100, 200], ['1', 300, 400]])
    assert regions.in_gene_regions('1', 150) is True
    assert regions.in_gene_regions('1', 350) is True
    assert regions.in_gene_regions('1', 250) is False

def test_startendregions_make_non_overlapping_regions():
    regions = StartEndRegions([['1', 100, 200], ['1', 150, 250], ['1', 300, 400]])
    non_overlapping = regions.make_non_overlapping_regions()
    assert len(non_overlapping.gene_regions) == 2
    assert non_overlapping.gene_regions[0].start == 100
    assert non_overlapping.gene_regions[0].end == 250
    assert non_overlapping.gene_regions[1].start == 300
    assert non_overlapping.gene_regions[1].end == 400

def test_startendregions_iteration():
    regions = StartEndRegions([['1', 100, 200], ['1', 300, 400]])
    iter_regions = iter(regions)
    assert next(iter_regions) == regions.gene_regions[0]
    assert next(iter_regions) == regions.gene_regions[1]
    with pytest.raises(StopIteration):
        next(iter_regions)

def test_startendregions_repr():
    regions = StartEndRegions([['1', 100, 200], ['1', 300, 400]])
    assert repr(regions) == "StartEndRegions with 2 regions across chromosome(s) ['1']"

def test_startendregions_contains():
    regions = StartEndRegions([['1', 100, 200], ['1', 300, 400]])
    region_to_check = StartEndRegion(['1', 150, 180])
    assert region_to_check in regions
    region_to_check = StartEndRegion(['1', 250, 280])
    assert region_to_check not in regions


"""
Now unit tests for the functions.
"""

def test_identify_regions():
    ## TODO, difficult to get right.
    pass

"""
Now, we perform unit tests on the MR-link-2 functions. to ensure that there is no creep of the 
First function v0
"""
def test_mr_link2_loglik_reference_v0_basic():
    th = np.array([0.1, 0.01, 0.01])
    lam = np.array([1, 2, 3])
    c_x = np.array([1, 1, 1])
    c_y = np.array([2, 2, 2])
    n_x = 100
    n_y = 100
    result = mr_link2_loglik_reference_v0(th, lam, c_x, c_y, n_x, n_y)
    assert np.isclose(result, 22.944216071347796)
    assert np.isclose(result, mr_link2_loglik_reference_v2(th, lam, c_x, c_y, n_x, n_y))


def test_mr_link2_loglik_reference_v0_zero_lam():
    th = np.array([0.1, 0.01, 0.01])
    lam = np.array([0, 0, 0])
    c_x = np.array([1, 1, 1])
    c_y = np.array([2, 2, 2])
    n_x = 100
    n_y = 100
    assert np.isnan(mr_link2_loglik_reference_v0(th, lam, c_x, c_y, n_x, n_y))


def test_mr_link2_loglik_reference_v0_large_values():
    th = np.array([0.5, 1e10, 1e10])
    lam = np.array([1e5, 1e5, 1e5])
    c_x = np.array([1e2, 1e2, 1e2])
    c_y = np.array([1e2, 1e2, 1e2])
    n_x = 1e6
    n_y = 1e6
    result = mr_link2_loglik_reference_v0(th, lam, c_x, c_y, n_x, n_y)
    assert np.isclose(result, 17617.166270043643,)
    assert np.isclose(result, mr_link2_loglik_reference_v2(th, lam, c_x, c_y, n_x, n_y))

def test_mr_link2_loglik_reference_v0_negative_values():
    th = np.array([0.1, -0.01, 0.01])
    lam = np.array([1, 2, 3])
    c_x = np.array([1, 1, 1])
    c_y = np.array([2, 2, 2])
    n_x = 100
    n_y = 100
    result = mr_link2_loglik_reference_v0(th, lam, c_x, c_y, n_x, n_y)
    assert np.isclose(result, 22.944216071347796)
    assert np.isclose(result, mr_link2_loglik_reference_v2(th, lam, c_x, c_y, n_x, n_y))


def test_mr_link2_loglik_reference_v2_basic():
    th = np.array([0.1, 0.01, 0.01])
    lam = np.array([1, 2, 3])
    c_x = np.array([1, 1, 1])
    c_y = np.array([2, 2, 2])
    n_x = 100
    n_y = 100
    result = mr_link2_loglik_reference_v2(th, lam, c_x, c_y, n_x, n_y)
    assert np.isclose(result, 22.944216071347796)


def test_mr_link2_loglik_reference_v2_zero_lam():
    th = np.array([0.1, 0.01, 0.01])
    lam = np.array([0, 0, 0])
    c_x = np.array([1, 1, 1])
    c_y = np.array([2, 2, 2])
    n_x = 100
    n_y = 100
    assert np.isnan(mr_link2_loglik_reference_v2(th, lam, c_x, c_y, n_x, n_y))


def test_mr_link2_loglik_reference_v2_large_values():
    th = np.array([0.5, 1e10, 1e10])
    lam = np.array([1e5, 1e5, 1e5])
    c_x = np.array([1e2, 1e2, 1e2])
    c_y = np.array([1e2, 1e2, 1e2])
    n_x = 1e6
    n_y = 1e6
    result = mr_link2_loglik_reference_v2(th, lam, c_x, c_y, n_x, n_y)
    assert isinstance(result, float)
    assert np.isclose(result, 17617.166270043643,)


def test_mr_link2_loglik_reference_v2_negative_values():
    th = np.array([0.1, -0.01, 0.01])
    lam = np.array([1, 2, 3])
    c_x = np.array([1, 1, 1])
    c_y = np.array([2, 2, 2])
    n_x = 100
    n_y = 100
    result = mr_link2_loglik_reference_v2(th, lam, c_x, c_y, n_x, n_y)
    assert isinstance(result, float)
    assert np.isclose(result, 22.944216071347796)
    assert np.isclose(result, mr_link2_loglik_reference_v0(th, lam, c_x, c_y, n_x, n_y))

"""
As we do use them as functions, we want to make sure that the h0 functions also work correctly.
"""
def test_mr_link2_loglik_alpha_h0_basic():
    th = np.array([0.01, 0.01])
    lam = np.array([1, 2, 3])
    cX = np.array([1, 1, 1])
    cY = np.array([2, 2, 2])
    nX = 100
    nY = 100
    result = mr_link2_loglik_alpha_h0(th, lam, cX, cY, nX, nY)
    assert result == mr_link2_loglik_reference_v2([0.0] + list(th), lam, cX, cY, nX, nY)


def test_mr_link2_loglik_alpha_h0_zero_lam():
    th = np.array([0.01, 0.01])
    lam = np.array([0, 0, 0])
    cX = np.array([1, 1, 1])
    cY = np.array([2, 2, 2])
    nX = 100
    nY = 100
    result = mr_link2_loglik_alpha_h0(th, lam, cX, cY, nX, nY)
    assert np.isnan(result)
    assert np.isnan(mr_link2_loglik_reference_v2([0.0] + list(th), lam, cX, cY, nX, nY))

def test_mr_link2_loglik_alpha_h0_large_values():
    th = np.array([1e10, 1e10])
    lam = np.array([1e5, 1e5, 1e5])
    cX = np.array([1e2, 1e2, 1e2])
    cY = np.array([1e2, 1e2, 1e2])
    nX = 1e6
    nY = 1e6
    result = mr_link2_loglik_alpha_h0(th, lam, cX, cY, nX, nY)
    assert isinstance(result, float)
    assert result == mr_link2_loglik_reference_v2([0.0] + list(th), lam, cX, cY, nX, nY)

def test_mr_link2_loglik_alpha_h0_negative_values():
    th = np.array([0.01, -0.01])
    lam = np.array([1, 2, 3])
    cX = np.array([1, 1, 1])
    cY = np.array([2, 2, 2])
    nX = 100
    nY = 100
    result = mr_link2_loglik_alpha_h0(th, lam, cX, cY, nX, nY)
    assert isinstance(result, float)
    assert result == mr_link2_loglik_reference_v2([0.0] + list(th), lam, cX, cY, nX, nY)

def test_mr_link2_loglik_alpha_h0_large_population():
    th = np.array([0.01, 0.01])
    lam = np.array([1, 2, 3])
    cX = np.array([1, 1, 1])
    cY = np.array([2, 2, 2])
    nX = 1e12
    nY = 1e12
    result = mr_link2_loglik_alpha_h0(th, lam, cX, cY, nX, nY)
    assert isinstance(result, float)
    assert result == mr_link2_loglik_reference_v2([0.0] + list(th), lam, cX, cY, nX, nY)

