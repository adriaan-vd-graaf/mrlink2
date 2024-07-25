import sys
import os
import hashlib
import subprocess
import tempfile
import numpy as np
import pandas as pd
import scipy
import pytest


sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from mr_link_2_standalone import *


""""
NB. These tests require us to have plink in your path. 
So unfortunately I will not test them in the github workflow.
"""

def test_plink_version():
    result = subprocess.run(['plink', '--version'], capture_output=True, text=True)
    assert result.returncode == 0, "PLINK is not installed or not in the PATH"
    assert 'PLINK v1.90' in result.stdout, "PLINK version is not 1.9"

"""
I guess these are to ensure that no regressions happen 
"""


def test_identification_of_regions_normal():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    exposure_loc = os.path.join(current_dir, '../example_files/yes_causal_exposure.txt')
    reference_bed = os.path.join(current_dir, '../example_files/reference_cohort')

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_dir = os.path.join(tmpdir, 'tmp_integration_testing_')

        try:
            regions = identify_regions(exposure_loc,
                                       reference_bed,
                                       5e-8,
                                       0.01,
                                       padding=250_000,
                                       plink_tmp_prepend=tmp_dir,
                                       r_sq_threshold=0.01,
                                       verbosity_level=1)  # Set verbosity_level to capture output
        except subprocess.CalledProcessError as e:
            print(f"Command '{e.cmd}' returned non-zero exit status {e.returncode}.")
            print(f"Output: {e.output}")
            print(f"Error: {e.stderr}")
            raise

    assert len(regions) == 1
    assert regions[0].chromosome == '2'
    assert regions[0].end == 103230976
    assert regions[0].start == 101773297


def test_identification_of_regions_small_padding():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    exposure_loc = os.path.join(current_dir, '../example_files/non_causal_exposure.txt')
    reference_bed = os.path.join(current_dir, '../example_files/reference_cohort')

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_dir = os.path.join(tmpdir, 'tmp_integration_testing_')
        try:
            regions = identify_regions(
                exposure_loc,
                reference_bed,
                5e-12,
                0.01,
                padding=5_000,
                plink_tmp_prepend=tmp_dir,
                r_sq_threshold=0.01,
                verbosity_level=1  # Increase verbosity
            )
        except subprocess.CalledProcessError as e:
            print(f"Command '{e.cmd}' returned non-zero exit status {e.returncode}.")
            print(f"Output: {e.output}")
            print(f"Error: {e.stderr}")
            raise
        except subprocess.TimeoutExpired as e:
            print(f"Command '{e.cmd}' timed out after {e.timeout} seconds.")
            raise

        reference_regions = [
            '2:101995210-102141235', '2:102143022-102307978',
            '2:102310933-102518306', '2:102519095-102529095',
            '2:102529412-102593218', '2:102595083-102659228',
            '2:102659423-103001345'
        ]
        assert len(regions) == 7
        assert [str(x) for x in regions] == reference_regions

""""
Integration / regression tests for MR-link-2
"""

"""
Going for byte for byte comparisons. Will probably not work on windows.
"""
def get_md5sum(file_path):
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def test_mr_link_2_integration_yes_causal():

    try:
        import numpy as np
    except ImportError:
        subprocess.run([sys.executable, "-m", "pip", "install", "numpy", "pandas", 'scipy'], check=True)
        import numpy as np

    current_dir = os.path.dirname(os.path.abspath(__file__))

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_loc = f'{tmpdir}/tmp_integration_testing_'

        command = [
            sys.executable, f"{current_dir}/../mr_link_2_standalone.py",
            "--reference_bed", f"{current_dir}/../example_files/reference_cohort",
            "--sumstats_exposure", f"{current_dir}/../example_files/yes_causal_exposure.txt",
            "--sumstats_outcome", f"{current_dir}/../example_files/yes_causal_outcome.txt",
            "--out", f'{tmp_loc}yes_causal.txt'
        ]

        result = subprocess.run(command, capture_output=True, text=True)
        data_frame = pd.read_csv(f'{tmp_loc}yes_causal.txt', sep='\t')

        assert result.returncode == 0, f"Command failed with return code {result.returncode}. Output: {result.stdout} Error: {result.stderr}"

        assert np.isclose(data_frame.alpha, 0.4193872544798938)
        assert np.isclose(data_frame['se(alpha)'], 0.0565878953774011)
        assert np.isclose(data_frame['p(alpha)'], 1.25110890260147e-13)

        assert np.isclose(data_frame.sigma_x,0.564994403645872)
        assert np.isclose(data_frame.sigma_y, 0.1549295595449322)
        assert np.isclose(data_frame['se(sigma_y)'], 0.0062419901963417)
        assert np.isclose(data_frame['p(sigma_y)'], 5.3812527742659674e-136)



def test_mr_link_2_integration_non_causal():

    current_dir = os.path.dirname(os.path.abspath(__file__))
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_loc = f'{tmpdir}/tmp_integration_testing_'

        command = [
            sys.executable, f"{current_dir}/../mr_link_2_standalone.py",
            "--reference_bed", f"{current_dir}/../example_files/reference_cohort",
            "--sumstats_exposure", f"{current_dir}/../example_files/non_causal_exposure.txt",
            "--sumstats_outcome", f"{current_dir}/../example_files/non_causal_outcome.txt",
            "--out", f'{tmp_loc}non_causal.txt'
        ]

        result = subprocess.run(command, capture_output=True, text=True)
        data_frame = pd.read_csv(f'{tmp_loc}non_causal.txt', sep='\t')

        assert result.returncode == 0, f"Command failed with return code {result.returncode}. Output: {result.stdout} Error: {result.stderr}"

        assert np.isclose(data_frame.alpha, -0.0330966574547156)
        assert np.isclose(data_frame['se(alpha)'], 0.0529061811094349)
        assert np.isclose(data_frame['p(alpha)'], 0.5315953131580522)

        assert np.isclose(data_frame.sigma_x,0.5641943699065054)
        assert np.isclose(data_frame.sigma_y, 0.1642151919013568)
        assert np.isclose(data_frame['se(sigma_y)'], 0.0066189399600698)
        assert np.isclose(data_frame['p(sigma_y)'], 7.011323342257353e-136)



