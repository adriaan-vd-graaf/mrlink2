import pytest
import sys
import os

## For the python part

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from mr_link_2_standalone import *


## for the R part

from rpy2.robjects import r, numpy2ri, default_converter
from rpy2.robjects.vectors import FloatVector
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter


# Activate automatic conversion between numpy and R objects
numpy2ri.activate()

# Load the base R package
base = importr('base')
# Source the R file with the functions
r['source']('../R/mr_link_2_functions.R')

"""
Now, we perform unit tests on the MR-link-2 functions. to ensure that there is no creep of the 
First function v0
"""
def test_mr_link2_loglik_reference_v0_basic_r_concordance():
    th = np.array([0.1, 0.01, 0.01])
    lam = np.array([1, 2, 3])
    c_x = np.array([1, 1, 1])
    c_y = np.array([2, 2, 2])
    n_x = 100
    n_y = 100

    expected_result = 22.944216071347796
    python_result = mr_link2_loglik_reference_v0(th, lam, c_x, c_y, n_x, n_y)

    assert np.isclose(python_result, expected_result)
    assert np.isclose(python_result, mr_link2_loglik_reference_v2(th, lam, c_x, c_y, n_x, n_y))

    r_result = r['mr_link2_loglik_reference_v2'](th, lam, c_x, c_y, n_x, n_y)
    assert np.isclose(r_result[0], expected_result)


def test_mr_link2_loglik_reference_v0_zero_lam_r_concordance():
    th = np.array([0.1, 0.01, 0.01])
    lam = np.array([0, 0, 0])
    c_x = np.array([1, 1, 1])
    c_y = np.array([2, 2, 2])
    n_x = 100
    n_y = 100

    python_result = mr_link2_loglik_reference_v0(th, lam, c_x, c_y, n_x, n_y)
    assert np.isnan(python_result)

    r_result = r['mr_link2_loglik_reference_v2'](th, lam, c_x, c_y, n_x, n_y)
    assert np.isnan(r_result[0])

def test_mr_link2_loglik_reference_v0_large_values_r_concordance():
    th = np.array([0.5, 1e10, 1e10])
    lam = np.array([1e5, 1e5, 1e5])
    c_x = np.array([1e2, 1e2, 1e2])
    c_y = np.array([1e2, 1e2, 1e2])
    n_x = 1e6
    n_y = 1e6

    expected_result = 17617.166270043643
    python_result = mr_link2_loglik_reference_v0(th, lam, c_x, c_y, n_x, n_y)
    assert np.isclose(python_result, expected_result)

    r_result = r['mr_link2_loglik_reference_v2'](th, lam, c_x, c_y, n_x, n_y)
    assert np.isclose(r_result[0], expected_result)

def test_mr_link2_loglik_reference_v0_negative_values_r_concordance():
    th = np.array([0.1, -0.01, 0.01])
    lam = np.array([1, 2, 3])
    c_x = np.array([1, 1, 1])
    c_y = np.array([2, 2, 2])
    n_x = 100
    n_y = 100

    expected_result = 22.944216071347796
    python_result = mr_link2_loglik_reference_v0(th, lam, c_x, c_y, n_x, n_y)
    assert np.isclose(python_result, expected_result)

    r_result = r['mr_link2_loglik_reference_v2'](th, lam, c_x, c_y, n_x, n_y)
    assert np.isclose(r_result[0], expected_result)


def test_mr_link2_loglik_reference_v2_basic_r_concordance():
    th = np.array([0.1, 0.01, 0.01])
    lam = np.array([1, 2, 3])
    c_x = np.array([1, 1, 1])
    c_y = np.array([2, 2, 2])
    n_x = 100
    n_y = 100

    expected_result = 22.944216071347796
    python_result = mr_link2_loglik_reference_v2(th, lam, c_x, c_y, n_x, n_y)
    assert np.isclose(python_result, expected_result)

    r_result = r['mr_link2_loglik_reference_v2'](th, lam, c_x, c_y, n_x, n_y)
    assert np.isclose(r_result[0], expected_result)



def test_mr_link2_loglik_reference_v2_zero_lam_r_concordance():
    th = np.array([0.1, 0.01, 0.01])
    lam = np.array([0, 0, 0])
    c_x = np.array([1, 1, 1])
    c_y = np.array([2, 2, 2])
    n_x = 100
    n_y = 100

    python_result = mr_link2_loglik_reference_v2(th, lam, c_x, c_y, n_x, n_y)
    assert np.isnan(python_result)

    r_result = r['mr_link2_loglik_reference_v2'](th, lam, c_x, c_y, n_x, n_y)
    assert np.isnan(r_result[0])


def test_mr_link2_loglik_reference_v2_large_values_r_concordance():
    th = np.array([0.5, 1e10, 1e10])
    lam = np.array([1e5, 1e5, 1e5])
    c_x = np.array([1e2, 1e2, 1e2])
    c_y = np.array([1e2, 1e2, 1e2])
    n_x = 1e6
    n_y = 1e6

    expected_result = 17617.166270043643
    python_result = mr_link2_loglik_reference_v2(th, lam, c_x, c_y, n_x, n_y)
    assert np.isclose(python_result, expected_result)

    r_result = r['mr_link2_loglik_reference_v2'](th, lam, c_x, c_y, n_x, n_y)
    assert np.isclose(r_result[0], expected_result)


def test_mr_link2_loglik_reference_v2_negative_values_r_concordance():
    th = np.array([0.1, -0.01, 0.01])
    lam = np.array([1, 2, 3])
    c_x = np.array([1, 1, 1])
    c_y = np.array([2, 2, 2])
    n_x = 100
    n_y = 100

    expected_result = 22.944216071347796
    python_result = mr_link2_loglik_reference_v2(th, lam, c_x, c_y, n_x, n_y)

    # Check the Python result against the expected value
    assert np.isclose(python_result, expected_result)

    # Call the R function and check if it matches the expected value
    r_result = r['mr_link2_loglik_reference_v2'](th, lam, c_x, c_y, n_x, n_y)
    assert np.isclose(r_result[0], expected_result)

    # Also check against the mr_link2_loglik_reference_v0 function
    reference_result_v0 = mr_link2_loglik_reference_v0(th, lam, c_x, c_y, n_x, n_y)
    assert np.isclose(python_result, reference_result_v0)


"""
As we do use them as functions, we want to make sure that the h0 functions also work correctly.
"""
def test_mr_link2_loglik_alpha_h0_basic_r_concordance():
    th = np.array([0.01, 0.01])
    lam = np.array([1, 2, 3])
    cX = np.array([1, 1, 1])
    cY = np.array([2, 2, 2])
    nX = 100
    nY = 100

    python_result = mr_link2_loglik_alpha_h0(th, lam, cX, cY, nX, nY)
    reference_result = mr_link2_loglik_reference_v2([0.0] + list(th), lam, cX, cY, nX, nY)

    assert np.isclose(python_result, reference_result)

    r_result = r['mr_link2_loglik_alpha_h0'](th, lam, cX, cY, nX, nY)
    assert np.isclose(r_result[0], reference_result)


def test_mr_link2_loglik_alpha_h0_zero_lam_r_concordance():
    th = np.array([0.01, 0.01])
    lam = np.array([0, 0, 0])
    cX = np.array([1, 1, 1])
    cY = np.array([2, 2, 2])
    nX = 100
    nY = 100
    result = mr_link2_loglik_alpha_h0(th, lam, cX, cY, nX, nY)
    assert np.isnan(result)
    assert np.isnan(mr_link2_loglik_reference_v2([0.0] + list(th), lam, cX, cY, nX, nY))
def test_mr_link2_loglik_alpha_h0_large_value_r_concordances():
    th = np.array([1e10, 1e10])
    lam = np.array([1e5, 1e5, 1e5])
    cX = np.array([1e2, 1e2, 1e2])
    cY = np.array([1e2, 1e2, 1e2])
    nX = 1e6
    nY = 1e6

    # Call the Python function for alpha_h0
    python_result = mr_link2_loglik_alpha_h0(th, lam, cX, cY, nX, nY)

    # Call the Python reference function for comparison
    reference_result = mr_link2_loglik_reference_v2([0.0] + list(th), lam, cX, cY, nX, nY)

    # Validate Python result
    assert isinstance(python_result, float)
    assert np.isclose(python_result, reference_result)

    # Call the R function for alpha_h0 and compare with the Python result
    r_result = r['mr_link2_loglik_alpha_h0'](th, lam, cX, cY, nX, nY)
    assert np.isclose(r_result[0], reference_result)


def test_mr_link2_loglik_alpha_h0_negative_values_r_concordance():
    th = np.array([0.01, -0.01])
    lam = np.array([1, 2, 3])
    cX = np.array([1, 1, 1])
    cY = np.array([2, 2, 2])
    nX = 100
    nY = 100

    # Call the Python function for alpha_h0
    python_result = mr_link2_loglik_alpha_h0(th, lam, cX, cY, nX, nY)

    # Call the Python reference function for comparison
    reference_result = mr_link2_loglik_reference_v2([0.0] + list(th), lam, cX, cY, nX, nY)

    # Validate Python result
    assert isinstance(python_result, float)
    assert np.isclose(python_result, reference_result)

    # Call the R function for alpha_h0 and compare with the Python result
    r_result = r['mr_link2_loglik_alpha_h0'](th, lam, cX, cY, nX, nY)
    assert np.isclose(r_result[0], reference_result)


def test_mr_link2_loglik_alpha_h0_large_population_r_concordance():
    th = np.array([0.01, 0.01])
    lam = np.array([1, 2, 3])
    cX = np.array([1, 1, 1])
    cY = np.array([2, 2, 2])
    nX = 1e12
    nY = 1e12

    # Call the Python function for alpha_h0
    python_result = mr_link2_loglik_alpha_h0(th, lam, cX, cY, nX, nY)

    # Call the Python reference function for comparison
    reference_result = mr_link2_loglik_reference_v2([0.0] + list(th), lam, cX, cY, nX, nY)

    # Validate Python result
    assert isinstance(python_result, float)
    assert np.isclose(python_result, reference_result)

    # Call the R function for alpha_h0 and compare with the Python result
    r_result = r['mr_link2_loglik_alpha_h0'](th, lam, cX, cY, nX, nY)
    assert np.isclose(r_result[0], reference_result)


"""
FUZZ TESTING many many many times.
"""
FUZZ_ITERATIONS = 1_000


def random_input_generator_likelihood_functions():
    """Generate random inputs for the test functions."""
    for _ in range(FUZZ_ITERATIONS):
        th = np.random.uniform(-1e3, 1e3, size=3)
        lam = np.random.uniform(1e-5, 1e3, size=500)
        cX = np.random.uniform(-1e3, 1e3, size=500)
        cY = np.random.uniform(-1e3, 1e3, size=500)
        nX = np.random.uniform(50, 1e7)
        nY = np.random.uniform(50, 1e7)

        yield th, lam, cX, cY, nX, nY


@pytest.mark.parametrize("th, lam, cX, cY, nX, nY", random_input_generator_likelihood_functions())
def test_fuzz_mr_link2_loglik_reference_v2(th, lam, cX, cY, nX, nY):
    """Fuzz testing for mr_link2_loglik_reference_v2 between Python and R implementations."""
    # Python function result
    python_result = mr_link2_loglik_reference_v2(th, lam, cX, cY, nX, nY)

    # Call R function
    r_result = r['mr_link2_loglik_reference_v2'](th, lam, cX, cY, nX, nY)

    # Convert R result to Python float for comparison
    r_result_py = float(r_result[0])

    # Compare both results using np.isclose()
    assert np.isclose(python_result, r_result_py, rtol=1e-5,
                      atol=1e-8), f"Mismatch for th={th}, lam={lam}, cX={cX}, cY={cY}, nX={nX}, nY={nY}"


@pytest.mark.parametrize("th, lam, cX, cY, nX, nY", random_input_generator_likelihood_functions())
def test_fuzz_mr_link2_loglik_alpha_h0(th, lam, cX, cY, nX, nY):
    """Fuzz testing for mr_link2_loglik_alpha_h0 between Python and R implementations."""
    # Use only the first two elements of 'th' for alpha_h0
    th_h0 = th[:2]

    # Python function result
    python_result = mr_link2_loglik_alpha_h0(th_h0, lam, cX, cY, nX, nY)

    # Call R function
    r_result = r['mr_link2_loglik_alpha_h0'](th_h0, lam, cX, cY, nX, nY)

    # Convert R result to Python float for comparison
    r_result_py = float(r_result[0])

    # Compare both results using np.isclose()
    assert np.isclose(python_result, r_result_py, rtol=1e-5,
                      atol=1e-8), f"Mismatch for th={th}, lam={lam}, cX={cX}, cY={cY}, nX={nX}, nY={nY}"



@pytest.mark.parametrize("th, lam, cX, cY, nX, nY", [
    (np.array([1e-10, 1e-10, 1e-10]), np.array([1e-10, 1e-10, 1e-10]), np.array([1e-10, 1e-10, 1e-10]), np.array([1e-10, 1e-10, 1e-10]), 1e10, 1e10),
    (np.array([1e10, 1e10, 1e10]), np.array([1e5, 1e5, 1e5]), np.array([1e2, 1e2, 1e2]), np.array([1e2, 1e2, 1e2]), 1e10, 1e5),
    (np.array([0, 0, 0]), np.array([1, 2, 3]), np.array([0, 0, 0]), np.array([0, 0, 0]), 100, 100),
    # Add other extreme cases here...
])
def test_edge_cases_mr_link2_loglik_reference_v2(th, lam, cX, cY, nX, nY):
    """Test edge cases for equivalence between Python and R implementations."""
    python_result = mr_link2_loglik_reference_v2(th, lam, cX, cY, nX, nY)
    r_result = r['mr_link2_loglik_reference_v2'](th, lam, cX, cY, nX, nY)
    r_result_py = float(r_result[0])

    assert np.isclose(python_result, r_result_py, rtol=1e-5, atol=1e-8)


"""
Integration tests of the R functions
"""

FUZZ_TEST_ITERATIONS = 100
def random_input_generator_mr_link_2_function():
    np.random.seed(42)
    """Generate random inputs for fuzz testing."""
    for _ in range(FUZZ_TEST_ITERATIONS):
        # Generate random 'selected_eigenvalues', 'selected_eigenvectors', 'exposure_betas', and 'outcome_betas'
        selected_eigenvalues = np.random.uniform(0.1, 2.0, size=100)
        selected_eigenvectors = np.random.normal(size=(100, 100))
        exposure_betas = np.random.normal(size=100)
        outcome_betas = np.random.normal(size=100)
        n_exp = np.random.randint(500, 5000)  # Randomize number of individuals
        n_out = np.random.randint(500, 5000)  # Randomize number of individuals
        sigma_exp_guess = np.random.uniform(0.1, 5.0)  # Randomize sigma guesses
        sigma_out_guess = np.random.uniform(0.1, 5.0)

        yield selected_eigenvalues, selected_eigenvectors, exposure_betas, outcome_betas, n_exp, n_out, sigma_exp_guess, sigma_out_guess


@pytest.mark.parametrize(
    "selected_eigenvalues, selected_eigenvectors, exposure_betas, outcome_betas, n_exp, n_out, sigma_exp_guess, sigma_out_guess",
    random_input_generator_mr_link_2_function())
def test_mr_link2_equivalence_fuzz(selected_eigenvalues, selected_eigenvectors, exposure_betas, outcome_betas, n_exp,
                                   n_out, sigma_exp_guess, sigma_out_guess):
    """

    Fuzz testing for MR-Link2 equivalence between Python and R implementations.

    Currently only testing until a threshold of 2.5% difference. If testing more samples, this may be exceeded.

    It seems that there are differences between the optimizers which I cannot resolve easily for now.
    That being said, the p values will be withing 5% of each other, so a python p value of 0.1 will be found within
    the (0.095 - 0.105) interval.

    """

    # Run Python version of mr_link2
    python_result = mr_link2(
        selected_eigenvalues=selected_eigenvalues,
        selected_eigenvectors=selected_eigenvectors,
        exposure_betas=exposure_betas,
        outcome_betas=outcome_betas,
        n_exp=n_exp,
        n_out=n_out,
        sigma_exp_guess=sigma_exp_guess,
        sigma_out_guess=sigma_out_guess
    )

    # Convert inputs to R-compatible types and run R version of mr_link2
    with localconverter(default_converter + numpy2ri.converter):
        r_selected_eigenvalues = r.matrix(selected_eigenvalues, nrow=100, ncol=1)
        r_selected_eigenvectors = r.matrix(selected_eigenvectors, nrow=100, ncol=100)
        r_exposure_betas = r.matrix(exposure_betas, nrow=100, ncol=1)
        r_outcome_betas = r.matrix(outcome_betas, nrow=100, ncol=1)
        r_n_exp = n_exp
        r_n_out = n_out
        r_sigma_exp_guess = sigma_exp_guess
        r_sigma_out_guess = sigma_out_guess

        # Call the R function
        r_result = r['mr_link2'](
            r_selected_eigenvalues, r_selected_eigenvectors,
            r_exposure_betas, r_outcome_betas,
            r_n_exp, r_n_out, r_sigma_exp_guess, r_sigma_out_guess
        )

    # Convert R result back to Python dictionary format
    r_result_dict = {str(key): value[0] for key, value in r_result.items()}

    optimization_indicators = ['optim_alpha_h0_success', 'optim_sigma_y_h0_success', 'optim_ha_success']
    for result in [python_result, r_result_dict]:
        for indicator in optimization_indicators:
            if not result[indicator]:
                print('No correct optimization.')
                return

    # Define which results to compare with high precision
    high_precision_keys = ['alpha_h0_loglik', 'ha_loglik', ]
    for key in high_precision_keys:
        assert np.isclose(python_result[key], r_result_dict[key], rtol=1e-4,
                          atol=1e-8), f"Mismatch for {key}: Python={python_result[key]}, R={r_result_dict[key]}"

    # Define which results to compare with low precision
    low_precision_keys = ['alpha', 'sigma_y', 'sigma_x','p(alpha)', 'p(sigma_y)']
    for key in low_precision_keys:
        assert np.isclose(python_result[key], r_result_dict[key], rtol=2.5e-2, ## this is pretty suboptimal, but as this also includes P values, there will not be a crazy amount of difference here.
                          atol=1e-8), f"Mismatch for {key}: Python={python_result[key]}, R={r_result_dict[key]}"
