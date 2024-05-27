#!/usr/bin/env python
"""
Calculates the exponent for gamma-distributed latent and infectious periods. 
"""
import sys
sys.path.insert(0,"../modules")

from math import isclose
from scipy.optimize import fsolve

# Laplace transform of latent period gamma(mu, kappa)
def _laplace_exposed(theta, mu, kappa):
    laplace = (1.0+theta*mu/kappa)**(-kappa)
    return laplace
# Laplace transform of uniform distribution on infectious period gamma(mu, kappa)
def _laplace_infectious(theta, mu, kappa):
    if theta <= 0.0:
        return 1.0
    laplace = 1.0-(1.0+theta*mu/kappa)**(-(kappa-1.0))
    laplace /= kappa-1.0
    laplace /= theta*mu/kappa
    return laplace
# Laplace transform of the random generation time
def theta_solve(e_mu, e_kappa, i_mu, i_kappa, r0):    
    def laplace_generation(theta):
        return _laplace_exposed(theta, e_mu, e_kappa)*_laplace_infectious(theta, i_mu, i_kappa)-1.0/r0
    theta0 = 0.1
    theta = fsolve(laplace_generation, theta0) # exponential rate (lambda)
    assert isclose(laplace_generation(theta), 0.0, abs_tol=1.0e-09)
    return theta
# Reads the string from test data, with columns COLS.
def test_theta_solve():
    (e_mu, e_kappa, i_mu, i_kappa, r0) = (3.5, 4, 5.5, 0.3, 2)
    theta = theta_solve(e_mu, e_kappa, i_mu, i_kappa, r0)
    assert isclose(theta, 0.14156035, abs_tol=1.0e-06)
    assert isclose(_laplace_exposed(theta, e_mu, e_kappa), 0.62682015, abs_tol=1.0e-06)
    assert isclose(_laplace_infectious(theta, i_mu, i_kappa), 0.79767698, abs_tol=1.0e-06)
    assert isclose(1.0/r0, 0.5, abs_tol=1.0e-06)

if __name__ == "__main__":
    test_theta_solve()
