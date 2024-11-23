#!/usr/bin/env python
"""
Detecting if a timeseries is stationary. 
In practice a simple filter implemented in detectScale followed by computing the statistical inneficiency ( i.e. correlation time delta ) avoids issues with spurious random peaks shifting the automatically detected region past them.

The below applies in part to statisticalInefficiency_multiscale.

Utilities for automatically detecting equilibrated region of molecular simulations.

John D. Chodera <john.chodera@choderalab.org>
Sloan Kettering Institute
Memorial Sloan Kettering Cancer Center

LICENSE

This is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation, either version 2.1
of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this software. If not, see <http://www.gnu.org/licenses/>.

"""
from __future__ import annotations

import numpy as np
from numba import njit
from uncertainties import ufloat


@njit  # type:ignore
def _multiscale_nb(dA_n, dB_n, N: int, sigma2_AB: float, mintime: int, fast: bool):
    """
    Accumulate the integrated correlation time by computing the normalized correlation time at increasing values of t.  Stop accumulating if the correlation function goes negative, since this is unlikely to occur unless the correlation function has decayed to the point where it is dominated by noise and indistinguishable from zero.
    """
    # Initialize statistical inefficiency estimate with uncorrelated value.
    g = 1.0
    t = 1
    increment = 1
    while t < N - 1:
        # compute normalized fluctuation correlation function at time t
        C: float = np.sum(
            dA_n[0 : (N - t)] * dB_n[t:N] + dB_n[0 : (N - t)] * dA_n[t:N]
        ) / (2 * (N - t) * sigma2_AB)
        # Terminate if the correlation function has crossed zero and we've computed the correlation
        # function at least out to 'mintime'.
        if (C <= 0.0) and (t > mintime):  # filters out noise
            break

        # Accumulate contribution to the statistical inefficiency.
        g += 2.0 * C * (1.0 - t / N) * increment

        # Increment t and the amount by which we increment t.
        t += increment

        # Increase the interval if "fast mode" is on.
        if fast:
            increment += 1

    # g must be at least unity
    if g < 1.0:
        g = 1.0

    # Return the computed statistical inefficiency.
    return g


def statisticalInefficiency_multiscale(
    A_n: np.ndarray,
    B_n: np.ndarray | None = None,
    fast: bool = True,
    mintime: int = 3,
):
    """
    Compute the (cross) statistical inefficiency of (two) timeseries using multiscale method from Chodera.

    Parameters
    ----------
    A_n : np.ndarray, float
        A_n[n] is nth value of timeseries A.  Length is deduced from vector.
    B_n : np.ndarray, float, optional, default=None
        B_n[n] is nth value of timeseries B.  Length is deduced from vector.
        If supplied, the cross-correlation of timeseries A and B will be estimated instead of the
        autocorrelation of timeseries A.
    fast : bool, optional, default=False
        f True, will use faster (but less accurate) method to estimate correlation
        time, described in Ref. [1] (default: False).
    mintime : int, optional, default=3
        minimum amount of correlation function to compute (default: 3)
        The algorithm terminates after computing the correlation time out to mintime when the
        correlation function first goes negative.  Note that this time may need to be increased
        if there is a strong initial negative peak in the correlation function.

    Returns
    -------
    g : np.ndarray,
        g is the estimated statistical inefficiency (equal to 1 + 2 tau, where tau is the correlation time).
        We enforce g >= 1.0.

    Notes
    -----
    The same timeseries can be used for both A_n and B_n to get the autocorrelation statistical inefficiency.
    The fast method described in Ref [1] is used to compute g.

    References
    ----------
    [1] J. D. Chodera, W. C. Swope, J. W. Pitera, C. Seok, and K. A. Dill. Use of the weighted
        histogram analysis method for the analysis of simulated and parallel tempering simulations.
        JCTC 3(1):26-41, 2007.

    See also: https://doi.org/10.1021/acs.jctc.5b00784

    Examples
    --------

    Compute statistical inefficiency of timeseries data with known correlation time.

    >>> from pymbar import testsystems
    >>> A_n = testsystems.correlated_timeseries_example(N=100000, tau=5.0)
    >>> g = statisticalInefficiency_multiscale(A_n, fast=True)

    """

    A_n = np.array(A_n)

    if B_n is not None:
        B_n = np.array(B_n)
    else:
        B_n = np.array(A_n)

    # Get the length of the timeseries.
    N = A_n.size

    # Be sure A_n and B_n have the same dimensions.
    if A_n.shape != B_n.shape:
        raise Exception("A_n and B_n must have same dimensions.")

    # Compute mean of each timeseries.
    mu_A = A_n.mean()
    mu_B = B_n.mean()

    # Make temporary copies of fluctuation from mean.
    dA_n = A_n.astype(np.float64) - mu_A  # type: ignore
    dB_n = B_n.astype(np.float64) - mu_B  # type: ignore

    # Compute estimator of covariance of (A,B) using estimator that will ensure C(0) = 1.
    sigma2_AB = (dA_n * dB_n).mean()  # standard estimator to ensure C(0) = 1

    # Trap the case where this covariance is zero, and we cannot proceed.
    if sigma2_AB == 0:
        return None
        # raise RuntimeError('Sample covariance sigma_AB^2 = 0 -- cannot compute statistical inefficiency')
    return _multiscale_nb(dA_n, dB_n, N, sigma2_AB, mintime, fast)


def series_ineff_2_ufloat(x: np.ndarray, ineff: float):
    return ufloat(x.mean(), x.std() * (ineff / len(x)) ** 0.5)


def series2ufloat(x: np.ndarray):
    assert not np.isnan(x).any()
    stat_ineff = statisticalInefficiency_multiscale(x)
    if stat_ineff is None:
        return None
    return series_ineff_2_ufloat(x, stat_ineff)
