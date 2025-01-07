"""Monotone smoothing splines implementation using penalized B-splines (aka. P-splines).

Solves the equation

    (B'B + λD3'D3 + κD1'VD1)α = B'y

where

    B'B: The least squares part
    λD3'D3: The smoothing part
    κD1'VD1: The monotonicity part
    α: The coefficients of the B-spline basis functions

The algorithm was introduced in [Eilers2005]

References:
-----------
[Eilers2005]: Eilers, P. H. C. 2005. Unimodal smoothing. Journal of Chemometrics
              19:317–328. DOI:10.1002/cem.935.
"""

from __future__ import annotations

import warnings

import numpy as np
from scipy.interpolate import BSpline


class MaxIterationWarning(UserWarning): ...


def create_pspline(
    x_train,
    y_train,
    bspline_degree=3,
    knot_segments=10,
    lambda_smoothing=0.1,
    kappa_penalty=10**6,
    maxiter: int = 30,
) -> BSpline:
    """
    Parameters
    ----------
    bspline_degree: int
        The degree of the B-spline (which is also the degree of the fitted spline
        function). The order of the splines is degree + 1.
    knot_segments: int
        number of inter-knot segments between min(x) and max(x). Defines the number of
        knots. Will be limited to len(x)-2. (giving too high number will result in
        len(x)-2).
    lambda_smoothing: float
        The smoothing parameter. Higher values will result in smoother curves.
    kappa_penalty: float
        The penalty parameter for enforcing monotonicity. Higher values will result in
        more monotonic curves.
    maxiter: int
        Maximum number of iterations for the algorithm. If the algorithm does not
        converge within this number of iterations, a warning is issued.
    """

    knot_interval = (max(x_train) - min(x_train)) / knot_segments

    # You need to add deg knots on each side of the interval. See, for example,
    # De Leeuw (2017) Computing and Fitting Monotone Splines
    # The basic interval is [min(x_train), max(x_train)], and
    # the extended interval is [min(knots), max(knots)].
    # You may only ask for values within the basic interval, as there are always m
    # (=deg+1) non-zero B-splines. Outside the basic interval, there are less B-splines
    # with non-zero values and the model is extrapolating.
    knots = np.linspace(
        min(x_train) - (bspline_degree + 1) * knot_interval,
        max(x_train) + (bspline_degree + 1) * knot_interval,
        bspline_degree * 2 + knot_segments + 1,
    )
    alphas = np.ones(len(x_train))

    # Introduced in scipy 1.8.0
    B = BSpline.design_matrix(x=x_train, t=knots, k=bspline_degree).toarray()
    n_base_funcs = B.shape[1]
    I = np.eye(n_base_funcs)
    D3 = np.diff(I, n=3, axis=0)
    D1 = np.diff(I, n=1, axis=0)

    # Monotone smoothing
    V = np.zeros(n_base_funcs - 1)

    B_gram = B.T @ B
    A = B_gram + lambda_smoothing * D3.T @ D3
    BTy = B.T @ y_train
    for _ in range(maxiter):

        W = np.diag(V * kappa_penalty)
        # The equation
        # (B'B + λD3'D3 + κD1'VD1)α = B'y
        alphas = np.linalg.solve(A + D1.T @ W @ D1, BTy)
        V_new = (D1 @ alphas < 0) * 1
        dv = np.sum(V != V_new)
        V = V_new
        if dv == 0:
            break
    else:
        warnings.warn(
            "Max iteration reached. The results are not reliable.", MaxIterationWarning
        )

    return BSpline(knots, alphas, bspline_degree, extrapolate=False)
