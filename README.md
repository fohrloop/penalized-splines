# Monotone Smoothing B-Splines

Monotone smoothing splines implementation using penalized B-splines (aka. P-splines).

Solves the equation

    (B'B + λD3'D3 + κD1'VD1)α = B'y

where

    B'B: The least squares part
    λD3'D3: The smoothing part
    κD1'VD1: The monotonicity part
    α: The coefficients of the B-spline basis functions

The algorithm was introduced in [^Eilers2005], and explained for example in
[^Hofner2011] and [^Hofner2014unified]. It is iterative at nature and gives *approximately*
monotone spline function as output which are "monotone enough" for most of the practical
applications. The monotonicity is enforced by a large penalty term κ. If you require
strictly monotone function, you should check I-splines. Sinthe this uses "regular"
P-splines, you might use any other added penalty terms to enforce other properties, like
cyclicity. The monotonicity part can be swapped to enforce convex/concave shapes
(U-shape or inverse U-shape).

[^Eilers2005]: Eilers, P. H. C. 2005. Unimodal smoothing. Journal of Chemometrics
              19:317–328. DOI:10.1002/cem.935.
[^Hofner2011]: Hofner et al. (2011) Monotonicity-constrained species distribution models.
            Link: https://esajournals.onlinelibrary.wiley.com/doi/pdf/10.1890/10-2276.1
[^Hofner2014unified]: Hofner et al (2014) A Unified Framework of Constrained Regression.
            Link: https://arxiv.org/abs/1403.7118