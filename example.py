import numpy as np
from matplotlib import pyplot as plt

from monotone_splines import create_pspline

y = np.array([1.0, 2.0, 2.5, 3.4, 3.0, 3.6, 3.33, 3.0])
x = np.array([1, 8, 15, 22, 30, 38, 46, 54])
x_valid = np.linspace(min(x), max(x), 1000)

knot_segments = 100

monotone_spline = create_pspline(
    x, y, lambda_smoothing=100000.0, knot_segments=knot_segments, kappa_penalty=10e6
)

plt.plot(x, y, label="raw values", marker="o", ls="", color="gray")
plt.plot(x_valid, monotone_spline(x_valid), label="spline", ls="--", color="tab:red")
plt.plot(
    x, monotone_spline(x), label="spline (train)", marker=".", ls="", color="tab:red"
)

plt.legend()
plt.grid(ls="--", lw=0.5, color="lightgray")
plt.show()
