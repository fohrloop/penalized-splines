import numpy as np

try:
    from matplotlib import pyplot as plt
except ImportError:
    print(
        "Please install matplotlib to run this example:\n\n    pip install matplotlib"
    )
    exit()


from monotone_splines import create_monotone_bspline

if __name__ == "__main__":
    y = np.array([1.0, 2.0, 2.5, 3.4, 3.0, 3.6, 3.33, 3.0])
    x = np.array([1, 8, 15, 22, 30, 38, 46, 54])
    x_valid = np.linspace(min(x), max(x), 1000)
    spline = create_monotone_bspline(x, y, lambda_smoothing=100000.0, knot_segments=120)

    plt.plot(x, y, label="raw values", marker="o", ls="", color="lightgray")
    plt.plot(x_valid, spline(x_valid), label="spline", ls="--", color="red")
    plt.plot(x, spline(x), label="spline (train)", marker="o", ls="", color="red")
    plt.legend()
    plt.grid(ls="--", lw=0.5, color="lightgray")
    plt.show()
