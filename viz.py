import numpy as np
from matplotlib import cm, pyplot as plt


def subplots_3d(nrows=1, ncols=1, **fig_kw):
    return plt.subplots(nrows, ncols, subplot_kw={"projection": "3d"}, **fig_kw)


def plot_solution(X, Y, u, ax=None, txt='Solution'):
    # Plot the solution of the stationary heat equation

    if ax is None:
        _, ax = subplots_3d()

    ax.plot_surface(X, Y, u, cmap=cm.coolwarm)
    ax.view_init(azim=30)  # Rotate the figure
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(txt)
    return plt


def loglogplot_error(h, e, ax=None):

    if ax is None:
        _, ax = plt.subplots()

    p = np.polyfit(np.log(h), np.log(e), 1)[0]

    ax.set_title("Error")
    ax.loglog(h, e, label=f"p={p:.3f}")
    ax.set_xlabel("h")
    ax.set_ylabel("e")
    ax.invert_xaxis()
    plt.legend()
    return plt

