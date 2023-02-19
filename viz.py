from matplotlib import cm, pyplot as plt


def subplots_3d(nrows=1, ncols=1, **fig_kw):
    return plt.subplots(nrows, ncols, subplot_kw={"projection": "3d"}, **fig_kw)


def plot_solution(ax, X, Y, u, txt='Solution'):
    # Plot the solution of the heat equation

    ax.plot_surface(X, Y, u, cmap=cm.coolwarm)
    ax.view_init(azim=30)  # Rotate the figure
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(txt)
    return plt
