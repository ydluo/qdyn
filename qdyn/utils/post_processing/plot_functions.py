"""
Useful functions for post-processing (plotting) 
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation, rc
from IPython.display import HTML
import matplotlib.cm as cm


t_day = 3600 * 24.0
t_year = 365 * t_day
figsize = (9, 8)


def timeseries(ot, ot_vmax):

    fig, axes = plt.subplots(nrows=4, ncols=1, figsize=figsize)

    axes[0].plot(ot["t"] / t_year, ot["tau"])
    axes[0].set_ylabel("tau [Pa]")

    axes[1].plot(ot["t"] / t_year, ot["theta"])
    axes[1].set_ylabel("state [s]")

    axes[2].plot(ot["t"] / t_year, ot["sigma"])
    axes[2].set_ylabel("sigma [Pa]")

    axes[3].plot(ot_vmax["t"] / t_year, ot_vmax["v"])
    axes[3].set_ylabel("max v [m/s]")
    axes[3].set_xlabel("time [yr]")
    axes[3].set_yscale("log")

    plt.tight_layout()
    plt.show()


def slip_profile(ox, warm_up=0, orientation="horizontal"):
    x_unique = ox["x"].unique()
    sort_inds = np.argsort(x_unique)
    x_unique = x_unique[sort_inds]
    t_vals = np.sort(ox["t"].unique())

    if warm_up > t_vals.max():
        print("Warm-up time > simulation time!")
        return

    ind_warmup = np.where(t_vals >= warm_up)[0][0]

    print(ind_warmup)

    Nx = len(x_unique)
    Nt = len(t_vals) - 1

    slice = np.s_[Nx * ind_warmup:Nx * Nt]
    data_shape = (Nt - ind_warmup, Nx)
    x = ox["x"][slice].values.reshape(data_shape)[:, sort_inds]
    slip = ox["slip"][slice].values.reshape(data_shape)[:, sort_inds]
    v = ox["v"][slice].values.reshape(data_shape)[:, sort_inds]

    slip -= slip[0]
    t = t_vals[ind_warmup:-1]
    t -= t[0]

    fig = plt.figure(figsize=figsize)

    if orientation == "horizontal":
        CS = plt.contourf(x, slip, np.log10(v), levels=200, cmap="magma")
        plt.xlabel("position [m]")
        plt.ylabel("slip [m]")
        CB = plt.colorbar(CS, orientation="horizontal")
        CB.ax.set_title("slip rate [m/s]")
    elif orientation == "vertical":
        x -= x.min()
        CS = plt.contourf(slip, x * 1e-3, np.log10(v), levels=200, cmap="magma")
        plt.ylabel("depth [km]")
        plt.xlabel("slip [m]")
        plt.gca().invert_yaxis()
        CB = plt.colorbar(CS, orientation="horizontal")
        CB.ax.set_title("slip rate [m/s]")
    else:
        print("Keyword 'orientation=%s' not recognised" % orientation)
        plt.close()
        return
    plt.tight_layout()
    plt.show()

def slip_profile_new(ox, dip, warm_up=0, orientation="horizontal"):
    """
    Plot slip profile along the fault trace.

    Parameters
    ----------
    ox : pandas DataFrame
        DataFrame containing ox output. 
        If fault is 0D or 1D --> ox=p.ox
        If fault is 2D, the x or z coordinate has to be provided
        --> for vertical cross-section: ox = p.ox[(p.ox["x"] == x_pos)] where x_pos is the x coordinate of the cross-section
        --> for horizontal cross-section: ox = p.ox[(p.ox["z"] == z_pos)] where z_pos is the z coordinate of the cross-section
        If there are multiple faults, the fault label has to be provided
        --> ox[(ox["fault_label"] == 1)]
        Example for slip profile of a vertical cross-section of Fault 1:
        p.ox[(p.ox["fault_label"] == 1) & (p.ox["x"] == x_pos)]
    dip : float
        Dip angle in degrees
    warm_up : float, optional
        Time in seconds to disregard in the plot (default is 0)
    orientation : str, optional
        Plot/cross section orientation. Can be 'horizontal' (default) or 'vertical'

    """
    
    # Check orientation parameter
    if orientation not in ["horizontal", "vertical"]:
        print("Error: orientation parameter must be either 'horizontal' or 'vertical'")
        return
    
    # Determine coordinate to plot along
    if orientation=="horizontal":
        coord="x"
    else:
        coord="z"
    
    # Get unique coordinates and sort them
    x_unique = ox[coord].unique()
    sort_inds = np.argsort(x_unique)
    x_unique = x_unique[sort_inds]
    
    # Get unique times and sort them
    t_vals = np.sort(ox["t"].unique())
    
    # Check that warm-up time is not greater than simulation time
    if warm_up > t_vals.max():
        print("Warm-up time > simulation time!")
        return
    
    # Determine index for warm-up time
    ind_warmup = np.where(t_vals >= warm_up)[0][0]
    
    # Get the number of unique coordinates and number of time steps
    Nx = len(x_unique)
    Nt = len(t_vals) - 1
    
    # Define the slice and shape for the data
    slice = np.s_[Nx * ind_warmup:Nx * Nt]
    data_shape = (Nt - ind_warmup, Nx)
    
    # Get the data values and reshape them
    x = ox[coord][slice].values.reshape(data_shape)[:, sort_inds]
    slip = ox["slip"][slice].values.reshape(data_shape)[:, sort_inds]
    v = ox["v"][slice].values.reshape(data_shape)[:, sort_inds]
    
    # Adjust slip values to start at zero
    slip -= slip[0]
    
    # Adjust time values to start at zero
    t = t_vals[ind_warmup:-1]
    t -= t[0]

    # Calculate downdip distance for vertical profile
    if orientation=="vertical":
        ddd = x / np.sin(np.deg2rad(dip))
    
    # Create figure and plot the slip profile
    fig = plt.figure(figsize=figsize)

    if orientation == "horizontal":
        CS = plt.contourf(x, slip, np.log10(v), levels=200, cmap="magma")
        plt.xlabel("position [m]")
        plt.ylabel("slip [m]")
        CB = plt.colorbar(CS, orientation="horizontal")
        CB.ax.set_title("slip rate [m/s]")
    elif orientation == "vertical":
        ddd -= ddd.min()
        CS = plt.contourf(slip, ddd * 1e-3, np.log10(v), levels=200, cmap="magma")
        plt.ylabel("downdip distance [km]")
        plt.xlabel("slip [m]")
        plt.gca().invert_yaxis()
        CB = plt.colorbar(CS, orientation="horizontal")
        CB.ax.set_title("slip rate [m/s]")
    else:
        print("Keyword 'orientation=%s' not recognised" % orientation)
        plt.close()
        return

    plt.tight_layout()
    plt.show()

def animation_slip(ox, warm_up=0, orientation="horizontal"):

    x_unique = ox["x"].unique()
    sort_inds = np.argsort(x_unique)
    x_unique = x_unique[sort_inds]
    t_vals = np.sort(ox["t"].unique())

    if warm_up > t_vals.max():
        print("Warm-up time > simulation time!")
        return

    ind_warmup = np.where(t_vals >= warm_up)[0][0]

    Nx = len(x_unique)
    Nt = len(t_vals) - 1

    slice = np.s_[Nx * ind_warmup:Nx * Nt]
    data_shape = (Nt - ind_warmup, Nx)
    slip = ox["slip"][slice].values.reshape(data_shape)[:, sort_inds]
    v = ox["v"][slice].values.reshape(data_shape)[:, sort_inds]

    slip -= slip[0]
    t = t_vals[ind_warmup:-1]
    t -= t[0]

    plt.ioff()

    if orientation == "horizontal":
        fig, ax = plt.subplots(nrows=2, ncols=1, figsize=figsize)

        ax[0].set_xlim((x_unique.min(), x_unique.max()))
        ax[0].set_ylim((slip.min(), slip.max()))
        ax[0].set_ylabel("slip [m]")

        ax[1].set_xlim((x_unique.min(), x_unique.max()))
        ax[1].set_ylim((v.min(), v.max()))
        ax[1].set_yscale("log")
        ax[1].set_ylabel("slip rate [m/s]")
        ax[1].set_xlabel("position [m]")

        line1, = ax[0].plot([], [], lw=2)
        line2, = ax[1].plot([], [], lw=2)
        lines = (line1, line2)

        def init_plot():
            for line in lines:
                line.set_data([], [])
            return lines,

        def animate(i):
            line1.set_data(x_unique, slip[i])
            line2.set_data(x_unique, v[i])

    elif orientation == "vertical":
        x_unique -= x_unique.min()
        x_unique *= 1e-3
        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=figsize)

        ax[0].set_ylim((x_unique.min(), x_unique.max()))
        ax[0].set_xlim((slip.min(), slip.max()))
        ax[0].set_ylabel("depth [km]")
        ax[0].set_xlabel("slip [m]")
        ax[0].invert_yaxis()
        ax[0].xaxis.tick_top()
        ax[0].xaxis.set_label_position("top")

        ax[1].set_ylim((x_unique.min(), x_unique.max()))
        ax[1].set_xlim((v.min(), v.max()))
        ax[1].set_xscale("log")
        ax[1].set_xlabel("slip rate [m/s]")
        ax[1].invert_yaxis()
        ax[1].xaxis.tick_top()
        ax[1].xaxis.set_label_position("top")

        line1, = ax[0].plot([], [], lw=2)
        line2, = ax[1].plot([], [], lw=2)
        lines = (line1, line2)

        def init_plot():
            for line in lines:
                line.set_data([], [])
            return lines,

        def animate(i):
            line1.set_data(slip[i], x_unique)
            line2.set_data(v[i], x_unique)

    anim = animation.FuncAnimation(fig, animate, init_func=init_plot,
                                   frames=len(t), interval=30, blit=True)
    plt.tight_layout()
    HTML(anim.to_html5_video())
    rc("animation", html="html5")
    plt.close()
    return anim

def plot_snapshot_3d(ox, set_dict, t_snapshot=0, prop="v", fault_labels = [1], scaling="Normalize"):
    
    """
    Plot a 3D snapshot.

    Parameters:
        ox (DataFrame): p.ox containing the snapshots.
        set_dict (dict): p.set_dict containing settings and properties.
        t_snapshot (float): The time value (s) for the desired snapshot (default is 0s).
        prop (str): The property to visualize (default is "v").
        fault_labels (list): A list of fault labels to filter the data (default is [1]).
        scaling (str): The scaling method for the colormap (default is "Normalize").

    Raises:
        ValueError: If an invalid list of fault labels or scaling option is provided.

    Returns:
        None: This function displays the 3D plot using Matplotlib.
        
    Note:
        It is possible to plot dsigma (sigma_0 - sigma) even if it's not originally included
        in the ox DataFrame, since the calculation is internally handled by this function.

    Example:
        To plot a 3D snapshot of fault label 2 with property "dsigma" using a normalized colormap:
        >>> plot_snapshot_3d(p.ox, p.set_dict, t_snapshot=1.0, fault_labels=[2], prop="dsigma", scaling="Normalize")
    """
    
    # Handle errors
    unique_labels = np.unique(ox["fault_label"])
    if not all(label in unique_labels for label in fault_labels):
        raise ValueError("Invalid list of fault labels")
                              
    # Index of elements of snapshot
    inds = (ox["t"] == t_snapshot) & (ox["fault_label"].isin(fault_labels))

    # filter snapshot
    mesh_snapshot = ox[inds]

    # Draw canvas
    plt.close("all")
    fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True, subplot_kw={"projection": "3d"})

    cmap = cm.magma

    # Decorate axes
    ax.set_xlabel("x [km]")
    ax.set_ylabel("y [km]")
    ax.set_zlabel("depth [km]")
    ax.set_title("t =" + str(t_snapshot/t_year) + " yr\n" + "fault " + str(fault_labels), pad = 20)

    # Set initial viewing angle
    ax.set_aspect("equal")
    ax.view_init(elev=51, azim=-60)

    # Select quantities to plot
    x = mesh_snapshot["x"]
    y = mesh_snapshot["y"]
    z = mesh_snapshot["z"]
    if prop == "dsigma":
        sigma_0 = set_dict["SIGMA"]
        col = sigma_0 - mesh_snapshot["sigma"]
    elif prop == "v":
        col = np.log10(mesh_snapshot[prop])
    else:
        col = mesh_snapshot[prop]

        
    # Colour scale normalisation
    if prop == "dsigma":
        cmap = cm.bwr
    else:
        cmap = cm.magma
    
    if scaling == "Normalize":
        norm = cm.colors.Normalize(col.min(), col.max())
    else: 
        raise ValueError("Invalid colormap scaling option. Choose from: Normalized")

    # Plot snapshot
    sc = ax.scatter(x, y, z,
                    c=col, norm=norm, cmap=cmap, s=10)
    # colorbar
    cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.5)
    cbar.set_label(str(prop), rotation=90, labelpad=20)
    if prop == "v":
        cbar.set_label("logv", rotation=90, labelpad=20)
    
    return fig