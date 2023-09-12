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

def slip_profile_new(ox, dip, warm_up=0, orientation="horizontal", figsize=figsize):
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
    figsize : tuple, optional
        Figure size. Default (9,8)

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



def timestep_profile(ox, dip, warm_up=0, orientation="horizontal", figsize=figsize):
    """
    Plot timestep profile along the fault trace with slip rate cotouring.

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
    figsize : tuple, optional
        Figure size. Default (9,8)
    
    
    Example
    -------
    To plot an horizontal profile in the middle of a 2D dipping fault with label 1 from a model with two faults:
    >>>> timestep_profile(p.ox[(p.ox["z"] == z_pos) &  (p.ox["fault_label"] == 1)], dip=60, warm_up=50*t_yr, orientation="horizontal")

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
    timestep = ox["step"][slice].values.reshape(data_shape)[:, sort_inds]
    v = ox["v"][slice].values.reshape(data_shape)[:, sort_inds]
    
    # max and min values
    
    
    # Adjust timestep values to start at zero
    # slip -= slip[0]
    timestep -= timestep[0]
    
    # Adjust time values to start at zero
    t = t_vals[ind_warmup:-1]
    t -= t[0]

    # Calculate downdip distance for vertical profile
    if orientation=="vertical":
        ddd = x / np.sin(np.deg2rad(dip))
    
    # Create figure and plot the slip profile
    fig = plt.figure(figsize=figsize)

    if orientation == "horizontal":
        CS = plt.contourf(timestep, x, np.log10(v), levels=100, cmap="magma")
        plt.ylabel("position [m]")
        # plt.ylabel("slip [m]")
        plt.xlabel("Timestep")
        CB = plt.colorbar(CS, orientation="horizontal")
        CB.ax.set_title("slip rate [m/s]")
    elif orientation == "vertical":
        ddd -= ddd.min()
        CS = plt.contourf(ddd * 1e-3, timestep, np.log10(v), levels=100, cmap="magma")
        plt.xlabel("downdip distance [km]")
        # plt.xlabel("slip [m]")
        plt.ylabel("Timestep")
        plt.gca().invert_yaxis()
        CB = plt.colorbar(CS, orientation="horizontal")
        CB.ax.set_title("slip rate [m/s]")
    else:
        print("Keyword 'orientation=%s' not recognised" % orientation)
        plt.close()

    plt.tight_layout()
    plt.show()
        
    return fig

# CRP: Experimental!
def timestep_profile_all(ox, dip, val, sigma_0=None, warm_up=0, orientation="horizontal", figsize=figsize):
    
    """
    Plot timestep profile along the fault trace with a value from a col of ox as contouring.
    
    Additional libraries
    --------
    cmcrameri

    Parameters
    ----------
    ox : pandas DataFrame
        DataFrame containing ox output.
        If the fault is 0D or 1D --> ox = p.ox
        If the fault is 2D, the x or z coordinate has to be provided:
        - For a vertical cross-section: ox = p.ox[(p.ox["x"] == x_pos)], where x_pos is the x coordinate of the cross-section.
        - For a horizontal cross-section: ox = p.ox[(p.ox["z"] == z_pos)], where z_pos is the z coordinate of the cross-section.
        If there are multiple faults, the fault label has to be provided:
        - ox[(ox["fault_label"] == 1)]
        Example for plotting a slip profile of a vertical cross-section of Fault 1:
        p.ox[(p.ox["fault_label"] == 1) & (p.ox["x"] == x_pos)]
    dip : float
        Dip angle in degrees.
    val : str
        The value to plot. It can be one of the following:
        - "v" for slip rate.
        - "dsigma" for change in stress (requires sigma_0 parameter).
        - "theta" for the angle parameter.
        - "tau" for shear stress.
        - "tau_dot" for shear stress rate.
        - "sigma" for normal stress.
        - "tau_sigma" for the ratio of shear stress to normal stress.
        - "slip" for cumulative slip.
    sigma_0 : float, optional
        Initial stress value. Required if val is "dsigma". Default is None.
    warm_up : float, optional
        Time in seconds to disregard in the plot. Default is 0.
    orientation : str, optional
        Plot/cross-section orientation. Can be "horizontal" (default) or "vertical" (for 2D fault).
        Note that it will raise an error if trying to plot a "vertical" cross section in a fault other than in 2D.
    figsize : tuple, optional
        Figure size. Default is (9, 8).

    Examples
    --------
    To plot a horizontal profile in the middle of a 2D dipping fault with label 1 from a model with two faults:
    >>>> timestep_profile_all(p.ox[(p.ox["z"] == z_pos) & (p.ox["fault_label"] == 1)], dip=60, val="v", warm_up=50*t_yr, orientation="horizontal")

    """
    import matplotlib as mpl
    import cmcrameri.cm as cmc
    
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
        raise ValueError("Warm-up time > simulation time!")
    
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
    timestep = ox["step"][slice].values.reshape(data_shape)[:, sort_inds]
    
    if val=="v":
        v = ox["v"][slice].values.reshape(data_shape)[:, sort_inds]
    elif val=="dsigma":
        # check sigma
        if sigma_0 ==None:
            raise ValueError("You have to set the value for sigma_0")         
        dsigma = sigma_0 - ox["sigma"][slice].values.reshape(data_shape)[:, sort_inds]      
    elif val=="theta":
        theta = ox["theta"][slice].values.reshape(data_shape)[:, sort_inds]
    elif val=="tau":
        tau = ox["tau"][slice].values.reshape(data_shape)[:, sort_inds]
    elif val=="tau_dot":
        tau_dot = ox["tau_dot"][slice].values.reshape(data_shape)[:, sort_inds]
    elif val=="sigma":
        sigma = ox["sigma"][slice].values.reshape(data_shape)[:, sort_inds]
    elif val=="tau_sigma":
        # check sigma
        if sigma_0 ==None:
            raise ValueError("You have to set the value for sigma_0")           
        tau_sigma = ox["tau"][slice].values.reshape(data_shape)[:, sort_inds] / sigma_0
    elif val=="slip":
        slip = ox["slip"][slice].values.reshape(data_shape)[:, sort_inds]
    else:
        raise ValueError("val not included in ox")
        
    # Adjust timestep values to start at zero
    # slip -= slip[0]
    timestep -= timestep[0]
    
    # Adjust time values to start at zero
    t = t_vals[ind_warmup:-1]
    t -= t[0]

    # Calculate downdip distance for vertical profile
    if orientation=="vertical":
        ddd = x / np.sin(np.deg2rad(dip))
    
    # Create figure and plot the slip profile
    fig = plt.figure(figsize=figsize)
    
    # Special colorscales
    if val == "v":
        # make a colormap that has creep and dyn clearly delineated
        # make a colormap that has land and ocean clearly delineated and of the
        # same length (256 + 256)
        colors_creep = cmc.batlowW(np.linspace(0, 0.4, 128))
        colors_dyn = cmc.batlowW(np.linspace(0.6, 1, 128))
        all_colors = np.vstack((colors_creep, colors_dyn))
        vcmap = mpl.colors.LinearSegmentedColormap.from_list(
            'vc_map', all_colors)

        # norm
        vmin = np.log10(1e-14)
        vmax = np.log10(1e1)
        vnorm = mpl.colors.TwoSlopeNorm(vmin=vmin, vcenter=np.log10(1e-2), vmax=vmax)

        #ticks (adjust number of elements according to vmin and vmax)
        ticks = np.linspace(vmin, vmax, 16, endpoint=True)

        # Define the contour levels with a break at v=1e-2
        levels = [-np.inf, np.log10(1e-2), np.inf]
        
    
    elif val == "dsigma":
        vcmap = cmc.vik
        vnorm = mpl.colors.CenteredNorm()
        
        #ticks (adjust number of elements according to vmin and vmax) --> not used for now
        dsigma_ox = sigma_0 - ox["sigma"]
        
        vmax = np.max([abs(max(dsigma_ox)), abs(max(dsigma_ox))])
        vmin = -vmax
        
        ticks = np.linspace(vmin, vmax, 10, endpoint=True)
    
    
    if orientation == "horizontal":
        if val == "v":
            CS = plt.contourf(timestep, x, np.log10(v), levels=200, cmap=vcmap, norm=vnorm)
            CB = plt.colorbar(plt.cm.ScalarMappable(cmap=vcmap, norm=vnorm), orientation="horizontal", extend='both', ticks=ticks, pad=0.2)
            CB.ax.set_title("log v [m/s]")
        elif val=="dsigma":
            CS = plt.contourf(timestep, x, dsigma/1e6, levels=200, cmap=vcmap, norm=vnorm)
            CB = plt.colorbar(CS, orientation="horizontal", pad=0.2)
            CB.ax.set_title("dsigma [MPa]")
        elif val=="sigma":
            CS = plt.contourf(timestep, x, sigma/1e6, levels=200, cmap=cmc.batlow)
            CB = plt.colorbar(CS, orientation="horizontal", pad=0.2)
            CB.ax.set_title("sigma [MPa]")        
        elif val=="tau":
            CS = plt.contourf(timestep, x, tau/1e6, levels=200, cmap=cmc.batlow)
            CB = plt.colorbar(CS, orientation="horizontal", pad=0.2)
            CB.ax.set_title("tau [MPa]")
        elif val == "tau_sigma":
            CS = plt.contourf(timestep, x, tau_sigma, levels=200, cmap=cmc.lipari)
            CB = plt.colorbar(CS, orientation="horizontal", pad=0.2)
            CB.ax.set_title("tau/sigma")
        elif val=="theta":
            CS = plt.contourf(timestep, x, theta, levels=200, cmap=cmc.lipari)
            CB = plt.colorbar(CS, orientation="horizontal", pad=0.2)
            CB.ax.set_title("theta")
        elif val == "tau_dot":
            CS = plt.contourf(timestep, x, tau_dot, levels=200, cmap=cmc.lipari)
            CB = plt.colorbar(CS, orientation="horizontal", pad=0.2)
            CB.ax.set_title("tau_dot")
        elif val == "slip":
            CS = plt.contourf(timestep, x, slip, levels=20, cmap=cmc.lipari)
            CB = plt.colorbar(CS, orientation="horizontal", pad=0.2)
            CB.ax.set_title("cumulative slip (m)")

        plt.ylabel("position (m)")
        plt.xlabel("Timestep")

    elif orientation == "vertical":
        ddd -= ddd.min()
        
        if val == "v":
            CS = plt.contourf(ddd, timestep, levels=200, cmap=vcmap, norm=vnorm)
            CB = plt.colorbar(plt.cm.ScalarMappable(cmap=vcmap, norm=vnorm), orientation="horizontal", extend='both', ticks=ticks, pad=0.2)
            CB.ax.set_title("log v [m/s]")
        elif val=="dsigma":
            CS = plt.contourf(ddd, timestep, dsigma/1e6, levels=200, cmap=vcmap, norm=vnorm)
            CB = plt.colorbar(CS, orientation="horizontal", pad=0.2)
            CB.ax.set_title("dsigma [MPa]")
        elif val=="sigma":
            CS = plt.contourf(ddd, timestep, sigma/1e6, levels=200, cmap=cmc.batlow)
            CB = plt.colorbar(CS, orientation="horizontal", pad=0.2)
            CB.ax.set_title("sigma [MPa]")        
        elif val=="tau":
            CS = plt.contourf(ddd, timestep, tau/1e6, levels=200, cmap=cmc.batlow)
            CB = plt.colorbar(CS, orientation="horizontal", pad=0.2)
            CB.ax.set_title("tau [MPa]")
        elif val == "tau_sigma":
            CS = plt.contourf(ddd, timestep, tau_sigma, levels=200, cmap=cmc.lipari)
            CB = plt.colorbar(CS, orientation="horizontal", pad=0.2)
            CB.ax.set_title("tau/sigma")
        elif val=="theta":
            CS = plt.contourf(ddd, timestep, theta, levels=200, cmap=cmc.lipari)
            CB = plt.colorbar(CS, orientation="horizontal", pad=0.2)
            CB.ax.set_title("theta")
        elif val == "tau_dot":
            CS = plt.contourf(ddd, timestep, tau_dot, levels=200, cmap=cmc.lipari)
            CB = plt.colorbar(CS, orientation="horizontal", pad=0.2)
            CB.ax.set_title("tau_dot")
        elif val == "slip":
            CS = plt.contourf(ddd, timestep, slip, levels=50, cmap=cmc.lipari)
            CB = plt.colorbar(CS, orientation="horizontal", pad=0.2)
            CB.ax.set_title("cumulative slip (m)")
        
        plt.xlabel("downdip distance (m)")
        plt.ylabel("Timestep")
        plt.gca().invert_yaxis()

    else:
        print("Keyword 'orientation=%s' not recognised" % orientation)
        plt.close()

    plt.tight_layout()
    plt.show()
        
    return fig





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

def plot_vmax_fault(fault):
    """
    Plot time series of vmax for each fault
    """
    
    fig,ax= plt.subplots(ncols=1, nrows=2, figsize=figsize, squeeze=False)
    
    labels_fault = [str(i) for i in np.arange(1,len(fault)+1)]

    for i in np.arange(0,len(p.fault)):
        ax[i,0].plot(fault[i]["t"]/t_year, fault[i]["vmax_fault"])
        ax[i,0].set_xlabel("time (yr)")
        ax[i,0].set_ylabel("v (m/s)")
        ax[i,0].set_title("Fault "+labels_fault[i])
    
    fig.tight_layout()

    return fig