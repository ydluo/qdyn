import numpy as np
import pandas as pd
import os
import sys
import pickle

def nucleation_point(mesh_dict, events_dict):

    """
    Script to calculate location of nucleation point and its frequency for each fault

    Author: Constanza Rodriguez Piceda
    Date: 27.07.2023

    Required libraries: qdyn, numpy, os, sys, pandas and pickle


    Usage:

    For an example workflow on how to use this function check the Notebook compute_events.py

    Preparation:

    The following variables should be created prior to calling the function (see example Notebook on how to create them)
    - mesh_dict: "p.mesh_dict"
    - events_dict: output dictionary of compute_events function

    Returns:
    - Dictionary "model_dict" storing the additional following DataFrames:
        - model_dict["np"][<fault_label>]:  store nucleation point information  of events
                            (iindex (n event), step (timestep of nucleation), ivmax (index of element with peak velocity at nucleation), 
                            v (nucleation velocity), theta, tau, dtau_dt, slip (cumulative slip from start of simualtion), sigma, fault_label,
                            deltat (time duration of event), n_event, and coordinates x y z of nucleation point.)
        - model_dict["vmax_ev"][<fault_label>]: store information about absolute frequency of nucleation location:
                                                (coordinates x y z, fault_label and count_np (absolute frequency)). Coordinates where no event
                                                occurs have "count_np"= NaN

    Note: This function is only compatible with versions > 3.0.0

    """

    """
    Calculate nucleation point of each event as location where the highest slip rate occurs in the first timestep of each event
    """

    model_dict = events_dict

    # create empty dictionary to store nucleation points
    model_dict["np"] = {}

    # create an empty dictionary where point frequency and coordinates of nucleation points will be stored
    model_dict["count_np"] = {}

    # fault labels

    fault_labels = np.unique(mesh_dict["FAULT_LABEL"])

    #########################
    #### FOR EACH FAULT #####
    #########################
    for i_f, fault_label in enumerate(fault_labels):

        df_ev_f = model_dict["ev"][fault_label]
        ot_ev_f = model_dict["vmax_ev"][fault_label]

        # Create empty data frame where info of nucleation points will be stored
        df_np_f = pd.DataFrame()

        # number of events in Fault 1
        events_f = df_ev_f.index

        # Loop over number of events
        for n_event in events_f:

            #Snapshots of one event 
            df_event_f = ot_ev_f[ot_ev_f["n_event"]==n_event]
            df_event_f=df_event_f.reset_index(drop=True) # reset indices to start counting from 0.

            # indices of rows with tmin in event
            i_tmin_f = df_event_f.index[df_event_f["t"]==df_event_f["t"].min()]

            # row with v max in snapshot with tmin (1st snapshot of the event)
            df_vmax_f = df_event_f.iloc[i_tmin_f][df_event_f.iloc[i_tmin_f]["v"] == df_event_f.iloc[i_tmin_f]["v"].max()].copy()

            df_np_f = df_np_f.append(df_vmax_f)

        # reset indices to real event numbers
        df_np_f = df_np_f.reset_index(drop=True)
        df_np_f.index = range(1, len(events_f) + 1)
        df_np_f = df_np_f.rename_axis("n_event")

        ## Coordinates of nucleation points 
        # (uses the ivmax from the ot_vmax output to search for the coordinates in the mesh_dict array. We have to substract 1 to the index because python starts counting at 0, while fortran starts counting at 1.
        df_np_f["x"] = mesh_dict["X"][df_np_f["ivmax"].values - 1]
        df_np_f["y"] = mesh_dict["Y"][df_np_f["ivmax"].values - 1]
        df_np_f["z"] = mesh_dict["Z"][df_np_f["ivmax"].values - 1]
        
        
        ## assign Mw to DataFrame df_np_f
        df_np_f["Mw"] = model_dict["ev"][fault_label]["Mw"]
        
        # save nucleation points to dictionary
        model_dict["np"][fault_label] = df_np_f
        

        """
        Calculate nucleation point frequency
        """
        
        # series with element coordinates from mesh_dict
        x_coord = mesh_dict["X"]
        y_coord = mesh_dict["Y"]
        z_coord = mesh_dict["Z"]

        # Pandas data Frame with coordinates
        coords = pd.DataFrame(list(zip(x_coord, y_coord, z_coord)), columns = ["x", "y", "z"])
        coords["fault_label"] = mesh_dict["FAULT_LABEL"]
        
        # Count of nucleation points fault
        df_count_np_f = df_np_f.value_counts(["x", "y", "z"]).to_frame()[0]
        df_count_np_f = df_count_np_f.reset_index() # convert indices of coordinates x, y and z to columns
        df_count_np_f.rename({0:"count_np"}, axis=1, inplace=True) # rename column with frequency of np

        # row with frequency
        #coords["count_np"] = np.nan
        df_count_np_f["i_np"] = np.nan


        # Find indices of nucleation points of F1 in the coordinate mesh
        list_i_np_f = [] # empty list where nucleation points of F1 will be stored
        for row in np.arange(0,len(df_count_np_f)):
            i_np_f = coords.index[(coords['x'] == df_count_np_f["x"][row]) & (coords['y'] == df_count_np_f["y"][row]) & (coords['z'] == df_count_np_f["z"][row])].tolist()[0]
            list_i_np_f.append(i_np_f)
            df_count_np_f.loc[row, "i_np"] = int(i_np_f)

        # Set indices of nucleation points as index of the DataFrame for F1
        df_count_np_f.index = list_i_np_f
        df_count_np_f.index.name = "i_np"
        
        # create a column with count_np in the DataFrame with coordinates
        coords["count_np"] = pd.Series(dtype=float)
        coords.loc[coords["fault_label"] == fault_label, "count_np"] = df_count_np_f["count_np"]
        
        # DataFrame with coordinates of one individual fault
        coords_f = coords[coords["fault_label"]==fault_label].copy().reset_index(drop=True)

        # save DataFrame in dictionary
        model_dict["count_np"][fault_label] = coords_f

    return model_dict