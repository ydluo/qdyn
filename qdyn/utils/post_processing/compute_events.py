import numpy as np
import pandas as pd
import os
import sys
import pickle
def compute_events(set_dict, mesh_dict, ot_fault, vmax=0.01,tmin=0,tmax=None, save_output=True):
        
    """
    Script to compute seismic events

    Author: Constanza Rodriguez Piceda
    Date: 12.09.2023

    Required libraries: qdyn, numpy, os, sys, pandas and pickle


    Example Usage
    -------------
    from extract_events import compute_events
    $ compute_events(p.set_dict, p.mesh_dict, p.fault, vmax=0.01,tmin=0,tmax=None, save_output=True)


    Arguments
    ---------
    vmax = velocity threshold (m/s) (default = 0.01 ms
    tmin = minimum cut-off simulation time to start considering events (yr) (default = 0yr)
    tmax = maximum cut-off simulation time to stop considering events (yr) (default = tmax simulation)
    save_output = boolean to save output dictionary as a binary "events.pkl" that can be open with python (default=True)

    Preparation
    -----------
    - The following variables should be created prior to calling the function (see example Notebook on how to create them)
        - set dict and mesh_dict: "p.set_dict" and "p.mesh_dict"
        - outputs ot_fault: "p.fault"
    - vmax: velocity threshold to define seismic events
    - t_min and t_max: min and max time of simulation to consider seismic events
    - save_output: save the output dictionary as a pkl binary file

    Returns
    -------
    - Dictionary "model_dict" storing the following DataFrames:
        - model_dict["ev"][<fault_label>]:  store event information 
                            (index (n event), min_t, max_t, cum_slip, cum_potency, peak_v, dt_event, t_event, t_interevent_intrafault, Mw, fault_label)
        - model_dict["vmax_ev"][<fault_label>]: store time series with vmax of events (index (n simulation time-step),
                                                step, t, potcy_fault, pot_rate_fault, vmax_fault, ivmax_fault, deltat, n_event, cum_slip_fault)
        - model_dict["seq"]: store information of sequences of events in a Dictionary. A "sequence" refers to a collection of events within a
        specific fault that occur during a time period with no events in another fault of the network.
            (key (seqin_seqfin), index(seq), min_t         max_t  t_interevent_interfault)

    Note
    -----
      - This function is only compatible with version that includes vmax_fault and ivmax_fault
      - For now src code is only outputing vmax_fault and ivmax_fault for each fault (not tau, sigma, theta as for ot_vmax). Thus the output dict doesn't contain those
      additional quantities.

    """


    """
    INITIAL PARAMETERS
    """

    # Predefine parameters
    t_yr = 3600 * 24 * 365.0    # Seconds per year

    # # Instantiate the QDYN class object
    # p = qdyn()


    """
    OPEN OUTPUTS
    """

    # Create empty dictionary where DataFrames of outputs will be stored
    model_dict={}

    #  assign dictionary with settings of meshdict to dictionary storing models
    model_dict["mesh_dict"] = mesh_dict

    #  assign dictionary with settings of meshdict to dictionary storing models
    model_dict["set_dict"] = set_dict

    """
    COMPUTE EVENT QUANTITIES
    To compute the quantities per events we will work with the output output_fault. 
    First we have to extract the time-steps with seismic events for each fault. 
    To do so, we use the output file output_fault to identify the time-steps 
    where the element with highest slip rate exceeds a given threshold, in this case 0.01m/s.
    This way we obtain the initial and final times of each event.
    Then we use the output_fault to calculate the slip and potency per event
    the interevent time within each fault, and the Mw.
    We also extract the peak velocity and its corresponding time and the time interval of the event
    Finally, if working with multiple faults we calculate the interevent time between faults
    (called t_interevent_interfault)
    """

    # empty dictionary to store events    
    model_dict["ev"] = {} # store information events
    model_dict["seq"] = {} # store information of sequences of events
    model_dict["vmax_ev"] = {} # store information of vmax of events

    # fault labels
    fault_labels = np.unique(model_dict["mesh_dict"]["FAULT_LABEL"])

    #########################
    #### FOR EACH FAULT #####
    #########################
    ### Events have to be calculated separately for each fault ###

    for i_f, fault_label in enumerate(fault_labels):

        # filter timesteps with velocity and time threshold in DataFrame output ot_fault
        if tmax!= None:
            ot_ev_f = ot_fault[i_f][(ot_fault[i_f]["vmax_fault"] >= vmax) & (ot_fault[i_f]["t"] >= tmin*t_yr) & (ot_fault[i_f]["t"] <= tmax*t_yr)].copy()
        else:
            ot_ev_f = ot_fault[i_f][(ot_fault[i_f]["vmax_fault"] >= vmax) & (ot_fault[i_f]["t"] >= tmin*t_yr)].copy()
        
        # Number of elements individual Fault
        n_f = len(model_dict["mesh_dict"]["FAULT_LABEL"][model_dict["mesh_dict"]["FAULT_LABEL"]==fault_label])

        # Calculate the time difference between previous and following time-step
        ot_ev_f["deltat"] = ot_ev_f['t'].diff()

        # calculate the number of event for each row (one event is considered independent if it occurs >1s after the previous time-step)
        ot_ev_f["n_event"] = (ot_ev_f["deltat"] >=1).cumsum() + 1

        # Group rows of vmax DataFrame by event number
        vmax_ev_f = ot_ev_f.groupby("n_event")

        # approach 1: Filter time-steps of fault outputs where there are seismic events

        ### calculate CUMULATIVE SLIP per fault per timestep ###

        ## 1st calculate fault area
        # get fault labels
        fault_lbl = model_dict["mesh_dict"]["FAULT_LABEL"]

        # Create a boolean mask to identify the positions where fault_labels == fault_label
        mask = fault_lbl == fault_label

        # Use the mask to select the corresponding Z values
        z_values = model_dict["mesh_dict"]["Z"][mask]

        # calculate fault distance along z
        Dz = max(z_values) - min(z_values)

        # extract dip
        dip = model_dict["mesh_dict"]["DIP_W"][0]

        # calculate width along dip (without half lengths first)
        Dw = Dz / np.cos(np.radians(dip))

        # Access the "DW" values from the mesh_dict for the specified model_name to get the half length of top and bottom elements
        dw_values = model_dict["mesh_dict"]["DW"][mask]

        # calculate half length of top and bottom elements
        half_dw_bottom  = dw_values[0]/2
        half_dw_top  = dw_values[-1]/2

        # calculate width along dip
        Dw = Dw + half_dw_bottom + half_dw_top

        # # Use the mask to select the corresponding X values
        x_values = model_dict["mesh_dict"]["X"][mask]

        # calculate fault distance along z
        Dx = max(x_values) - min(x_values)

        # array with unique x coordinates
        x_unique = np.unique(x_values)

        # calculate half length of right and left elements
        half_dx_right = x_unique[1] - x_unique[0]
        half_dx_left = x_unique[-1] - x_unique[-2]

        # calculate length along strike
        Dx = Dx + half_dx_right + half_dx_left

        # calculate area
        area1 = Dx*Dw

        ## calculate cumulative slip per fault (this is slip since t0)
        ot_ev_f["cum_slip_fault"] = ot_ev_f["potcy_fault"] / area1

        ### end calculation cumulative slip per fault ###

        # calculate the time difference between previous and following time-step
        ot_ev_f["deltat"] = ot_ev_f['t'].diff()

        # calculate the number of event for each row (one event is considered independent if it occurs >1s after the previous time-step)
        ot_ev_f["n_event"] = (ot_ev_f["deltat"] >=1).cumsum() + 1

        # group rows of fault DataFrame by number of event
        ev_f = ot_ev_f.groupby("n_event")

        # Calculate start and end time of each event
        df_ev_f = ev_f["t"].agg(["min", "max"]) # calculate start and end time of event and create DataFrame with values
        df_ev_f = df_ev_f.rename(columns={"min": "min_t", "max": "max_t"}) # rename columns
        df_ev_f = df_ev_f.reindex(columns=['min_t', 'max_t']) # reorder columns

        ## Calculate cumulative slip of each event
        #df_ev1["cum_slip"] = ev1["slip_dt_fault"].sum()/n_f1 # old wrapper
        # extract the index corresponding to the max and min time of each event
        idx_t_max_f = ev_f["t"].idxmax()
        idx_t_min_f = ev_f["t"].idxmin()

        # Calculate cumulative slip of each event (as the substraction between time-step at tmax and tmin)
        df_ev_f["cum_slip"] = (ot_ev_f.loc[idx_t_max_f, 'cum_slip_fault'].reset_index(drop=True) - ot_ev_f.loc[idx_t_min_f, 'cum_slip_fault'].reset_index(drop=True)).values


        # Calculate cumulative potency of each event
        # df_ev1["cum_potency"] = ev1["potcy_fault"].sum() # old wrapper
        df_ev_f["cum_potency"] = (ot_ev_f.loc[idx_t_max_f, 'potcy_fault'].reset_index(drop=True) - ot_ev_f.loc[idx_t_min_f, 'potcy_fault'].reset_index(drop=True)).values

        # Calculate peak velocity for each event
        df_ev_f["peak_v"] = ot_ev_f["vmax_fault"].max()

        # Calculate duration of each event
        df_ev_f["dt_event"] = df_ev_f["max_t"] - df_ev_f["min_t"]

        # Calculate time of event based on the time of occurence of the max peak velocity
        t_event_f = ot_ev_f["t"].loc[ot_ev_f.groupby("n_event").vmax_fault.idxmax()] # find the time-step where the peak velocity occurs for each event
        t_event_f.index = df_ev_f.index.to_list() # rename the indices so they match the number of event
        df_ev_f["t_event"] = t_event_f # add the column to the DataFrame storing the information about the events

        # Calculate interevent-time within fault (recurrence time)
        df_ev_f["t_interevent_intrafault"] = df_ev_f['max_t'] - df_ev_f['min_t'].shift(1)

        """
        MAGNITUDE
        """

        """
        Calculate seismic moment and moment magnitude for each event
        """
        G = model_dict["set_dict"]["MU"] # shear modulus

        # Calculate seismic moment Mo of each event in the fault as P*G
        Mo_f = df_ev_f["cum_potency"]*G
        Mo_f=Mo_f.rename("Mo") # rename Series

        # Calculate moment magnitude Mw as 2/3*np.log10(Mo1) - 6.06
        Mw_f = 2/3*np.log10(Mo_f) -6.06
        Mw_f=Mw_f.rename("Mw") # rename Series

        # Assign Mw to column in event DataFrame 
        df_ev_f["Mw"] = Mw_f

        ## end magnitude calculation ##

        # assign fault label to Dataframe with events
        df_ev_f["fault_label"] = fault_label

        # save events information in a DataFrame that can be accessed through the fault label
        model_dict["ev"][fault_label] = df_ev_f
        model_dict["vmax_ev"][fault_label] = ot_ev_f

    #########################################
    #### INTEREVENT TIME BETWEEN FAULTS #####
    #########################################
    # Get the list of keys in model_dict["ev"]
    keys = list(model_dict["ev"].keys())

    # Do not consider the simulation case with only 1 fault

    if (len(keys) != 1):    
        # Iterate over each pair of entries
        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):
                # Concatenate DataFrames of events
                df_ev = pd.concat([model_dict["ev"][keys[i]], model_dict["ev"][keys[j]]], axis=0)

                # Reorder dataframe according to min_t of events
                df_ev = df_ev.sort_values(by="min_t")

                # Define "sequences" that group consecutive events corresponding to one fault
                df_ev["seq"] = (df_ev["fault_label"] != df_ev["fault_label"].shift(1)).cumsum()
                seq = df_ev.groupby("seq") # group events by the sequences

                df_seq = seq["t_event"].agg(["min", "max"]) # calculate start and end time of each sequence and create DataFrame with values
                df_seq = df_seq.rename(columns={"min": "min_t", "max": "max_t"}) # rename columns
                df_seq = df_seq.reindex(columns=['min_t', 'max_t']) # reorder columns
                df_seq["t_interevent_interfault"] = df_seq['min_t'] - df_seq['min_t'].shift(1) # calculate interevent time between faults

                # Save the result in the model_dict["seq"] dictionary with appropriate key
                model_dict["seq"][f"{keys[i]}_{keys[j]}"] = df_seq

    ######################
    #### SAVE OUTPUT #####
    ######################
    # save relevant entries to binary pkl file
    # Get the relevant keys and values from p.model_dict
    if (len(keys) != 1): 
        output_keys = ["ev", "vmax_ev", "seq"]
    else:
        output_keys = ["ev", "vmax_ev"]

    output_dict = {key: model_dict[key] for key in output_keys}

    if save_output:
        # create a binary pickle file (output file)
            f = open("events.pkl","wb")
            # write the python object (dict) to pickle file
            pickle.dump(output_dict,f)
        
            # close file
            f.close()

    return model_dict