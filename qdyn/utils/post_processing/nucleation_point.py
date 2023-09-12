import numpy as np
import pandas as pd
import os
import sys
import pickle

def nucleation_point(mesh_dict, events_dict, save_output=True):
    
    """
    Calculate the location of nucleation points and their frequency for each fault.

    This function calculates the nucleation point of each event as the location 
    where the highest slip rate occurs in the first timestep of each event. 
    It also computes the frequency of nucleation points for each fault.

    Args
    ----
    - mesh_dict (dict)
    - events_dict (dict): A dictionary containing event information, including event dataframes (from run compute_events).
    - save_output (bool, optional): If True, save the output to a binary pickle file (default is True).

    Returns
    -------
    - dict: A dictionary containing nucleation point information and point frequency for each fault.

    Example
    -------
    - A complete example can be found in the Notebook compute_events.ipynb
    ```
    result = nucleation_point(p.mesh_dict, events_dict, save_output=True)
    ```

    Note:
    - The function stores nucleation point information in the "np" key of the output dictionary.
    - The function also calculates and stores point frequency in the "count_np" key of the output dictionary.
    - The function is only compatible with src that outputs ivmax_fault and vmax_fault

    Author: Constanza Rodriguez Piceda
    Date: 12.09.2023
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
        df_np_f = pd.DataFrame(columns=ot_ev_f.columns)
    
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
            df_vmax_f = df_event_f.iloc[i_tmin_f][df_event_f.iloc[i_tmin_f]["vmax_fault"] == df_event_f.iloc[i_tmin_f]["vmax_fault"].max()].copy()
    
            df_np_f = pd.concat([df_np_f, df_vmax_f])
    
        # reset indices to real event numbers
        df_np_f = df_np_f.reset_index(drop=True)
        df_np_f.index = range(1, len(events_f) + 1)
        df_np_f = df_np_f.rename_axis("n_event")
    
        ## Coordinates of nucleation points 
        # (uses the ivmax from the fault output to search for the coordinates in the mesh_dict array. We have to substract 1 to the index because python starts counting at 0, while fortran starts counting at 1.
        if len(df_np_f)!=0:
            # if the catalog is not empty
            df_np_f["x"] = mesh_dict["X"][(df_np_f["ivmax_fault"].values - 1).astype(int)]
            df_np_f["y"] = mesh_dict["Y"][(df_np_f["ivmax_fault"].values - 1).astype(int)]
            df_np_f["z"] = mesh_dict["Z"][(df_np_f["ivmax_fault"].values - 1).astype(int)]
        else:
            # if catalog is empty
            df_np_f["x"] = np.nan
            df_np_f["y"] = np.nan
            df_np_f["z"] = np.nan
        
        
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
        if len(df_np_f)!=0:
            df_count_np_f = df_np_f.value_counts(["x", "y", "z"]).to_frame()
            df_count_np_f = df_count_np_f.reset_index() # convert indices of coordinates x, y and z to columns
            df_count_np_f.rename({0:"count_np"}, axis=1, inplace=True) # rename column with frequency of np
        else:
            df_count_np_f = pd.DataFrame(columns = ["x", "y", "z", "count_np"])
            df_count_np_f["x"] = np.nan
            df_count_np_f["y"] = np.nan
            df_count_np_f["z"] = np.nan
            df_count_np_f["count_np"] = np.nan
    
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
        print(df_count_np_f)
        coords["count_np"] = pd.Series(dtype=float)
        coords.loc[coords["fault_label"] == fault_label, "count_np"] = df_count_np_f["count_np"]
        
        # DataFrame with coordinates of one individual fault
        coords_f = coords[coords["fault_label"]==fault_label].copy().reset_index(drop=True)
    
        # save DataFrame in dictionary
        model_dict["count_np"][fault_label] = coords_f
        

        if save_output:
            output_dict = model_dict
            # create a binary pickle file (output file)
            f = open("count_np.pkl","wb")
            # write the python object (dict) to pickle file
            pickle.dump(output_dict,f)
            
             # close file
            f.close()
    
    return model_dict