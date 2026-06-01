
"""
Basic analysis tools for the fixed-potential gadget snapshots.

This module provides utilities for loading and processing Gadget-4 simulation
snapshots in HDF5 format. It handles coordinate unit conversions and particle
sorting for consistent time-series analysis.
"""

#-------

# Section: Imports
# ---
import sys
import numpy as np
import h5py
import csv
import os

import units as u

FloatType = np.float64

# Section: Functions
# ---

def load_snapshot(fname, sort=True, cgs_units=True):
    """
    Load a Gadget-4 HDF5 snapshot and extract particle data.
    
    Parameters
    ----------
    fname : str
        Path to the HDF5 snapshot file.
    sort : bool, optional
        If True, sort particles by ID for consistent ordering (default: True).
    cgs_units : bool, optional
        If True, convert from code units to CGS units (default: True).
    
    Returns
    -------
    time : float
        Simulation time at snapshot.
    Pos : ndarray, shape (N, 3)
        Particle positions [cm if cgs_units=True, code units otherwise].
    Vel : ndarray, shape (N, 3)
        Particle velocities [cm/s if cgs_units=True, code units otherwise].
    Mass : ndarray, shape (N,)
        Particle masses [g if cgs_units=True, code units otherwise].
    ParticleIDs : ndarray, shape (N,)
        Particle identifiers (not sorted).
    """
    try:
        data = h5py.File(fname, "r")
    except:
        print("Could not open file: " + fname + " !")
        sys.exit(1)

    time = FloatType(data["Header"].attrs["Time"])
    Pos = np.array(data["PartType1"]["Coordinates"], dtype = FloatType) 
    Vel = np.array(data["PartType1"]["Velocities"], dtype = FloatType)
    Mass = np.array(data["PartType1"]["Masses"], dtype = FloatType)
    ParticleIDs = np.array(data["PartType1"]["ParticleIDs"])
    
    # If this flag is set, convert from code units to cgs units
    if (cgs_units):
        time *= u.Tcode
        Pos  *= u.Lcode
        Vel  *= u.Vcode
        Mass *= u.Mcode
    
    # The order of the particles in the snapshots is not conserved,
    # so let's sort the particles so that the i-th entry in the
    # array always corresponds to the same particle
    if sort == True:
        sortargs = np.argsort(ParticleIDs)
    else:
        sortargs = np.arange(len(ParticleIDs))
    return time, Pos[sortargs,:], Vel[sortargs,:], Mass[sortargs], ParticleIDs

# Load a list of snapshots (labelled by i_list) and concatenate
# them into a big array of positions with dimensions (N_snapshots, N_particles, 3)
# Also outputs the list of times of the snapshots
def load_all_snapshots(fname_root, i_list, cgs_units = True):
    """
    Load multiple snapshots and stack them into time-series arrays.
    
    Parameters
    ----------
    fname_root : str
        Root directory containing snapshot files named 'snapshot_XXX.hdf5'.
    i_list : array-like
        List of snapshot indices to load (e.g., [0, 1, 2, ...]).
    cgs_units : bool, optional
        If True, convert to CGS units (default: True).
    
    Returns
    -------
    ts : ndarray, shape (N_snapshots,)
        Simulation times at each snapshot.
    Pos_full : ndarray, shape (N_snapshots, N_particles, 3)
        Stacked particle positions for all snapshots.
    Vel_full : ndarray, shape (N_snapshots, N_particles, 3)
        Stacked particle velocities for all snapshots.
    """
    for i, i_file in enumerate(i_list):
        filename = os.path.join(fname_root, 'snapshot_%03d.hdf5' % i_file)
        t, Pos, Vel, _, _ = load_snapshot(filename, sort=True, cgs_units=cgs_units)
        if i == 0:
            ts = []
            N_part = len(Pos)
            Pos_full = np.zeros((len(i_list), N_part, 3))
            Vel_full = np.zeros((len(i_list), N_part, 3))

        ts.append(t)
        Pos_full[i,:,:] = Pos
        Vel_full[i,:,:] = Vel
    ts = np.array(ts)
    return ts, Pos_full, Vel_full