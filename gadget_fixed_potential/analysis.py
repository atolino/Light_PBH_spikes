
"""
Basic analysis tools for the fixed-potential gadget snapshots.

Usage
-----
    python analysis.py
"""

import sys
import numpy as np
import h5py
import csv
import os

import units as u

FloatType = np.float64

def load_snapshot(fname, sort=True, cgs_units=True):
    """
    Load a Gadget-4 HDF5 snapshot and extract particle data.
 
    """
    try:
        data = h5py.File(fname, "r")
    except:
        print("Could not open file: " + fname + " !")
        sys.exit(1)

    time = FloatType(data["Header"].attrs["Time"])
    Pos = np.array(data["PartType1"]["Coordinates"], dtype=FloatType)
    Vel = np.array(data["PartType1"]["Velocities"], dtype=FloatType)
    Mass = np.array(data["PartType1"]["Masses"], dtype=FloatType)
    ParticleIDs = np.array(data["PartType1"]["ParticleIDs"])

    # If this flag is set, convert from code units to cgs units
    if cgs_units:
        time *= u.Tcode
        Pos *= u.Lcode
        Vel *= u.Vcode
        Mass *= u.Mcode

    # The order of the particles in the snapshots is not conserved,
    # so let's sort the particles so that the i-th entry in the
    # array always corresponds to the same particle
    if sort is True:
        sortargs = np.argsort(ParticleIDs)
    else:
        sortargs = np.arange(len(ParticleIDs))
    return time, Pos[sortargs, :], Vel[sortargs, :], Mass[sortargs], ParticleIDs

def load_all_snapshots(fname_root, i_list, cgs_units=True):
    """
    Load multiple snapshots and stack them into time-series arrays.
  
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