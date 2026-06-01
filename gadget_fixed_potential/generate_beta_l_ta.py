"""
Compute beta and angular-momentum diagnostics from Gadget-4 outputs.

This script supports two related tasks:

1) Simulation-derived angular momentum and torque time series:
   - Computes L = r x v from snapshots
   - Computes tau = dL/dt
   - Extracts the first turning point ("closest approach") per particle
   - Saves:
        l_sim_min_<mr>.npy
        r_closest_<mr>.npy
        r_in_<mr>.npy
        torque_<mr>.npy
        torque_vec_<mr>.npy

2) Initial-condition beta parameter:
   - Computes beta = tau_i / tau_c at t = t0 for each particle
   - Saves:
        beta_in_<mr>_new.npy

Usage
-----
    python generate_beta_l_ta.py

The script processes snapshot data from directories named output_<mr>_ta_001,
where <mr> is formatted in scientific notation (e.g., 1e-05 for mr=1e-5).
"""

#-------

# Section: Imports
# ---
import glob
import os
import re
from typing import Callable, Tuple

import numpy as np

import units as u
from analysis import load_all_snapshots


# Section: load helper
# ---
def load_or_compute(filename: str, compute_func: Callable[[], np.ndarray]) -> np.ndarray:
    """
    Load data from file or compute and save it interactively.
    
    Checks if a file exists. If yes, prompts user to either load or recompute.
    If no, computes the data and saves it automatically.
    
    Parameters
    ----------
    filename : str
        Path to the .npy file to load or save.
    compute_func : callable
        Function that computes the data if file doesn't exist or user chooses to recompute.
        Should return an ndarray.
    
    Returns
    -------
    ndarray
        Either loaded or newly computed data.
    """
    if os.path.exists(filename):
        choice = input(f"{filename} already exists. Replace it? (y/n): ").strip().lower()
        if choice == "y":
            print(f"Recomputing {filename} ...")
            data = compute_func()
            np.save(filename, data)
        else:
            print(f"Loading {filename} ...")
            data = np.load(filename)
    else:
        print(f"Computing {filename} ...")
        data = compute_func()
        np.save(filename, data)
    return data


# Section: load snapshots
# ---
def load_data(mr: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Load and organize all snapshots for a given mass ratio.
    
    Scans the output directory for all snapshot files matching the mass ratio tag,
    loads them in order, and returns stacked time-series data.
    
    Parameters
    ----------
    mr : float
        Mass ratio (light PBH mass / heavy PBH mass).
    
    Returns
    -------
    t : ndarray, shape (N_snapshots,)
        Simulation time at each snapshot [CGS units].
    Pos : ndarray, shape (N_snapshots, N_particles, 3)
        Particle positions at each snapshot [CGS units].
    Vel : ndarray, shape (N_snapshots, N_particles, 3)
        Particle velocities at each snapshot [CGS units].
    r : ndarray, shape (N_snapshots, N_particles)
        Particle distances from origin.
    Npart : int
        Number of particles.
    """
    """Load positions/velocities from all snapshots for a given mass ratio tag."""
    mr_tag = f"{mr:.0e}"
    outputs_dir = f"./output_{mr_tag}_ta_001"
    snapshot_files = glob.glob(f"{outputs_dir}/snapshot_*.hdf5")

    snap_nums = []
    for f in snapshot_files:
        m = re.search(r"snapshot_(\d+)\.hdf5$", f)
        if m:
            snap_nums.append(int(m.group(1)))
    snap_nums = np.sort(snap_nums)

    t, Pos, Vel = load_all_snapshots(outputs_dir, snap_nums, cgs_units=True)

    Nsnap, Npart, _ = Pos.shape
    print(f"Loaded {Nsnap} snapshots, {Npart} particles")

    r = np.linalg.norm(Pos, axis=2)
    return t, Pos, Vel, r, Npart


# Section: Simulation analysis (L, tau, closest approach)
# ---
def generate_lsim_data(mr: float):
    """
    Compute angular momentum and torque diagnostics from snapshots.
    
    Calculates time-dependent angular momentum L = r x v and torque tau = dL/dt
    for each particle across all snapshots. Identifies the closest approach point
    (first turning point) for each particle's trajectory. Saves all derived quantities.
    
    Parameters
    ----------
    mr : float
        Mass ratio (used for data file naming).
    
    Returns
    -------
    l_min : ndarray, shape (N_particles,)
        Angular momentum at closest approach for each particle.
    r_init : ndarray, shape (N_particles,)
        Initial distance from origin for each particle.
    r_closest : ndarray, shape (N_particles,)
        Distance at closest approach for each particle.
    tau_mod : ndarray, shape (N_snapshots, N_particles)
        Torque magnitude time series for all particles.
    
    Saves
    -----
    l_sim_min_<mr>.npy, r_closest_<mr>.npy, r_in_<mr>.npy, torque_<mr>.npy, torque_vec_<mr>.npy
    """
    t, Pos, Vel, r, Npart = load_data(mr)
    r_init = r[0]

    # Angular momentum: L = r x v
    l_vec = np.cross(Pos, Vel)                # (Nt, Npart, 3)
    l_mod = np.linalg.norm(l_vec, axis=2)     # (Nt, Npart)

    # Torque: tau = dL/dt
    tau_vec = np.gradient(l_vec, t, axis=0)   # (Nt, Npart, 3)
    tau_mod = np.linalg.norm(tau_vec, axis=2) # (Nt, Npart)

    # Closest approach: first turning point where r starts increasing again
    def find_min_approach(r_i, l_i):
        for j in range(len(r_i) - 1):
            if r_i[j + 1] > r_i[j]:
                return l_i[j], r_i[j]
        return l_i[-1], r_i[-1]

    l_min = np.zeros(Npart)
    r_closest = np.zeros(Npart)
    for p in range(Npart):
        l_min[p], r_closest[p] = find_min_approach(r[:, p], l_mod[:, p])

    mr_tag = f"{mr:.0e}"
    np.save(f"l_sim_min_{mr_tag}.npy", l_min)
    np.save(f"r_closest_{mr_tag}.npy", r_closest)
    np.save(f"r_in_{mr_tag}.npy", r_init)
    np.save(f"torque_{mr_tag}.npy", tau_mod)
    np.save(f"torque_vec_{mr_tag}.npy", tau_vec)

    print("Saved l_min, r_closest, r_in, torque")
    return l_min, r_init, r_closest, tau_mod


# Section: Beta computation at initial conditions
# ---
def compute_beta(mr: float) -> np.ndarray:
    """
    Compute beta = tau_i / tau_c for each particle at initial conditions.
    
    Evaluates the dimensionless torque parameter beta at t=0 by computing:
    - tau_i: actual torque from particle-particle interactions in the ICs
    - tau_c: characteristic torque from analytical scaling relations
    
    Uses the spike density profile n(r) ∝ r^{-9/4} for the characteristic torque.
    
    Parameters
    ----------
    mr : float
        Mass ratio (light PBH mass / heavy PBH mass).
    
    Returns
    -------
    beta_p : ndarray, shape (N_particles,)
        Dimensionless torque parameter beta for each particle.
    
    Saves
    -----
    beta_in_<mr>_new.npy : Computed beta values (if not already present).
    """
    _, Pos, _, r, Npart = load_data(mr)

    Pos0 = Pos[0]   # initial positions (CGS units)
    mH = u.Msun
    mL = mr * mH

    # Spike number-density normalization n(r) ∝ r^{-9/4}
    L = 1e-2 * u.pc
    p_exp = 9 / 4
    Vtilde = 4 * np.pi * (L ** (3 - p_exp)) / (3 - p_exp)
    ntilde = Npart / Vtilde

    def beta_computation() -> np.ndarray:
        beta_p = np.zeros(Npart)

        for p in range(Npart):
            xp = Pos0[p]
            dr = Pos0 - xp
            rnorm = np.linalg.norm(dr, axis=1)

            # Exclude self-interaction term (r = 0)
            mask = rnorm > 0
            dr = dr[mask]
            rnorm = rnorm[mask]

            # Torque contribution from each other particle
            cross_terms = np.cross(xp, dr)
            tau_vecs = -u.G * mL**2 * cross_terms / rnorm[:, None] ** 3
            torque_p = tau_vecs.sum(axis=0)
            tau_abs = np.linalg.norm(torque_p)

            # Characteristic torque (analytical)
            xpnorm = np.linalg.norm(xp)
            tau_c = 2.6 * u.G * mL**2 * (ntilde ** (2 / 3)) * xpnorm ** (-0.5)

            beta_p[p] = tau_abs / tau_c

        return beta_p

    mr_tag = f"{mr:.0e}"
    return load_or_compute(f"beta_in_{mr_tag}_new.npy", beta_computation)


# Section: CLI entry point
# ---
def main() -> None:
    mr = 1e-5
    compute_beta(mr)


if __name__ == "__main__":
    main()