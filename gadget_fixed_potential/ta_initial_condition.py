"""
Generates initial conditions for a spike-like distribution of light PBHs in a
fixed central potential, with n(r) ∝ r^{-9/4}.

This script creates an HDF5 file with initial particle positions and velocities
for Gadget-4 simulations. Particles are distributed with a power-law density
profile characteristic of gravitational capture in a spike around a central
black hole. The distribution is sampled uniformly in the transformed coordinate
x = r^(3/4) to maintain the desired density profile.

Usage
-----
    python ta_initial_condition.py [mass_ratio] [num_particles]
    
Arguments
---------
    mass_ratio : float, optional
        Ratio of light PBH mass to central heavy PBH mass (default: 0.001).
    num_particles : int, optional
        Number of particles to generate (default: 1000).

Output
------
    IC_newton.hdf5 : Gadget-4 initial condition file
    t_collapse.txt : Collapse time in code units and seconds
"""

#-------
# ---
import numpy as np
import h5py
from matplotlib import pyplot as plt
import sys
import os

sys.path.append("../..")
import units as u

FloatType = np.float64
IntType = np.int32

# Section: Parameters
# ---
filename = "IC_newton.hdf5"

# Command-line arguments: mr & N of particles
# python ta_initial.py mr N
if len(sys.argv) > 1:
    mass_ratio = float(sys.argv[1])
else:
    mass_ratio = 0.001

if len(sys.argv) > 2:
    number_particles = int(sys.argv[2])
else:
    number_particles = 1000

if number_particles < 2:
    raise ValueError("number_particles must be >= 2")

# Section: Physical setup
# ---
M_heavy = 1.0 * u.Msun # Heavy PBH of solar mass
m = mass_ratio * M_heavy # N light PBHs; mr is a command-line argument

r_s = 2 * u.G * M_heavy / u.c**2

L = 1e-2 * u.pc
r_min = 300.0 * r_s
r_max = L

if r_min >= r_max:
    raise ValueError("r_min must be smaller than r_max")

# Section: Allocate arrays
# ---
Pos = np.zeros((number_particles, 3), dtype=FloatType)
Vel = np.zeros((number_particles, 3), dtype=FloatType)
Mass = np.zeros((number_particles, 1), dtype=FloatType)
ids = np.arange(number_particles, dtype=IntType)

Mass[:] = m

# Section: Shell-based sampling enforcing n(r) ∝ r^{-9/4}
# ---
np.random.seed(42)

# Transform variable x = r^{3/4}
x_min = r_min**(3.0/4.0)
x_max = r_max**(3.0/4.0)

Nshell = 100

x_edges = np.linspace(x_min, x_max, Nshell + 1)

N_base = number_particles // Nshell
remainder = number_particles % Nshell

idx = 0

for i in range(Nshell):
    Ni = N_base + (1 if i < remainder else 0)

    x_lo = x_edges[i]
    x_hi = x_edges[i+1]

    # sample uniformly in x
    U = np.random.rand(Ni)
    x = x_lo + U * (x_hi - x_lo)
    r = x**(4.0/3.0)

    # isotropic angles
    mu = 2.0 * np.random.rand(Ni) - 1.0
    phi = 2.0 * np.pi * np.random.rand(Ni)

    sin_theta = np.sqrt(1.0 - mu**2)

    Pos[idx:idx+Ni, 0] = r * sin_theta * np.cos(phi)
    Pos[idx:idx+Ni, 1] = r * sin_theta * np.sin(phi)
    Pos[idx:idx+Ni, 2] = r * mu

    idx += Ni

# Safety check: no nans 
if not np.isfinite(Pos).all():
    raise ValueError("Non-finite positions generated")

# Section: Time of collapse
# ---
radius = np.sqrt(np.sum(Pos**2, axis=1))
r_max_actual = radius.max()

V = (4*np.pi/3) * r_max_actual**3
rho_uniform = (m * number_particles) / V

Tcoll = float(
    np.sqrt(3*np.pi / (32 * u.G * rho_uniform)) / u.Tcode
)

# Section: Write HDF5 IC file
# ---
IC = h5py.File(filename, "w")

header = IC.create_group("Header")
part0 = IC.create_group("PartType0")
part1 = IC.create_group("PartType1")

NumPart = np.array([0, number_particles], dtype=IntType)

header.attrs.create("NumPart_ThisFile", NumPart)
header.attrs.create("NumPart_Total", NumPart)
header.attrs.create("NumPart_Total_HighWord", np.zeros(2, dtype=IntType))
header.attrs.create("MassTable", np.zeros(2, dtype=IntType))
header.attrs.create("Time", 0.0)
header.attrs.create("Redshift", 0.0)
header.attrs.create("BoxSize", 0)
header.attrs.create("NumFilesPerSnapshot", 1)

header.attrs.create("Omega0", 0.0)
header.attrs.create("OmegaB", 0.0)
header.attrs.create("OmegaLambda", 0.0)
header.attrs.create("HubbleParam", 1.0)

header.attrs.create(
    "Flag_DoublePrecision",
    1 if Pos.dtype == np.float64 else 0
)

part1.create_dataset("ParticleIDs", data=ids)
part1.create_dataset("Coordinates", data=Pos / u.Lcode)
part1.create_dataset("Masses", data=Mass / u.Mcode)
part1.create_dataset("Velocities", data=Vel / u.Vcode)

IC.close()

# Section: Save collapse time
# ---
t_collapse_code = float(Tcoll)
t_collapse_seconds = float(Tcoll * u.Tcode)

with open("t_collapse.txt", "w") as f:
    f.write(f"{t_collapse_code:.12e}\n")
    f.write(f"{t_collapse_seconds:.12e}\n")
