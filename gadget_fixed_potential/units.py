"""
CGS + Gadget-4 code units for fixed-potential simulations of light PBHs.

Usage
-----
    python units.py
"""

cm = 1.0
g = 1.0
s = 1.0

h = 0.678

Lcode = 3.085678e24 * cm / h  # length unit [cm] = (Mpc/h) in cm
Mcode = 1.989e43 * g / h      # mass unit [g] = (1e10 Msun/h) in g
Vcode = 1e5 * cm / s          # velocity [cm/s] = (km/s) in cm/s
Tcode = Lcode / Vcode         # time unit [s]

m = 1e2 * cm
km = 1e3 * m

pc = 3.0857e16 * m
kpc = 1e3 * pc
Mpc = 1e6 * pc

au = 1.495978707e11 * m

kg = 1e3 * g
Msun = 1.989e33 * g

yr = 60 * 60 * 24 * 365.25 * s
kyr = 1e3 * yr
Myr = 1e6 * yr

G = 6.67430e-11 * m ** 3 / kg / s ** 2
c = 299_792_458 * m / s