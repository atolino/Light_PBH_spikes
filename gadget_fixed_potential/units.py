"""
CGS + Gadget-4 code units for fixed-potential simulations of light PBHs.

This module defines unit conversions between CGS physical units and internal
Gadget-4 code units. All unit definitions follow Gadget-4 v4 standards.

Units System
============
Gadget-4 internal units are set to:
    - Length: 1 Mpc/h
    - Mass: 1e10 Msun/h  
    - Velocity: 1 km/s
    - Time: derived from length/velocity

To convert a quantity from physical units to code units, divide by the corresponding
code unit constant (e.g., position_code = position_physical / Lcode).

To convert back to physical units, multiply by the code unit constant
(e.g., position_physical = position_code * Lcode).

Physical Constants
==================
- Gravitational constant G: 6.67430e-11 m^3/(kg s^2)
- Speed of light c: 299,792,458 m/s
"""

#-------

# Section: Base CGS units
# ---
cm = 1.0
g = 1.0
s = 1.0

# Section: Gadget code units
# ---
h = 0.678

Lcode = 3.085678e24 * cm / h   # length unit [cm]  = (Mpc/h) in cm
Mcode = 1.989e43 * g / h       # mass unit   [g]   = (1e10 Msun/h) in g
Vcode = 1e5 * cm / s           # velocity    [cm/s]= (km/s) in cm/s
Tcode = Lcode / Vcode          # time unit    [s]

# Section: Convenience units (CGS)
# ---
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

# Section: Physical constants (SI expressed via CGS primitives above)
# ---
G = 6.67430e-11 * m**3 / kg / s**2
c = 299_792_458 * m / s