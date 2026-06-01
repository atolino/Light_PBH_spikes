"""
This script computes the distribution of the dimensionless parameter epsilon = l/A,
where l is the simulated angular momentum of particles in the GADGET4 simulation and A is a
characteristic torque-related coefficient. The distribution of epsilon is then compared to
the analytical prediction g_numeric(epsilon) from the paper.
"""

#-------

# Section: Imports
# ---
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.integrate import quad
import units as u

# Section: Matplotlib style
# ---
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{mathpazo}'
plt.rcParams['axes.linewidth'] = 2
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Section: Analytical g_numeric
# ---
def g_numeric(gamma):    
    """
    Compute the analytical distribution g_numeric(gamma).
    
    This is the predicted dimensionless distribution from torque statistics theory.
    
    Parameters
    ----------
    gamma : float or ndarray
        Dimensionless parameter epsilon = l/A.
    
    Returns
    -------
    float or ndarray
        Probability density at gamma.
    """    
    I = 1 + (gamma**(9/4) / (np.sqrt(3)*3))
    return gamma * 0.595 * (I**(-9/4))

# Section: Compute epsilon = l/A and A values
# ---
def compute_epsilon(mr):
    """
    Compute epsilon = l/A from simulation data.
    
    Loads precomputed angular momentum and orbital parameters, calculates the
    characteristic torque coefficient A using integral and particle statistics,
    and computes the dimensionless ratio epsilon for each particle.
    
    Parameters
    ----------
    mr : float
        Mass ratio (used to name data files).
    
    Returns
    -------
    epsilon_vals : ndarray
        Dimensionless parameters epsilon for all particles (valid data only).
    A_vals : ndarray
        Characteristic torque coefficients A for all particles.
    """
    l_sim = np.load(f"l_sim_min_{mr}.npy")
    r_in = np.load(f"r_in_{mr}.npy")
    r_closest = np.load(f"r_closest_{mr}.npy")
    Npart = len(l_sim)
    mH = u.Msun
    mL = mr * mH
    fPBH = 1
    L = 1e-2*u.pc
    Vtilde = 4 * np.pi * (L ** (3-9/4)) / (3-9/4)
    ntilde = Npart / Vtilde
    epsilon_vals = []
    A_vals = []
    for i in range(Npart):
        rin = r_in[i]
        rmin = r_closest[i]
        l = l_sim[i]
        if rmin <= 0 or rmin >= rin or l <= 0:
            continue
        eps_orbit = rmin / rin
        # characteristic torque
        tau_c = 2.6 * u.G * mL**2 * ntilde**(2/3) * rin**(-0.5)
        # Integral I(alpha, eps)
        I_val = quad(lambda u: u**0.5 * (1-u)**(-0.5), eps_orbit, 1.0)[0]
        # A coefficient
        A = (1 / (np.sqrt(2*u.G*mH)*mL)) * tau_c * rin**(3/2) * I_val
        if A > 0:
            epsilon_vals.append(l / A)
            A_vals.append(A)
    return np.array(epsilon_vals), np.array(A_vals)

# Section: Plot P(epsilon) vs g_numeric
# ---
def plot_epsilon_distribution(mr):
    """
    Generate and save histogram of epsilon distribution vs analytical prediction.
    
    Creates a comparison plot showing the simulated epsilon distribution as a
    histogram overlaid with the analytical g_numeric(epsilon) curve.
    
    Parameters
    ----------
    mr : float
        Mass ratio (used for data file names and plot title).
    
    Generates
    ----------
    epsilon_distribution.pdf : Publication-quality PDF figure
    """
    epsilon_vals, A_vals = compute_epsilon(mr)
    # Histogram of epsilon
    epsmax = 10
    nbins = 120
    eps_vals = epsilon_vals[epsilon_vals<=epsmax]
    bins = np.linspace(eps_vals.min(), eps_vals.max(), nbins)
    epsl = np.linspace(0,epsmax, 100)
    g_vals = np.array([g_numeric(eps) for eps in epsl])
    hist, bin_edges = np.histogram(eps_vals, bins=bins, density=True)
    # Plot
    fig, ax = plt.subplots(figsize=(10,8), constrained_layout=True)
    ax.step(bin_edges[:-1], hist, where='post', linewidth=2, color="royalblue", label = "Simulated")
    plt.plot(epsl, g_vals, lw=2.5, color="crimson", label="Analytical")   
    plt.xlabel(r"$\epsilon = \ell / \mathcal{A}$", fontsize=20)
    plt.ylabel(r"$\mathcal P (\ell / \mathcal{A})$", fontsize=20)    
    ax.tick_params(which='major', direction='in', width=2, length=10, top=True, right=True, labelsize=18)
    ax.tick_params(which='minor', direction='in', width=1, length=7, top=True, right=True, pad=11)    
    ax.set_title(f"$m_r = 10^{{{np.log10(mr):.0f}}}$", fontsize=21)
    plt.xlim(0,10)
    plt.ylim(0,0.5)
    plt.legend(frameon=False, fontsize=15)
    plt.savefig("epsilon_distribution.pdf", dpi=300)
    

# Section: Run
# ---
mr = 1e-5
plot_epsilon_distribution(mr)