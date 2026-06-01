"""
PDF of beta from initial conditions, compared to the analytical prediction.

Definition
----------
beta = tau_i / tau_c

where:
- tau_i is the (vector) torque from the discrete PBH distribution in the ICs
- tau_c is the characteristic torque from the analytical estimate

This script computes and visualizes the distribution of the dimensionless torque
parameter beta, comparing simulated values (from particle-particle interactions)
against the analytical prediction from torque statistics.

Outputs
-------
- beta_distribution_<mr>.pdf : Publication-quality PDF plot
- beta_distribution_<mr>.png : High-resolution PNG plot (300 dpi)

Usage
-----
    python pdf_beta_ta.py
"""

#-------

# Section: Imports
# ---
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from scipy.special import j0

# Section: Matplotlib style
# ---
matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{mathpazo}"
plt.rcParams["axes.linewidth"] = 2
plt.rc("text", usetex=True)
plt.rc("font", family="serif")


# Section: Analytical model
# ---
def g_numeric(gamma: float) -> float:
    """
    Compute the analytical distribution of beta from torque statistics.
    
    Evaluates the predicted probability density for the dimensionless torque
    parameter beta using numerical integration.
    
    Parameters
    ----------
    gamma : float
        Dimensionless torque parameter beta.
    
    Returns
    -------
    float
        Probability density P(beta) at the given gamma.
    """
    """Analytical prediction g_numeric(gamma) used for comparison."""
    f = lambda s: s * j0(s * gamma) * np.exp(-s ** (3 / 2))
    I, _ = quad(f, 0, np.inf)
    return gamma * I


# Section: Plotting
# ---
def generate_beta_plot(mr: float) -> None:
    """
    Generate and save histogram of beta distribution vs analytical prediction.
    
    Creates a comparison plot showing the simulated beta distribution as a
    step histogram overlaid with the analytical g_numeric(beta) curve.
    Particles with extreme beta values are filtered out before plotting.
    
    Parameters
    ----------
    mr : float
        Mass ratio (used for data file names and plot title).
    
    Generates
    ----------
    beta_distribution_{mr}.pdf : PDF figure
    beta_distribution_{mr}.png : PNG figure (300 dpi)
    """
    """Load beta samples for `mr` and plot histogram vs analytical curve."""
    beta_p = np.load(f"beta_in_{mr}_new.npy")

    beta_max_plot = 10
    mask = beta_p <= beta_max_plot

    n_total = len(beta_p)
    n_kept = int(np.sum(mask))
    n_removed = n_total - n_kept

    print(f"Keeping {n_kept}/{n_total} particles")
    print(f"Removed {n_removed} particles with beta > {beta_max_plot}")

    beta_plot = beta_p[mask]

    nbins = 200
    bins = np.linspace(beta_plot.min(), beta_plot.max(), nbins)
    betas = np.linspace(0, beta_max_plot, 100)
    numerical = np.array([g_numeric(beta) for beta in betas])

    hist, bin_edges = np.histogram(beta_plot, bins=bins, density=True)

    fig, ax = plt.subplots(constrained_layout=True, figsize=(10, 8))
    ax.step(bin_edges[:-1], hist, where="post", linewidth=2, color="royalblue", label="Simulated")
    ax.plot(betas, numerical, lw=2, color="crimson", label="Analytical")

    ax.set_xlabel(r"$\beta$", fontsize=20)
    ax.set_ylabel(r"$\mathcal P (\beta)$", fontsize=20)

    ax.tick_params(which="major", direction="in", width=2, length=10, top=True, right=True, labelsize=18, pad=11)
    ax.tick_params(which="minor", direction="in", width=1, length=7, top=True, right=True, pad=11)

    ax.set_xlim(0, 10)
    ax.set_ylim(0, 0.5)
    ax.set_title(f"$m_r = 10^{{{np.log10(mr):.0f}}}$", fontsize=21)
    ax.legend(frameon=False, fontsize=15)

    plt.savefig(f"beta_distribution_{mr}.pdf")
    plt.savefig(f"beta_distribution_{mr}.png", dpi=300)


# Section: CLI entry point
# ---
def main() -> None:
    generate_beta_plot(1e-5)


if __name__ == "__main__":
    main()