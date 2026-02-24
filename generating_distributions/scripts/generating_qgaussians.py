import argparse
import pickle
import numpy as np
from scipy.special import gamma
from scipy.integrate import cumulative_trapezoid, quad
import matplotlib
import matplotlib.pyplot as plt

import xobjects as xo
import xdeps as xd
import xtrack as xt
import xpart as xp
import xfields as xf
import xcoll as xc
import xsuite as xs


parser = argparse.ArgumentParser()
parser.add_argument(
    "mode", 
    type=str, 
    help="'4D' for non-factorisable, '2D' for factorisable"
)
parser.add_argument(
    "q", 
    type=float, 
    help="q value of q-Gaussian distribution"
)
parser.add_argument(
    "beta", 
    type=float, 
    help="beta value of q-Gaussian distribution"
)


def Cq(q):
    """
    Normalisation for q-Gaussian function.

    Parameters
    ----------
    q : float
        q parameter (1 < q < 3).

    Returns
    -------
    float
        Normalisation factor.
    """

    if q <= 1 or q >= 3:
        raise ValueError("This implementation only supports 1 < q < 3!")
    else:
        return np.sqrt(np.pi) * gamma((3-q)/(2*(q-1))) / (np.sqrt(q-1) * gamma(1/(q-1)))
    

def eq(x, q):
    """
    q-exponention function.

    Parameters
    ----------
    x : float | ndarray
        Value(s) to calculate q-exponential for.
    q : float
        q parameter (1 < q < 3).

    Returns
    -------
    float | ndarray
        Value(s) of q-exponential for given x.
    """

    if q <= 1 or q >= 3:
        raise ValueError("This implementation only supports 1 < q < 3!")
    else:
        return np.where(1+(1-q)*x>0, (1+(1-q)*x)**(1/(1-q)), 0)


def q_Gaussian(x, q, beta):
    """
    q-Gaussian function.

    Parameters
    ----------
    x : float | ndarray
        Value(s) to calculate q-Gaussian for.
    q : float
        q parameter (1 < q < 3).
    beta : float
        Scale parameter.
    
    Returns
    -------
    float | ndarray
        Value(s) of q-Gaussian for given x.
    """

    return np.sqrt(beta) * eq(-beta*x**2, q) / Cq(q)


def f2D_r(r, q, beta):
    """
    2D density function obtained from a 1D q-Gaussian profile. Based on Eq 20 in
    https://cds.cern.ch/record/2912366/files/document.pdf

    Parameters
    ----------
    r : float | ndarray
        Radius in 2D space.
    q : float
        q parameter (1 < q < 3).
    beta : float
        Scale parameter.

    Returns
    -------
    float | ndarray
    """

    if not (1 < q < 3):
        raise ValueError("Require 1 < q < 3")

    A = - beta**1.5 * (q - 3) * np.sqrt(q - 1) * r / (2 * np.pi)
    B = (1 + 1/(beta*(q-1)*r**2))**((q+1)/(2-2*q))
    C = (beta*(q-1)*r**2)**(q/(1-q))

    return A * B * C


def f4D_m(m, q, beta):
    """
    4D density function obtained from a 1D q-Gaussian profile. Based on Eq 21 in
    https://cds.cern.ch/record/2912366/files/document.pdf

    Parameters
    ----------
    m : float | ndarray
        Radius in 4D space.
    q : float
        q parameter (1 < q < 3).
    beta : float
        Scale parameter.

    Returns
    -------
    float | ndarray
    """

    if not (1 < q < 3):
        raise ValueError("Require 1 < q < 3")

    pref = (
        -beta * (q-3) * (q+1)
        * (1/(beta*(q-1)))**(1/(1-q))
        / (4*np.pi**2 * m**3 * gamma(1/(q-1)))
        * (beta*(q-1))**((q+1)/(2-2*q))
        * gamma(q/(q-1))
    )

    power1 = (1 + 1/(beta*m**2*(q-1)))**(1/(1-q) - 1.5)
    power2 = (beta*m**2*(q-1))**(1/(1-q))

    return pref * power1 * power2


def sample_radius(pdf_func, q, beta, D, N, rmax=50, gridsize=200000):
    """
    Sampling function to sample radius from custom PDF.

    Parameters
    ----------
    pdf_func : function
        Custom PDF to sample from.
    q : float
        q parameter (1 < q < 3).
    beta : float
        Scale parameter.
    D : int
        Dimension (2 or 4).
    N : int
        Number of points to sample.
    rmax : float
        Maximum radius. Default is ``50`` sigma.
    gridsize : int
        Number of grid points. Default is ``200000``.
    

    Returns
    -------
    ndarray
        Array of samples.
    """

    r = np.linspace(1e-10, rmax, gridsize)

    # radial probability density p(r) proportional to r^(D-1) f_D(r)
    if D == 2:
        surface = 2*np.pi*r
    elif D == 4:
        surface = 2*np.pi**2 * r**3
    else:
        raise ValueError("Only D=2 or D=4 supported")

    pdf = surface * pdf_func(r, q, beta)
    pdf[pdf < 0] = 0

    cdf = cumulative_trapezoid(pdf, r, initial=0)
    cdf /= cdf[-1]

    u = np.random.rand(N)

    return np.interp(u, cdf, r)


def random_direction(D, N):
    v = np.random.normal(size=(N, D))
    v /= np.linalg.norm(v, axis=1)[:, None]
    
    return v


def generate_4D_nonfactorizable(N, q, beta):
    m = sample_radius(f4D_m, q, beta, D=4, N=N)
    directions = random_direction(4, N)
    res = directions * m[:, None]

    return res[:, 0], res[:, 1], res[:, 2], res[:, 3]


def generate_2D_factorizable(N, q, beta):

    r_x = sample_radius(f2D_r, q, beta, D=2, N=N)
    r_y = sample_radius(f2D_r, q, beta, D=2, N=N)

    theta_x = 2*np.pi*np.random.rand(N)
    theta_y = 2*np.pi*np.random.rand(N)

    x  = r_x * np.cos(theta_x)
    px = r_x * np.sin(theta_x)

    y  = r_y * np.cos(theta_y)
    py = r_y * np.sin(theta_y)

    return x, px, y, py


def sample_q_gaussian_nd(N, q=1.2, beta=1.0, D=4):
    """
    Generate N samples from a D-dimensional q-Gaussian (q>1).
    
    Parameters
    ----------
    N : int
        Number of particles.
    q : float
        q parameter (1 < q < 3).
    beta : float
        Scale parameter (similar to inverse variance).
    D : int
        Dimension (4 for full 4D, 2 for one transverse plane at a time).
        
    Returns
    -------
    samples : ndarray shape (N, D)
    """
    
    if q <= 1:
        raise ValueError("This implementation assumes q > 1.")

    # Degrees of freedom equivalent (Student-t representation)
    nu = 2.0 / (q - 1.0) - D
    if nu <= 0:
        raise ValueError("Invalid q for given dimension (nu <= 0). Increase q or reduce D.")
    
    # Gaussian core
    Z = np.random.normal(size=(N, D))
    
    # Chi-square variable
    S = np.random.chisquare(df=nu, size=N)
    
    # Scaling
    scale = np.sqrt(nu / S)[:, None]
    
    samples = Z * scale
    
    # Apply beta scaling
    samples /= np.sqrt(beta)
    
    return samples


def generate_beam(N, q=1.2, beta=1.0, mode="4D"):
    """
    Generate normalised x, px, y, py particle coordinates.
    
    Parameters
    ----------
    N : int
        Number of particles.
    q : float
        q parameter (1 < q < 3).
    beta : float
        Scale parameter (similar to inverse variance).
    mode : str
        Can be "4D" for non-factorizable distribution (single 4D q-Gaussian) or 
        "2D" for a factorizable distribution (two independent 2D q-Gaussians).
        
    Returns
    -------
    x : ndarray shape (N,)
    px : ndarray shape (N,)
    y: ndarray shape (N,)
    py: ndarray shape (N,)
    """

    if mode == "4D":
        # Single 4D q-Gaussian
        samples = sample_q_gaussian_nd(N, q=q, beta=beta, D=4)
        return samples[:, 0], samples[:, 1], samples[:, 2], samples[:, 3]

    elif mode == "2D":
        # Two independent 2D q-Gaussians
        xpx = sample_q_gaussian_nd(N, q=q, beta=beta, D=2)
        ypy = sample_q_gaussian_nd(N, q=q, beta=beta, D=2)
        
        return xpx[:, 0], xpx[:, 1], ypy[:, 0], ypy[:, 1]
    
    else:
        raise ValueError("Mode must be '4D' or '2D'!")


def main():
    '''
    Print xsuite package versions
    '''
    print(f"xobjects: {xo.__version__}")
    print(f"xdeps: {xd.__version__}")
    print(f"xtrack: {xt.__version__}")
    print(f"xpart: {xp.__version__}")
    print(f"xfields: {xf.__version__}")
    print(f"xcoll: {xc.__version__}")
    print(f"xsuite: {xs.__version__}")


    '''
    Parse arguments
    '''
    args = parser.parse_args()
    mode = args.mode
    q = args.q
    beta = args.beta
    N = 100000
    bins = 100
    plot_range = [-8, 8]
    xx = np.linspace(plot_range[0], plot_range[1], N)


    '''
    Generate coordinates
    '''
    if mode == "4D":
        x, px, y, py = generate_4D_nonfactorizable(N, q, beta)
    elif mode == "2D":
        x, px, y, py = generate_2D_factorizable(N, q, beta)


    '''
    Plot
    '''
    cmap = matplotlib.colormaps["viridis"]
    cmap.set_under('white')

    fs = 20

    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    ax00 = axes[0, 0]
    ax01 = axes[0, 1]
    ax02 = axes[0, 2]
    ax03 = axes[0, 3]
    ax10 = axes[1, 0]
    ax11 = axes[1, 1]
    ax12 = axes[1, 2]
    ax13 = axes[1, 3]

    ax00.hist2d(
        x, y, 
        bins=bins, 
        range=[plot_range, plot_range], 
        cmap=cmap, 
        vmin=1
    )
    ax00.set_xlabel(r"$x$ [$\sigma$]", fontsize=fs)
    ax00.set_ylabel(r"$y$ [$\sigma$]", fontsize=fs)
    ax00.tick_params(axis='both', labelsize=fs-2)

    ax01.hist2d(
        x, px, 
        bins=bins, 
        range=[plot_range, plot_range], 
        cmap=cmap, 
        vmin=1
    )
    ax01.set_xlabel(r"$x$ [$\sigma$]", fontsize=fs)
    ax01.set_ylabel(r"$p_x$ [$\sigma$]", fontsize=fs)
    ax01.tick_params(axis='both', labelsize=fs-2)

    ax02.hist2d(
        y, py, 
        bins=bins, 
        range=[plot_range, plot_range], 
        cmap=cmap, 
        vmin=1
    )
    ax02.set_xlabel(r"$y$ [$\sigma$]", fontsize=fs)
    ax02.set_ylabel(r"$p_y$ [$\sigma$]", fontsize=fs)
    ax02.tick_params(axis='both', labelsize=fs-2)

    ax03.hist2d(
        px, py, 
        bins=bins, 
        range=[plot_range, plot_range], 
        cmap=cmap, 
        vmin=1
    )
    ax03.set_xlabel(r"$p_x$ [$\sigma$]", fontsize=fs)
    ax03.set_ylabel(r"$p_y$ [$\sigma$]", fontsize=fs)
    ax03.tick_params(axis='both', labelsize=fs-2)

    ax10.hist(
        x, 
        bins=bins, 
        range=plot_range, 
        color="black", 
        histtype="bar", 
        density=True
    )
    ax10.plot(
        xx, 
        q_Gaussian(xx, q, beta), 
        color="red", 
        label=r"$q=$" + f"{q}\n" + r"$\beta=$" + f"{beta}"
    )
    ax10.set_xlabel(r"$x$ [$\sigma$]", fontsize=fs)
    ax10.set_ylabel("density", fontsize=fs)
    ax10.tick_params(axis='both', labelsize=fs-2)
    ax10.legend(loc="lower center", fontsize=fs-2)

    ax11.hist(
        px, 
        bins=bins, 
        range=plot_range, 
        color="black", 
        histtype="bar", 
        density=True
    )
    ax11.plot(
        xx, 
        q_Gaussian(xx, q, beta), 
        color="red", 
        label=r"$q=$" + f"{q}\n" + r"$\beta=$" + f"{beta}"
    )
    ax11.set_xlabel(r"$p_x$ [$\sigma$]", fontsize=fs)
    ax11.set_ylabel("density", fontsize=fs)
    ax11.tick_params(axis='both', labelsize=fs-2)
    ax11.legend(loc="lower center", fontsize=fs-2)

    ax12.hist(
        y, 
        bins=bins, 
        range=plot_range, 
        color="black", 
        histtype="bar", 
        density=True
    )
    ax12.plot(
        xx, 
        q_Gaussian(xx, q, beta), 
        color="red", 
        label=r"$q=$" + f"{q}\n" + r"$\beta=$" + f"{beta}"
    )
    ax12.set_xlabel(r"$y$ [$\sigma$]", fontsize=fs)
    ax12.set_ylabel("density", fontsize=fs)
    ax12.tick_params(axis='both', labelsize=fs-2)
    ax12.legend(loc="lower center", fontsize=fs-2)

    ax13.hist(
        py, 
        bins=bins, 
        range=plot_range, 
        color="black", 
        histtype="bar", 
        density=True
    )
    ax13.plot(
        xx, 
        q_Gaussian(xx, q, beta), 
        color="red", 
        label=r"$q=$" + f"{q}\n" + r"$\beta=$" + f"{beta}"
    )
    ax13.set_xlabel(r"$p_y$ [$\sigma$]", fontsize=fs)
    ax13.set_ylabel("density", fontsize=fs)
    ax13.tick_params(axis='both', labelsize=fs-2)
    ax13.legend(loc="lower center", fontsize=fs-2)

    fig.subplots_adjust(wspace=0.5, hspace=0.3)

    plt.savefig(f"./results/{mode}_q{q}_beta{beta}.png", dpi=150, format='png', bbox_inches='tight')


if __name__ == "__main__":
    main()