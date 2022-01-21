import time
import numpy as np
import scipy.stats
from scipy.spatial import cKDTree


def get_volume_fraction(segmentation):
    particles_mask = ((segmentation > 0.0) | (np.bitwise_not(np.isnan(segmentation)))).astype(np.uint8)
    return np.sum(particles_mask) / np.prod(segmentation.shape)


def get_areas(segmentation):
    areas = []
    for inst in np.unique(segmentation):
        if inst == 0 or not np.isnan(inst):
            inst_mask = (segmentation == inst).astype(np.uint8)
            areas.append(np.sum(inst_mask))
    return np.array(areas)


def rdf2d_with_3d_vol(particles, dr, rho=None, rcutoff=0.8, eps=1e-15, progress=False):
    """
    Computes 3D radial distribution functions g(r) from a set of particles of shape (N, 3)
    that are sitting on a plane. Adapted from rdfpy.

    Parameters
    ----------
    particles : (N, 3) np.array
        Set of particle from which to compute the radial distribution function g(r).
    dr : float
        Delta r. Determines the spacing between successive radii over which g(r) is computed.
    rho : float, optional
        Number density. If left as None, box dimensions will be inferred from the
        particles and the number density will be calculated accordingly.
    rcutoff : float
        radii cutoff value between 0 and 1. The default value of 0.8 means the independent
        variable (radius) over which the RDF is computed will range from 0 to 0.8*r_max. 
        This removes the noise that occurs at r values close to r_max, due to fewer valid 
        particles available to compute the RDF from at these r values.
    eps : float, optional
        Epsilon value used to find particles less than or equal to a distance in KDTree.
    progress : bool, optional
        Set to False to disable progress readout.

    FROM COMMIT: https://github.com/by256/rdfpy/blob/c766dc8cd5543ffe6d1703405dc8ee4ad421d149/rdfpy/rdfpy.py


    Returns
    -------
    g_r : (n_radii) np.array
        radial distribution function values g(r).
    radii : (n_radii) np.array
        radii over which g(r) is computed
    """

    start = time.time()

    # translate particles such that the particle with min coords is at origin
    particles = particles - np.min(particles, axis=0)
    min_x, min_y = np.min(particles, axis=0)
    max_x, max_y = np.max(particles, axis=0)

    # dimensions of box
    w, h = (max_x-min_x), (max_y-min_y)

    r_max = (np.min([w, h]) / 2)*rcutoff
    radii = np.arange(dr, r_max, dr)
    g_r = np.zeros(shape=(len(radii)))

    N = len(particles)
    if not rho:
        l = np.sqrt(w*h)
        rho = N / (w*h*l)  # number density

    # create a KDTree for fast nearest-neighbor lookup of particles
    tree = cKDTree(particles)

    for r_idx, r in enumerate(radii):
        # find all particles that are at least r + dr away from the edges of the box
        valid_idxs = (particles[:, 0]-(r+dr) >= min_x) & (particles[:, 0]+(r+dr) <= max_x) \
            & (particles[:, 1]-(r+dr) >= min_y) & (particles[:, 1]+(r+dr) <= max_y)
        valid_particles = particles[valid_idxs]

        # compute n_i(r) for valid particles.
        for particle in valid_particles:
            n = tree.query_ball_point(particle, r+dr-eps, return_length=True) - \
                tree.query_ball_point(particle, r, return_length=True)
            g_r[r_idx] += n

        # normalize
        n_valid = len(valid_particles)
        shell_vol = (4/3)*np.pi*((r+dr)**3 - r**3)
        g_r[r_idx] /= n_valid*shell_vol*rho

        if progress:
            print('Computing RDF     Radius {}/{}    Time elapsed: {:.3f} s'.format(
                r_idx+1, len(radii), time.time()-start), end='\r', flush=True)

    return g_r, radii, rho


def compute_structure_factor_3d(g_r, radii, rho, qmax=1.0, n=512):
    """
    Compute structure factor S(q) from a radial distribution function g(r).
    
    Parameters
    ----------
    
    g_r : np.ndarray
        Radial distribution function, g(r).
    radii : np.ndarray
        Independent variable of g(r).
    rho : float
        Average number density of particles.
    qmax : float
        Maximum value of momentum transfer (the independent variable of S(q)).
    n : int
        Number of points in S(q).
        
    Returns
    -------
    S_q : np.ndarray
        Structure factor
    Q : np.ndarray
        Momentum transfer (the independent variable of S(q)).
    
    """
    n_r = len(g_r)
    Q = np.linspace(0.0, qmax, n)
    S_q = np.zeros_like(Q)
    
    dr = radii[1] - radii[0]
    h_r = g_r - 1
    
    for q_idx, q in enumerate(Q):
        if q_idx == 0:
            S_q[q_idx] = np.sum([4*np.pi*(r**2)*h_r[r_idx]*dr for r_idx, r in enumerate(radii)])
        else:
            S_q[q_idx] = np.sum([4*np.pi*(r**2)*h_r[r_idx]*dr*np.sin(q*r)/(q*r) for r_idx, r in enumerate(radii)])
    
    S_q = 1 + rho * S_q / n_r
    
    return S_q, Q


def spherical_formfactor(q, r):
    return 3*(np.sin(q*r) - (q*r)*np.cos(q*r))/(q*r)**3


def polydisperse_spherical_formfactor(q, ave_r, std_r, vol=False, n=4096):
    ff = np.zeros_like(q)
    size_dist = scipy.stats.norm(loc=ave_r, scale=std_r)

    min_r = ave_r - (6*std_r)
    max_r = ave_r + (6*std_r)

    steps = np.linspace(min_r, max_r, n)
    step = steps[1]-steps[0]

    for r_i in steps:
        w_i = size_dist.pdf(r_i) * step
        if vol:
            current_vol = (4/3)*np.pi*r_i**3
            ff += w_i*spherical_formfactor(q, r_i)**2*current_vol
        else:
            ff += w_i*spherical_formfactor(q, r_i)**2
    return ff

