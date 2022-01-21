import scipy
import numpy as np


def compute_structure_factor_2d(radii, g_r, number_density, qmax=0.1):
    dr = radii[1]-radii[0]
    Q = np.linspace(0.0, qmax, 512)
    S_q = np.zeros_like(Q)
    
    h_r = g_r - 1
    
    for i, q, in enumerate(Q):
        S_q[i] = np.sum([2*np.pi*r*h_r[j]*scipy.special.j0(q*r)*dr for j, r in enumerate(radii)])
    
    S_q = 1 + number_density*S_q
    return Q, S_q


def infinite_length_cylindrical_formfactor(q, r):
    return 2*scipy.special.jv(1, q*r) / (q*r)


def polydisperse_infinite_length_cylindrical_formfactor(q, ave_r, std_r, n=2048):
    ff = np.zeros_like(q)
    size_dist = scipy.stats.norm(loc=ave_r, scale=std_r)

    min_r = ave_r - (4*std_r)
    max_r = ave_r + (4*std_r)

    steps = np.linspace(min_r, max_r, n)
    step = steps[1]-steps[0]

    for r_i in steps:
        w_i = size_dist.pdf(r_i) * step
        ff += w_i * infinite_length_cylindrical_formfactor(q, r_i)**2
    return ff
    