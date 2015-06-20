# -*- coding: utf-8 -*-
# Copyright 2014-2015 Victor Amin, http://vamin.net/
"""
This script fits data for bound, free, and dR to a dual langmuir model.

Victor Amin 2014-2015
"""

import sys
import numpy as np
from scipy import optimize


def bnd_v_free(Ka1, Ka2, Ns, m, free):
    """Sum of two langmuirs. Returns number of ligands bound given the
    parameters:
        Ka1: first adsorption constant
        Ka2: second adsorption constant
        Ns: the number of sites for the first adsorption
        m: multiple of sites for the second adsorption relative to the sites
           for the first adsorption (i.e. Ns2 = m * Ns1)
        free: the concentration of free ligands
    """
    return (Ka1 * (Ns * free) / (1. + Ka1 * free)
            + (m * Ka2) * (Ns * free) / (1. + Ka2 * free)
            )


def dR_v_bound(Ka1, Ka2, Ns, m, dr, bound):
    """Returns the overall dR given the parameters:
        Ka1: first adsorption constant
        Ka2: second adsorption constant
        Ns: the number of sites for the first adsorption
        m: multiple of sites for the second adsorption relative to the sites
           for the first adsorption (i.e. Ns2 = m * Ns1)
        dr: apparent increase in radius per ligand for the first adsorption
            (dr for the second adsorption is assumed to be 0)
        bound: the number of ligands bound
    """
    x = bound
    return (dr
            * (1. / (2. * (Ka1 - Ka2)))
            * (x * (Ka1-Ka2)+Ka1*(Ns)+Ka2*(Ns*m)
               - ((-x * Ka1 + x * Ka2 + Ka1 * Ns)**2
               + 2. * Ka2 * (x * Ka1 - x * Ka2 + Ka1 * Ns)
               * Ns * m + Ka2**2 * ((Ns * m)**2))**0.5)
            )


def dual_site_global(fitting_parameters, bound, free, delta_R):
    """Returns the overall error of the fitting parameters relative to the
    provided data for bnd_v_free and dR_v_bound as a 1-d array.

        fitting_parameters: tuple containing Ka1, Ka2, m, Ns, and dr
        bound: number of bound ligands
        free: concentration of free ligands
        delta_R: observed increase in apparent radius
    """
    # parse fitting parameters
    Ka1, Ka2, m = fitting_parameters[:3]
    Ns = fitting_parameters[3:(len(fitting_parameters)-3)/2+3]
    dr = fitting_parameters[(len(fitting_parameters)-3)/2+3:]
    # calculate error of fit versus provided values
    bnd_v_free_err = bound - bnd_v_free(Ka1, Ka2, Ns, m, free)
    dR_v_bound_err = delta_R - dR_v_bound(Ka1, Ka2, Ns, m, dr, bound)
    # normalize bnd_v_free errors so they don't dominate over dR_v_bound errors
    err_correction_factor = bound / delta_R
    bnd_v_free_err = bnd_v_free_err / err_correction_factor
    # return errors as 1-d array
    err = np.nan_to_num(np.concatenate((bnd_v_free_err, dR_v_bound_err)))
    return np.concatenate(err)

# initialize fitting parameters
Ka1_init = 10000.
Ka2_init = 100.
m_init = 1.
Ns_init = (100., 100., 100., 100., 100., 100., 100.)
dr_init = (.001, .001, .001, .001, .001, .001, .001)
print "Initial Parameters:"
print "Ka1 = %f, Ka2 = %f, m = %f" % (Ka1_init, Ka2_init, m_init)
print "Ns ="
print Ns_init
print "dr ="
print dr_init

# read in experimental data
bound = np.genfromtxt('data/bound.csv', delimiter=',')
free = np.genfromtxt('data/free.csv', delimiter=',')
delta_R = np.genfromtxt('data/delta_R.csv', delimiter=',')

# fit data
p0 = (Ka1_init, Ka2_init, m_init) + Ns_init + dr_init
np.seterr(all='ignore')  # else, complaints about negative and missing values
fits, red_cov, info, err, flag = optimize.leastsq(dual_site_global,
                                                  p0,
                                                  args=(bound, free, delta_R),
                                                  maxfev=100000,
                                                  full_output=1)
if flag > 0 and flag < 5:
    sys.stderr.write('\nFit succeeded.\n\n')
else:
    sys.stderr.write('\nFit FAILED.\n\n')
    sys.exit()

# calculate error:
s_sq = (dual_site_global(fits,
                         bound,
                         free,
                         delta_R)**2).sum()/(np.size(delta_R)-len(p0))
cov = red_cov * s_sq
error = []
for i in xrange(len(fits)):
    try:
        error.append(np.absolute(cov[i][i])**0.5)
    except:
        error.append(0.0)

# output fits
Ka1_out, Ka2_out, m_out = fits[:3]
Ka1_err, Ka2_err, m_err = error[:3]
Ns_out = fits[3:(len(fits)-3)/2+3]
Ns_err = error[3:(len(fits)-3)/2+3]
dr_out = fits[(len(fits)-3)/2+3:]
dr_err = error[(len(fits)-3)/2+3:]
print "\nFinal Parameters:"
print "Ka1 = %f +/- %f" % (Ka1_out, Ka1_err)
print "Ka2 = %f +/- %f" % (Ka2_out, Ka2_err)
print "m = %f +/- %f" % (m_out, m_err)
print "Ns ="
print '\n'.join(map(lambda x: '%f +/- %f' % x, zip(Ns_out, Ns_err)))
print "dr ="
print '\n'.join(map(lambda x: '%f +/- %f' % x, zip(dr_out, dr_err)))

print
print 'Reduced chi-sqared: %f' % s_sq
print

print 'Covariance matrix (csv):'
for i, row in enumerate(cov):
    for j in xrange(len(p0)):
        sys.stdout.write('%f' % cov[i,j])
        if j < len(p0)-1:
            print ',',
print '\n'

print 'Correlation matrix (csv):'
for i, row in enumerate(cov):
    for j in xrange(len(p0)):
        sys.stdout.write('%f' % (cov[i,j]/np.sqrt(cov[i,i]*cov[j,j])))
        if j < len(p0)-1:
            print ',',
    print
