# practice using galario to fit 2d gaussian
import numpy as np
import galario
# galario.double.threads(16) # to set number of threads for parallelization
from galario.double_cuda import sampleImage
from radmc3dPy import *
import h5py
import matplotlib.pyplot as plt

# band 3 uvhdf5
uvfile = '/scratch/zdl3gk/data/alma/hh212_2017b3/2017.1.00712.S/science_goal.uid___A001_X1284_X17bc/group.uid___A001_X1284_X17bd/member.uid___A001_X1284_X17be/calibrated/RiverDir/data.hdf5'

fid = h5py.File(uvfile, "r")
freqs = fid["freqs"][:] # [Hz]
uu = fid["uu"][:,:] # [kilolam]
vv = fid["vv"][:,:] # [kilolam]
real = fid["real"][:,:] # [Jy]
imag = fid["imag"][:,:] # [Jy]
weight = fid["weight"][:,:] #[1/Jy^2]
flag = fid["flag"]

fid.close()

# take out the zero weigting data points
reg = weight != 0.
uu = uu[reg]
vv = vv[reg]
real = real[reg]
imag = imag[reg]
weight = weight[reg]

# determine image size
nxy, dxy = galario.double.get_image_size(uu, vv, verbose=True)

# set up easy gaussian
def GaussianProfile(f0, sigma, Rmin, dR, nR):
    """ Gaussian brightness profile. """
    # radial grid
    R = np.linspace(Rmin, Rmin + dR*nR, nR, endpoint=False)
    return f0 * np.exp(-0.5*(R/sigma)**2)

def lnpostfn(p, p_ranges, Rmin, dR, nR, nxy, dxy, u, v, Re, Im, w):
    """ Log of posterior probability function """
    lnprior = lnpriorfn(p, p_ranges)  # apply prior
    if not np.isfinite(lnprior):
        return -np.inf
    # unpack the parameters
    f0, sigma, inc, PA, dRA, dDec = p
    f0 = 10.**f0        # convert from log to real space
    # convert to radians
    sigma *= arcsec
    Rmin *= arcsec
    dR *= arcsec
    inc *= deg
    PA *= deg
    dRA *= arcsec
    dDec *= arcsec
    # compute the model brightness profile
    f = GaussianProfile(f0, sigma, Rmin, dR, nR)
    chi2 = chi2Profile(f, Rmin, dR, nxy, dxy, u, v, Re, Im, w,
                       inc=inc, PA=PA, dRA=dRA, dDec=dDec)
    return -0.5 * chi2 + lnprior

import emcee
# radial grid parameters
Rmin = 1e-4
dR = 0.01
nR = 2000

# parameter space domain
p_ranges = [[1, 20],
            [0., 8.],
            [0., 90.],
            [0., 180.],
            [-2., 2.],
            [-2., 2.]]

ndim = len(p_ranges)
nwalkers = 12
nthreads = 1
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnpostfn, 
            args=[p_ranges, Rmin, dR, nR, nxy, dxy, uu, vv, real, imag, weight],
            threads=nthreads)


from galario import deg, arcsec
from galario.double import chi2Profile

def lnpriorfn(p, par_ranges):
    """ Uniform prior probability function """
    for i in range(len(p)):
        if p[i] < par_ranges[i][0] or p[i] > par_ranges[i][1]:
            return -np.inf
    jacob = -p[0]       # jacobian of the log transformation
    return jacob

## prepare to run mcmc
nsteps = 300

# initial guess for the parameters
p0 = [10, 0.5, 70., 23., 0., 0.] 
#  3 parameters for the model + 4 (inc, PA, dRA, dDec)

# initialize the walkers with an ndim-dimensional Gaussian ball
pos = [p0 + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

# execute the MCMC
#pos, prob, state = sampler.run_mcmc(pos, nsteps, rstate0=state, lnprob0=prob)
sampler.run_mcmc(pos, nsteps)

# plot the resulting MCMC
import corner
samples = sampler.chain[:, -1000:, :].reshape((-1, ndim))
fig = corner.corner(samples, labels=["$f_0$", "$\sigma$", r"$i$", r"PA", r"$\Delta$RA", r"$\Delta$Dec"],
                    show_titles=True, quantiles=[0.16, 0.50, 0.84], label_kwargs={'labelpad':20, 'fontsize':0}, fontsize=8)
fig.savefig("triangle_example.png")
