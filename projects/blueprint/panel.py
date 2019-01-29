# setup_river2.py
# setup to emcee the sample observations
import os
import matplotlib.pyplot as plt
from radmc3dPy import *
import river
import numpy as np
import emcee
import galario
import fntools

rundir = 'run1/'
projdir = '/scratch/zdl3gk/coding/python/river_code/river_v2/projects/blueprint/'
os.system('rm -rf '+rundir)
os.system('mkdir '+rundir)

# read in the observation data
obsfile = 'testing2/sampleObserved.hdf5'
obsuv = river.UVdata.UVdata()
obsuv.readHDF5(filename=obsfile)

# determine the image pixel size
nxy, dxy = galario.double.get_image_size(obsuv.uu, obsuv.vv, verbose=True)

# set up parameters for emcee
par = { #on/off, value, sigma, lower, upper, unit
    'rho0': [ 0, 1e-15, 1e-15, 1e-17, 1e-13, ''],
    'T0':   [ 1, 95., 20., 10., 160., ''], 
    'radius': [ 1, 25., 5., 1, 60., 'au'], 
    'inc':  [ 0, 0., 0., 0., 90., ''], 
    'PA' :  [ 0, 0., 5., -30., 30., ''], 
    'dRA': [
            [ 0, 0., 0.01, -0.05, 0.05, '']
           ],
    'dDec': [
             [0, 0., 0.01, -0.05, 0.05, '']
            ]
      }

riverpar = river.params.RiverPar(rundir=rundir, projdir=projdir)
riverpar.getDefaultRadpar()
riverpar.dis = 400.
riverpar.opactype = 'Beck'
riverpar.specs = ['Sil_Draine']
riverpar.modelname = 'simple_1'

waterpar = river.params.WaterPar()
waterpar.par = par
waterpar.initiatePar()
ndim = len(waterpar.parC)

mill = river.emceeTool.emceeTool()
mill.waterpar = waterpar
mill.riverpar = riverpar
mill.steps = 300
mill.nwalkers = 10
mill.nthreads = 16
obsgrid = [[obsuv]]
mill.getPball()
mill.doSampler(obsgrid)
sampler = mill.sampler

# outputs
fntools.zylwritevec(sampler.chain, rundir+'dat.chain')
samples = sampler.chain.reshape((-1, ndim))

# obtained values
vals = np.zeros([ndim, 3], dtype=np.float64)
for ii in range(ndim):
    valsii = np.percentile(samples[:,ii], [16, 50, 84])
    vals[ii, :] = valsii

boat = river.resultsTool.resultsTool(riverpar=riverpar)
boat.waterpar = waterpar
boat.sampler = sampler
boat.getSamples(inxcut=2)
fig = boat.getTriangle(title_fmt='.2e')

