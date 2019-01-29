# setup_river3.py
# setup to emcee the real observations
import os
import matplotlib.pyplot as plt
from radmc3dPy import *
import river
import numpy as np
import emcee
import galario
import fntools

rundir = 'testing4/'
projdir = '/scratch/zdl3gk/coding/python/river_code/river_v2/projects/blueprint/'
os.system('rm -rf '+rundir)
os.system('mkdir '+rundir)

# read in the observation data
obsfile = 'avgall0.hdf5'
obsuv = river.UVdata.UVdata()
obsuv.readHDF5(filename=obsfile)

obsgrid = [[obsuv]]

# determine the image pixel size
nxy, dxy = galario.double.get_image_size(obsuv.uu*1e3, obsuv.vv*1e3, verbose=True)

# set up parameters for emcee
par = { #on/off, value, sigma, lower, upper, unit
    'inc':  [ 0, 86., 2., 0., 90., ''],
    'PA' :  [ 1, 113, 10., 0., 180, ''],
    'dRA': [
            [ 1, 0., 0.02, -0.15, 0.15, '']
           ],
    'dDec': [
             [1, 0., 0.02, -0.1, 0.1, '']
            ], 
    'mdisk': [1, 0.05, 0.01, 0.01, 0.8, 'ms'],
    'sigp': [1, 0.5, 0.05, 0., 1.5, ''],
    'Rsig': [1, 45, 5, 0, 100, 'au'],
    'Rt': [ 1, 20, 10, 0, 100, 'au'],
    'T0mid': [ 1, 105, 20, 0, 500, ''],
    'T0atm': [ 1, 130, 20, 0, 500, ''],
    'H0':    [ 1, 7.5, 2., 0., 25., 'au'],
    'Hd' :   [ 1, 5., 1., 0., 10., '']
      }

riverpar = river.params.RiverPar(rundir=rundir, projdir=projdir)
riverpar.getDefaultRadpar()
riverpar.dis = 400.
riverpar.opactype = 'Beck'
riverpar.specs = ['Sil_Draine']
riverpar.modelname = 'zyl1'

waterpar = river.params.WaterPar()
waterpar.par = par
waterpar.initiatePar()
ndim = len(waterpar.parC)

mill = river.emceeTool.emceeTool()
mill.waterpar = waterpar
mill.riverpar = riverpar
mill.steps = 100
mill.nwalkers = 40
mill.nthreads = 16
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
boat.getSamples(fdir=rundir)
fig = boat.getTriangle(title_fmt='.2e')
fig,axes = boat.getChiStep()
