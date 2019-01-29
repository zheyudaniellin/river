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

rundir = 'testing3/'
projdir = '/scratch/zdl3gk/coding/python/river_code/river_v2/projects/blueprint/'
os.system('rm -rf '+rundir)
os.system('mkdir '+rundir)

# read in the observation data
obsfile = 'testing2/sampleObserved.hdf5'
obsuv = river.UVdata.UVdata()
obsuv.readHDF5(filename=obsfile)

# determine the image pixel size
nxy, dxy = galario.double.get_image_size(obsuv.uu*1e3, obsuv.vv*1e3, verbose=True)

# set up parameters for emcee
par = { #on/off, value, sigma, lower, upper, unit
    'inc':  [ 1, 43., 2., 0., 90., ''],
    'PA' :  [ 1, 5, 5., -30., 30., ''],
    'dRA': [
            [ 1, -0.02, 0.02, -0.15, 0.15, '']
           ],
    'dDec': [
             [1, 0., 0.02, -0.1, 0.1, '']
            ], 
    'mdisk': [0, 0.22, 0.1, 0., 0.8, 'ms'],
    'sigp': [0, 0.5, 0.05, 0., 1.5, ''],
    'Rsig': [0, 45, 5, 0, 100, 'au'],
    'Rt': [ 0, 53, 10, 0, 100, 'au'],
    'T0': [ 0, 55, 10, 0, 100, ''],
    'qtemp':[ 0, -0.4, 0.1, -0.75, -0.1, '']
      }

riverpar = river.params.RiverPar(rundir=rundir, projdir=projdir)
riverpar.getDefaultRadpar()
riverpar.dis = 400.
riverpar.opactype = 'Beck'
riverpar.specs = ['Sil_Draine']
riverpar.modelname = 'disk_thin'

waterpar = river.params.WaterPar()
waterpar.par = par
waterpar.initiatePar()
ndim = len(waterpar.parC)

mill = river.emceeTool.emceeTool()
mill.waterpar = waterpar
mill.riverpar = riverpar
mill.steps = 300
mill.nwalkers = 8
mill.nthreads = 16
obsgrid = [[obsuv]]
mill.getPball()
mill.doSampler(obsgrid)
sampler = mill.sampler

#------------------ outputs ---------------------------
fntools.zylwritevec(sampler.chain, rundir+'dat.chain')
samples = sampler.chain.reshape((-1, ndim))

boat = river.resultsTool.resultsTool(riverpar=riverpar)
boat.waterpar = waterpar

# read in chain and samples
boat.getSamples(fdir=rundir)

# obtained values
boat.getFitval()

# triangle
fig = boat.getTriangle(title_fmt='.2e')

# values per step
fig,axes = boat.getChiStep()

# get the final best model
chi, cup, im, net = boat.getBestModel(obsgrid)

# get the grid of model uv
modelgrid = net.toUVdata()

# plotting uv-amplitude
boat.plotUVamp(obsgrid, modelgrid, fdir=rundir, uvbin_size=100e3)

# plot the image residuals
