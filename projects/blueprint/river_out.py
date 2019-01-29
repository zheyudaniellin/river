# river_out.py
# output plots

import os
import matplotlib.pyplot as plt
from radmc3dPy import *
import river
import numpy as np
import emcee
import galario
import fntools

rundir = 'samplerun4/'
projdir = '/scratch/zdl3gk/coding/python/river_code/river_v2/projects/blueprint/'

# read in the observation data
obsfile = 'testing2/sampleObserved.hdf5'
obsuv = river.UVdata.UVdata()
obsuv.readHDF5(filename=obsfile)

# read in riverpar
# read in waterpar

#------------------ outputs ---------------------------
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

# get the final value
chi, cup, im, net = boat.getBestModel(obsgrid)

# get the grid of model uv
modelgrid = net.toUVdata()

# plotting uv-amplitude
boat.plotUVamp(obsgrid, modelgrid, fdir=rundir, uvbin_size=100e3)

# plot the image residuals

