# setup_river1.py
# setup sample observations of gaussian blob
import os
import matplotlib.pyplot as plt
from radmc3dPy import *
import river
import numpy as np

rundir = 'sampleObs/'
projdir = '/scratch/zdl3gk/coding/python/river_code/river_v2/projects/blueprint/'
os.system('rm -rf '+rundir)
os.system('mkdir '+rundir)

par = { #on/off, value, sigma, lower, upper, unit
    'rho0': [ 0, 1e-15, 5e-16, 3e-16, 5e-15, ''],
    'T0':   [ 1, 101., 10., 80., 120., ''], 
    'radius':[ 1, 19., 5., 15, 30., 'au'], 
    'inc':  [ 0, 0., 0., 0., 90., '']
      }
radpar = {
        'scattering_mode_max':0, 'alignment_mode':0, 
        'setthreads':16, 
        'incl_lines':0, 'incl_dust':1,
        'nphot':1e6, 'nphot_scat':1e6, 'modified_random_walk':1,
        'mc_scat_maxtauabs':15
         }

waterpar = river.params.WaterPar()
waterpar.par = par
waterpar.initiatePar()

cup = river.radmc3dTool.radmc3dTool(rundir=rundir, projdir=projdir)
cup.modelname = 'simple_1'
cup.getDefaultPar()
cup.updatemcPar(ppar=radpar)
cup.updatemcPar(ppar=waterpar.ppar)
cup.getOpacPar(specs=['Sil_Draine'])
cup.opactype='Beck'
op = cup.getmcOpac(waterpar=waterpar)
model = cup.getmcModel()
cup.getmcImage(waterpar=waterpar)


# get the real data
#uvfile = '/scratch/zdl3gk/data/alma/hh212_2017b3/2017.1.00712.S/science_goal.uid___A001_X1284_X17bc/group.uid___A001_X1284_X17bd/member.uid___A001_X1284_X17be/calibrated/RiverDir/data.hdf5'

#uvobs = river.UVdata.UVdata()
#uvobs.readHDF5(filename=uvfile)
#uvobs.applyflag()
#reg = (uvobs.weight > 10000.) & (uvobs.weight < 12000)
#uvobs.applyflag(flagreg=reg)
#uvobs.getAmpPhase()
#obstab = uvobs.river2UVTable()
#uvbin_size = 100e3
#axes = obstab.plot(label='Data', uvbin_size=uvbin_size)

# example u, v distribution
nuvpoints = 5000
meanuv = [0, 0]
cov = [[1., 0], [0, 1]]
#u, v = np.random.multivariate_normal(meanuv, cov, 5000).T * 1000.
uvradius = np.random.rand(nuvpoints) * (5000.-5.) + 5.
uvtheta = np.random.rand(nuvpoints) * np.pi * 2.
u = uvradius * np.cos(uvtheta) 
v = uvradius * np.sin(uvtheta)

# get the visibility sampling
im = image.readImage(fname=cup.imfname, dpc=400.)
net = river.galarioTool.galarioTool(rundir=rundir, projdir=projdir)
net.radmc3dImage = im
net.uu = [u]
net.vv = [v]
net.dRA = [0.]
net.dDec = [0.]
net.PA = 0.
net.getVis()
uvdat = net.toUVdata(istokes=0, iwav=0)
uvdat.getAmpPhase()
uvdat.weight = max(uvdat.amp) * 10.**4 * uvdat.weight
uvtab = uvdat.river2UVTable()

# output the sample data
uvdat.writeHDF5(filename=rundir+'sampleObserved.hdf5', objname='GaussianBlob')

# plotting
uvbin_size = 100e3              # uv-distance bin [lam]
axes = uvtab.plot(label='Data', uvbin_size=uvbin_size)
plt.show()

# to overplot data and model
#uvtab.plot(label='Model', uvbin_size=uvbin_size, axes=axes, yerr=False, linestyle='-', color='r')

