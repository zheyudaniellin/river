# setup_river1.py
# setup sample observations of gaussian blob
import os
import matplotlib.pyplot as plt
from radmc3dPy import *
import river
import numpy as np
import galario

rundir = 'testing2/'
projdir = '/scratch/zdl3gk/coding/python/river_code/river_v2/projects/blueprint/'
os.system('rm -rf '+rundir)
os.system('mkdir '+rundir)

par = { #on/off, value, sigma, lower, upper, unit
    'inc':  [ 1, 43., 0., 0., 90., ''],
    'mdisk': [0, 0.22, 0.05, 0.1, 0.5, 'ms'],
    'sigp': [1, 0.5, 0.05, 0.1, 1.5, ''],
    'Rsig': [0, 45, 10, 20, 70, 'au'],
    'Rt': [ 0, 53, 10, 10, 70, 'au'],
    'T0': [ 1, 55, 10, 20, 70, ''],
    'qtemp':[ 1, -0.4, 0.1, -0.75, -0.1, '']
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
cup.modelname = 'disk_thin'
cup.getDefaultPar()
cup.updatemcPar(ppar=radpar)
cup.updatemcPar(ppar=waterpar.ppar)
cup.getOpacPar(specs=['Sil_Draine'])
cup.opactype='Beck'
op = cup.getmcOpac(waterpar=waterpar)
model = cup.getmcModel()
cup.impar['sizeau'] = 200
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
nuvpoints = 1000
meanuv = [0, 0]
cov = [[1., 0], [0, 1]]
#u, v = np.random.multivariate_normal(meanuv, cov, 5000).T * 1000.
uvradius = np.random.rand(nuvpoints) * (18000 - 10) + 10
uvtheta = np.random.rand(nuvpoints) * np.pi * 2.
u = uvradius * np.cos(uvtheta) #[klam]
v = uvradius * np.sin(uvtheta) #[klam]

# estimate pixel number[#] and pixel size [cm]
nxy, dxy = galario.double.get_image_size(u*1e3, v*1e3, # in lamda
    PB=1.22*870e-4/16e5, verbose=True)

# get the visibility sampling
im = image.readImage(fname=cup.imfname, dpc=400.)
net = river.galarioTool.galarioTool(rundir=rundir, projdir=projdir)
net.radmc3dImage = im
net.uu = [u]
net.vv = [v]
net.dRA = np.array([-0.05]) /3600. * np.pi / 180.
net.dDec = np.array([0.03]) / 3600. * np.pi / 180.
net.PA = 10. * np.pi / 180.
net.getVis()
uvdat = net.toUVdata(istokes=0, iwav=0)
uvdat.getAmpPhase()
uvdat.weight = (1. / (uvdat.amp*0.01)**2) * uvdat.weight
uvtab = uvdat.river2UVTable()

# output the sample data
uvdat.writeHDF5(filename=rundir+'sampleObserved.hdf5', objname='GaussianBlob')

# plotting
uvbin_size = 100e3              # uv-distance bin [lam]
axes = uvtab.plot(label='Data', uvbin_size=uvbin_size)
plt.show()

# to overplot data and model
#uvtab.plot(label='Model', uvbin_size=uvbin_size, axes=axes, yerr=False, linestyle='-', color='r')

# plot the image
image.plotImage(image=im, au=True, bunit='Tb', cmap=plt.cm.jet)
plt.show()

