# setup_exampleUV.py
# use the uv coverage to sample the example image of a gaussian blob
import os
import numpy as np
from river import *
from radmc3dPy import *
import galario
import matplotlib.pyplot as plt
import uvplot

uvfile = '/scratch/zdl3gk/data/alma/hh212_2017b6/2017.1.00712.S/science_goal.uid___A001_X1284_X17b6/group.uid___A001_X1284_X17b7/member.uid___A001_X1284_X17b8/calibrated/RiverDir/data.hdf5'

imdata = 'myimage.i0.out'

# read in the uvfile
uvdat = UVdata.UVdata()
uvdat.readHDF5(filename=uvfile)
uvdat.applyflag()
uvdat.getAmpPhase()
reg = uvdat.weight<6300.
uvdat.applyflag(flagreg=reg)
reg = uvdat.weight>6400.
uvdat.applyflag(flagreg = reg)

im = image.readImage(fname=imdata)
b7 = im.image[:,:,0,0]
b7c = b7.copy(order='C') #can check with b7c.flags
dis = 400.
dxy = (im.sizepix_x / natconst.au / dis) / 3600.
dRA = 0.
dDec = 0.
PA = 0.* natconst.rad
vis = galario.double.sampleImage(b7c, dxy, uvdat.uu/1e3, uvdat.vv/1e3,
        dRA=dRA, dDec=dDec, PA=PA)

visdat = UVdata.UVdata()
visdat.freqs = np.array([870.], dtype=np.float64)
visdat.uu = uvdat.uu
visdat.vv = uvdat.vv
visdat.real = np.real(vis)
visdat.imag = np.imag(vis)
visdat.weight = np.ones(vis.shape)
visdat.flag = np.zeros(vis.shape, dtype=np.int)
visdat.writeHDF5(filename='sampleObserved.hdf5', objname='HH212')
vistab = visdat.river2UVTable()

uvbin_size = 100e3
axes = uvtab.plot(labvel='Data', uvbin_size=uvbin_size)
plt.show()
