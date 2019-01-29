import numpy as np
import matplotlib.pyplot as plt
import river
from radmc3dPy import *
import galario

uvfile = '/scratch/zdl3gk/data/alma/hh212_2017b3/2017.1.00712.S/science_goal.uid___A001_X1284_X17bc/group.uid___A001_X1284_X17bd/member.uid___A001_X1284_X17be/calibrated/RiverDir/data.hdf5'


uvdat = river.UVdata.UVdata()
uvdat.readhdf5(filename=uvfile)
uvdat.applyflag()
uvdat.getAmpPhase()

fname = 'myimage.i45.out'
im = image.readImage(fname=fname)
b7 = im.image[:,:,0,1]
b7c = b7.copy(order='C') #can check with b7c.flags
dis = 400.
dxy = (im.sizepix_x / natconst.au / dis) / 3600.
dRA = 0.
dDec = 0.
PA = 113.* natconst.rad
vis = galario.double.sampleImage(b7c, dxy, uvdat.uu/1e3, uvdat.vv/1e3, 
        dRA=dRA, dDec=dDec, PA=PA)

uvdist = np.sqrt(uvdat.uu**2. + uvdat.vv**2.)

