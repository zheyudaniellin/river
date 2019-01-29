import matplotlib.pyplot as plt
import numpy as np
import uvplot
import zylconst
import h5py
import river

#uvfile = '/scratch/zdl3gk/data/alma/hh212_2017b3/2017.1.00712.S/science_goal.uid___A001_X1284_X17bc/group.uid___A001_X1284_X17bd/member.uid___A001_X1284_X17be/calibrated/RiverDir/data.hdf5'
uvfile = 'avgall0.hdf5'

uvdat = river.UVdata.UVdata()
uvdat.readHDF5(filename=uvfile)
uvdat.applyflag()
uvdat.getAmpPhase()

# the uvtable object of uvplot
uvtab = uvdat.river2UVTable()

# some projection etc
dRA = 0.0 * uvplot.arcsec 	# Delta Right Ascension offset [rad]
dDEC = 0.0 * uvplot.arcsec	# Delta Declination offset [rad]

inc = zylconst.rad * 0.	# inclination [rad]
PA = zylconst.rad * 0. 	# position Angle [rad]

uvtab.apply_phase(dRA, dDEC)
uvtab.deproject(inc, PA)

# plotting
uvbin_size = 100e3              # uv-distance bin [lam]
axes = uvtab.plot(label='Data', uvbin_size=uvbin_size)
plt.savefig('uvtabplot.png')
plt.show()

# plot the uvdist to amp
uvdist = np.sqrt(uvdat.uu**2 + uvdat.vv**2)
plt.plot(uvdist, uvdat.amp, 'k.')
plt.savefig('uvdist2amp.png')
plt.show()
