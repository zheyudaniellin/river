# UVdata.py
# class to manipulate uv data
# assume we are reading/writing in uvhdf5
import os
import numpy as np
import h5py
import uvplot

class UVdata(object):
    """UVdat
    """
    def __init__(self, filename=None):
        """
        Initial UVdat object
        Attributes
        freqs 	: frequencies in Hz
        uu 	: u in klambda
        vv	: v in klambda
        real	: in Jy
        imag	: in Jy
        weight 	: in 1 / Jy^2
        flag	: in boolean
        ----------
        """
        self.freqs = np.array([], dtype=np.float64)
        self.uu = np.array([], dtype=np.float64)
        self.vv = np.array([], dtype=np.float64)
        self.real = np.array([], dtype=np.float64)
        self.imag = np.array([], dtype=np.float64)
        self.weight = np.array([], dtype=np.float64)
        self.flag = np.array([], dtype=np.float64)
        self.filename = []

        self.amp = np.array([], dtype=np.float64)
        self.phase = np.array([], dtype=np.float64)

    def readHDF5(self, filename=None):
        if os.path.isfile(filename) is False:
            raise ValueError('File does not exist')
        if filename is not None:
            fid = h5py.File(filename, "r")
            freqs = fid["freqs"][:]     # [Hz]
            uu = fid["uu"][:,:]         # [kilolam]
            vv = fid["vv"][:,:]         # [kilolam]
            real = fid["real"][:,:]     # [Jy]
            imag = fid["imag"][:,:]     # [Jy]
            weight = fid["weight"][:,:]         #[1/Jy^2]
            flag = fid["flag"][:,:]             # boolean mask
            attributes = fid.attrs
            fid.close()

        self.uu = np.squeeze(uu)
        self.vv = np.squeeze(vv)
        self.real = np.squeeze(real)
        self.imag = np.squeeze(imag)
        self.weight = np.squeeze(weight)
        self.flag = np.squeeze(flag)
        self.freqs = freqs
        self.filename = [filename]

    def writeHDF5(self, filename=None, objname=None):
        nfreq = len(self.freqs)

#        shape = (nfreq, nvis)
        # technically, need to be nfreq by nvis, 
        # but right now, collapse all to nfreq=1, until I think of
        # how to deal with flagging and still keep the 2d array
        nvis = len(self.uu)
        nfreq = 1
        shape = (nfreq, nvis)

        if filename is None:
            filename = "data.hdf5"
        if objname is None:
            objname = "My target"

        fid = h5py.File(filename, "w")

        fid.attrs["OBJECT"] = objname  # attributes are added like dictionary values in Python

        fid.create_dataset("freqs", (nfreq,), dtype="float64")[:] = self.freqs # [Hz]

        fid.create_dataset("uu", shape, dtype="float64")[:,:] = self.uu # [kilolambda]
        fid.create_dataset("vv", shape, dtype="float64")[:,:] = self.vv # [kilolambda]

        fid.create_dataset("real", shape, dtype="float64")[:,:] = self.real # [Jy]
        fid.create_dataset("imag", shape, dtype="float64")[:,:] = self.imag # [Jy]

        fid.create_dataset("weight", shape, dtype="float64")[:,:] = self.weight #[1/Jy^2]
        fid.create_dataset("flag", shape, dtype="bool")[:,:] = self.flag # Boolean array
        fid.close()

    def getAmpPhase(self):
        self.amp = np.sqrt(self.real**2. + self.imag**2.)
        self.phase = np.arctan2(self.imag, self.real)

    def applyflag(self, flagreg=None):
        """ 
        apply the flag. keep the good ones, and just ignore the bad ones
        flagreg : intarray, in size of uu,vv, etc
        """
        if flagreg is None:
            flagreg = self.flag == 1
        tokeep = ~flagreg

        self.uu = self.uu[tokeep]
        self.vv = self.vv[tokeep]
        self.real = self.real[tokeep]
        self.imag = self.imag[tokeep]
        self.weight = self.weight[tokeep]
        self.flag = self.flag[tokeep]

    def river2UVTable(self):
        """
        output this data into what uvplot wants
        """
        uvtab = uvplot.UVTable()
        uvtab._u = self.uu * 1e3 # uvtable is in lambda
        uvtab._v = self.vv * 1e3
        uvtab._re = self.real
        uvtab._im = self.imag
        uvtab._weights = self.weight
        uvtab.ndat = len(self.uu)
        return uvtab

    def plotUVhisto2d(self, ax=None):
        """
        plot the scattered data into 2d histograms, or density plots
        """
        if ax is None:
            ax = plt.gca()

        uvdist = np.sqrt(self.uu**2 + self.vv**2)
        ax.hist2d(uvdist, self.amp, 
data[:, 1], bins=100)
