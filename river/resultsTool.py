# resultsTool
"""
This module conatins the classes/function to plot the final results for inspection
- getTriangle = plot posterior probability
- getChiStep = plot values and chi for each step and every walker
- getBestModel = get the radmc3dImage, galarioTool, chi for the final best model
- getResidualImage = calculate the residual image 
- plot uvdist to amplitude, overlaying data to model 
"""
import corner
import numpy as np
import fntools
import pdb
import emceeTool
import matplotlib.pyplot as plt

class resultsTool(object):
    """
    used to calculate results
    """
    def __init__(self, riverpar=None):
        self.waterpar = []	# waterpar
        self.riverpar = riverpar# riverpar
        self.radmc3dTool = []	# radmc3dTool
        self.chain = []		# all the values
        self.samples = []	# the samples used for results, after burn-in
        self.fitval = []	# the best value and uncertainty

    def getSamples(self, fdir='', inxcut=None):
        """
        Parameters
        ----------
        inxcut : int
                 First index to consider in samples
        """
        if fdir is '':
            fname = 'dat.chain'
        else:
            if fdir[-1] is '/':
               fname = fdir+ 'dat.chain'
            else:
               fname = fdir + '/' + 'dat.chain'
        chain = fntools.zylreadvec(fname)
        self.chain = chain
        nwalkers, nsteps, ndim = chain.shape
        if inxcut is None:
            inxcut = int(np.floor(0.1 * nsteps))
        samples = chain[:,inxcut:,:].reshape((-1, ndim))
        self.samples = samples

    def getFitval(self, percentiles=[16,50,84]):
        """
        calculates the best values by percentiles
        Parameters
        ----------
        percentiles : 1 by 3 list, in percent
                      Default is [16,50,84]. Always keep the middle one 50
        """
        ntotsteps, ndim = self.samples.shape
        nval = len(percentiles)
        npercent = len(percentiles)
        fitval = np.zeros([ndim, npercent], dtype=np.float64)
        for ii in range(ndim):
            valsii = np.percentile(self.samples[:,ii], percentiles)
            fitval[ii,:] = valsii
        self.fitval = fitval

    def getTriangle(self, fname=None, **kwargs):
        """
        plots the triangle plot. 
        fname : string
                the full name of plot to be saved. Default is 'triangle.png' in rundir
        **kwargs : arguments to corner
        
        """
        if fname is None:
            fname = self.riverpar.rundir + 'triangle.png'
        fig = corner.corner(self.samples,
            labels=self.waterpar.parCL, truths=self.waterpar.par0, 
            show_titles=True, quantiles=[0.16, 0.50, 0.84], **kwargs)
        fig.savefig(fname)

        return fig

    def getChiStep(self, fname=None, **kwargs):
        """ 
        plots chi as a function of steps
        Parameters
        ----------
        fname : string
                the full name of plot to be saved. Default is 'chistep.png' in rundir
        **kwargs : arguments to plotting
        """
        if self.chain is []:
            raise ValueError('Chain must be read into resultsTool object first')
        if fname is None:
            fname = self.riverpar.rundir + 'chistep.png'
        nwalkers, nsteps, ndim = self.chain.shape

        nyplot = int(np.ceil(np.sqrt(ndim)))
        nxplot = nyplot

        fig, axes = plt.subplots(ncols=nxplot, nrows=nyplot, sharex='col')
        iax_y = 0
        for ic in range(ndim):
            iax_x = np.mod(ic, nxplot)
            iax_y = int(np.floor(ic / nyplot))
            for iw in range(nwalkers):
                axes[iax_x, iax_y].plot(self.chain[iw,:,ic])

            axes[iax_x, iax_y].set_ylabel(self.waterpar.parCL[ic])
        fig.tight_layout()
        fig.savefig(fname)
        return fig, axes

    def getBestModel(self, obsgrid):
        """
        use emceeTool.getChi to plot the best model
        """
        par = self.fitval[:,1]
        if par is []:
            raise ValueError('self.fitval must be obtained first. Run getFitval()')

        nparC = len(self.waterpar.parC)
        for ipar in range(nparC):
            self.waterpar.updateParC(self.waterpar.parCL[ipar], par[ipar])

        chi, cup, im, net = emceeTool.getChi(waterpar, riverpar, obsgrid)
        return chi, cup, im, net

    @staticmethod
    def plotUVamp(obsgrid, modelgrid, fdir='', uvbin_size=100e3):
        """
        plots the uvdistance to amplitude and phase. plot model to real data using UVTable
        Parameters
        ----------
        obsgrid : list of UVdata object

        modelgrid : list of UVdata object. can be obtained through galarioTool

        fdir    : string
                  directory to store the file
        uvbin_size : float
                     uv-distance bin in wavelength. argument for UVtable.plot()
        """
        stokestag = ['I', 'Q', 'U', 'V']
        freq = net.radmc3dImage.freq
        # translate to UVtable
        nstokes = len(obsgrid)
        if nstokes is 4:
            stokes = 1
        else:
            stokes = 0
        nwav = len(obsgrid[0])
        for istokes in range(nstokes):
            tabgrid[istokes] = range(nwav)
            for iwav in range(nwav):
                # convert model to UVtable
                modeltab = modelgrid[istokes][iwav].river2UVTable()

                # convert data to UVtable
                obstab = obsgrid[istokes][iwav].river2UVTable()

                # now plotting
                axes = obstab.plot(label='Data', uvbin_size=uvbin_size)
                modeltab.plot(label='Model', uvbin_size=uvbin_size, 
                    axes=axes, yerr=False, linestyle='-',color='r')
                
                if fdir is not '':
                    if fdir[-1] is not '/': fdir = fdir + '/'
                fname = fdir + 'uvamp_%s_f%d.png'%(stokestag[istokes],freq[iwav])
                axes[0].figure.savefig(fname)

    @staticmethod
    def getResidualImage(realim, im):
        """
        calculate the images residuals in pixel sizes, etc in realim
        Parameters
        ----------
        realim : radmc3dImage
                 The real image translated to radmc3dPy terms
        im     : radmc3dImage
                 The model image
        Return
        ------
        resim  : radmc3dImage
                 The residual image in radmc3dPy terms
        """


