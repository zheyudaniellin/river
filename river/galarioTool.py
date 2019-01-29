# galarioTool.py
# tool to call galario
import numpy as np
import galario
import UVdata
from radmc3dPy import natconst
import pdb

class galarioTool(object):
    """
    Attributes
    ----------
    rundir : string
              the directory where all the input/output files are
    projdir : string
              directory where the command file is

    radmc3dimage : radmc3d image object
                   only single object, though can be multi wavelength and stokes

    vis : 2d list
          the calculated visibilities of the image, in shape of stokes (I,Q,U,V) by 
          wavelength (in order of camera_wavelength grid)
          Ex: stokes I, at 2nd wavelength = vis[0][1]. Stokes U at 3rd wavelength = vis[2][2]
            
    uu  : 1d list
          u in klambda for each wavelength

    vv  : 1d list 
          v in klambda for each wavelength

    PA  : float in degrees

    dRA : 1d list in arcseconds
          offset for each wavelength

    dDec: 1d list in arcseconds
          offset in DEC for each wavelength

    chi2 : list
          chi squared in dimensions like vis. not normalized

    rchi2 : list
            reduced chi squared in dimensions like vis. 

    """
    def __init__(self, rundir=None, projdir=None):
        self.rundir = rundir
        self.projdir = projdir
        self.radmc3dImage = []
        self.vis = []
        self.uu = []	# [klambda]
        self.vv = []	# [klambda]
        self.dRA = [] 	# [arcsec]
        self.dDec = []	# [arcsec]
        self.PA = []	# [deg]
        self.chi2 = []	# not normalized
        self.rchi2 = []

    def getVis(self):
        # must have self.uu, self.vv ready
        # set up the dimensions
        nwav = len(self.radmc3dImage.wav)
        stokes = self.radmc3dImage.stokes
        if stokes:
            nstokes = 4
        else:
            nstokes = 1

        # calculate the vis
        vis = range(nstokes)
        for istokes in range(nstokes):
            vis[istokes] = range(nwav)
            for iwav in range(nwav):
                if stokes:
                    imii = self.radmc3dImage.imageJyppix[:,:,istokes, iwav].copy(order='C')
                    imii = np.squeeze(imii)
                else:
                    imii = self.radmc3dImage.imageJyppix[:,:,iwav].copy(order='C')
                    imii = np.squeeze(imii)
                # prepare the inputs and conversions
                uii = self.uu[iwav] * 1e3
                vii = self.vv[iwav] * 1e3
                dxy = self.radmc3dImage.sizepix_x / (self.radmc3dImage.dpc * natconst.pc)
                draii = self.dRA[iwav] / 3600. * np.pi / 180. #need [rad]
                ddecii = self.dDec[iwav] / 3600. * np.pi / 180. # need [rad]
                paii = self.PA * np.pi / 180. # input pa in [rad]
                visii = galario.double.sampleImage(imii, dxy, uii, vii,
                            PA=paii, dRA=draii, dDec=ddecii)
                vis[istokes][iwav] = visii
        self.vis = vis

    def getVisChi2(self, obsgrid):
        """ 
        organizes input to calculate chi2
        obsgrid : list, in UVdata object
                    should be in same shape as self.vis. nstokes by nwav
        """
        nstokes = len(self.vis)
        nwav = len(self.vis[0])
        chi2 = range(nstokes)
        rchi2 = range(nstokes)
        for istokes in range(nstokes):
            chi2[istokes] = range(nwav)
            rchi2[istokes] = range(nwav)
            for iwav in range(nwav):
                obsii = obsgrid[istokes][iwav]
                modii = self.vis[istokes][iwav]
                uii = self.uu[iwav]
                vii = self.vv[iwav]
                cii, rcii = self.getChi2(uii, vii, modii, obsii)

                chi2[istokes][iwav] = cii
                rchi2[istokes][iwav] = rcii
        self.chi2 = chi2
        self.rchi2 = rchi2

    def toUVdata(self, istokes=None, iwav=None):
        """
        Output all the visibilities to UVdata object
        Parameters
        ----------
        istokes : int, optional. Default all
                  the desired stokes: I, Q, U, V
        iwav    : int, optional. Default all
                  the desired wave
        Returns
        -------
        data  : list, in UVdata objects
                    or 
                    single UVdata object. if input istokes and iwav
        """
        # if no index given, output all 
        if (istokes is None) and (iwav is None):
            nwav = len(self.radmc3dImage.wav)
            stokes = self.radmc3dImage.stokes
            if stokes:
                nstokes = 4
            else:
                nstokes = 1
            data = range(nstokes)
            for istokes in range(nstokes):
                data[istokes] = range(nwav)
                for iwav in range(nwav):
                    datii = UVdata.UVdata()
                    visii = self.vis[istokes][iwav]
                    uii = self.uu[iwav]
                    vii = self.vv[iwav]
                    datii.real = np.real(visii)
                    datii.imag = np.imag(visii)
                    datii.uu = uii
                    datii.vv = vii
                    datii.weight = np.ones(visii.shape, dtype=int)
                    datii.flag = np.zeros(visii.shape, dtype=int)
                    datii.freqs = np.array([self.radmc3dImage.freq[iwav]])
                    data[istokes][iwav] = datii
        else:
            data = UVdata.UVdata()
            vis = self.vis[istokes][iwav].copy()
            data.uu = self.uu[iwav].copy()
            data.vv = self.vv[iwav].copy()
            data.real = np.real(vis)
            data.imag = np.imag(vis)
            data.weight = np.ones(vis.shape)
            data.flag = np.zeros(vis.shape, dtype=int)
            data.freqs = np.array([self.radmc3dImage.freq[0]])
        return data 

    @staticmethod
    def getChi2(uu, vv, vis, obsvis):
        """
        Calculate chi squared and reduced chi squared by 
            chi**2 = sum(weight* ((real(obs) - real(model))**2 + (imag(real)-imag(obs))**2))

        Parameters
        ----------
        uu : ndarray
             the u in klambda 
        vv : ndarra
             the v in klambda
        vis : ndarray, complex
              the visibilities of model data
        uvobs : UVdata object
                UVdata of the observed. its arrays should be in same number and dimensions
                as vis, and uu,vv should be the same as uvobs' uu, vv
        """
        # do some checking
        if isinstance(obsvis, UVdata.UVdata) is False:
            raise ValueError('observed vis input has to be an instance of UVdata')

        umodshape = uu.shape
        vmodshape = vv.shape
        if umodshape != vmodshape:
            raise ValueError('input uu is not in the same shape as vv')

        uobsshape = obsvis.uu.shape
        if umodshape != uobsshape:
            raise ValueError('input uu is not in the same shape as observation uu')

        if obsvis.weight.max() == 0:
            raise ValueError('maximum data point weight is 0')

        # do calculation
        nvis = len(vis)
        realmod = np.real(vis)
        imagmod = np.imag(vis)
        chi2 = np.sum(obsvis.weight * ((obsvis.real - realmod)**2. + (obsvis.imag - imagmod)**2.))
        rchi2 = chi2 / nvis

        return chi2, rchi2

def makeVis(waterpar, riverpar, im, obsgrid):
    """
    produce the visibilities based on radmc3d Image
    Parameter
    ---------
    waterpar : 
    riverpar : 
    im : radmc3dImage
    obsgrid : list of UVdata
              in shape of nstokes x nwav. 
    """
    net = galarioTool(rundir=riverpar.rundir, projdir=riverpar.projdir)
    net.radmc3dImage = im
    net.PA = waterpar.getValbyName('PA')
    net.uu = range(im.nwav)
    net.vv = range(im.nwav)
    net.dRA = range(im.nwav)
    net.dDec = range(im.nwav)
    for iwav in range(im.nwav):
        obsii = obsgrid[0][iwav] # uu,vv should be stokes independent
        net.uu[iwav] = obsii.uu
        net.vv[iwav] = obsii.vv
        net.dRA[iwav] = waterpar.getValbyIndex('dRA', iwav)
        net.dDec[iwav] = waterpar.getValbyIndex('dDec', iwav)
    net.getVis()
    net.getVisChi2(obsgrid)
    return net
