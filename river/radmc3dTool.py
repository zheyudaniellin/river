# radmc3dTool.py
# used to set up radmc3d parameters and files
# model specific parameters are in WaterPar()
# this will be general code executions, should not enter emcee calculations
import numpy as np
from radmc3dPy import *
import os
import pdb

class radmc3dTool(object):
    """ 
    Attributes
    ----------
    modelname   : string,
                  name of the radmc3dPy model
    inppar      : dictionary
                  used for input for radmc3d
    binary      : boolean, default False
                  True for binary output files
    domctherm   : boolean, default False
                  True to use radmc3d to calculate temperature
    impar  	: dictionary
                  inputs for making image
    radmc3dPar	: radmc3dPar object
                  the total parameters for the radmc3d model
    radmc3dGrid	: radmc3dGrid object
                  the grid object 
    rundir	: string
                  the local path name of run
    projdir 	: string
                  the global location of working directory
    imfname	: string
                  file name of radmc3dImage, with rundir in front
    """
    def __init__(self, rundir=None, projdir=None): 
        self.modelname = 'simple_1'
        self.opacpar = {}
        self.opactype = 'Beck'
        self.binary=False
        self.domctherm = False
        self.impar = {
            'npix':300,
            'phi':0,
            'sizeau':200.,
            'stokes':False,
            'secondorder':0,
            'nostar':False,
            'camwav':[870.] #in microns
            }
        self.radmc3dPar = []
        self.radmc3dGrid = []
        self.rundir = rundir
        self.projdir = projdir
        self.imfname = []
        

    def getDefaultPar(self, fname=None):
        """
        Load default parameter file and write into current directory
        Or given fname, read the file and write into current directory
        """
        par = params.radmc3dPar()
        if fname is None:
            par.loadDefaults(model=self.modelname)
        else:
            if os.path.isfile(fname):
                par.readPar(fname=fname)
            else:
                raise ValueError('Cannot find parameter file: %s'%(fname))

        self.radmc3dPar = par
        self.radmc3dPar.writeParfile(fdir=self.rundir)

    @staticmethod 
    def ppar2string(pparin):
        pparout = pparin.copy()
        for ikey in pparin.keys():
            if isinstance(pparin[ikey], float):
                pparout[ikey] = ("%.7e"%pparin[ikey])
            elif isinstance(pparin[ikey], int):
                pparout[ikey] = ("%d" %pparin[ikey])
            elif isinstance(pparin[ikey], str):
                pparout[ikey] = pparin[ikey]
            elif isinstance(pparin[ikey], list):
                dum = '['
                for ii in range(len(pparin[ikey])):
                    if type(pparin[ikey]) is float:
                        dum = dum + (".7e" % pparin[ikey][ii])
                    elif type(pparin[ikey]) is int:
                        dum = dum + ("%d" % pparin[ikey][ii])
                    elif type(pparin[ikey]) is str:
                        dum = dum + (pparin[ikey][ii])
                    else:
                        raise TypeError('Unknown data type in '+ikey)
                    if ii < len(pparin[ikey])-1:
                        dum = dum + ', '
                dum = dum + ']'
                pparout[ikey] = dum
        return pparout

    def updatemcPar(self, ppar=None):
        """ updates the radmc3dPar and its file by ppar dictionary
            can use waterpar.ppar for waterpar
        """
        pparstr = self.ppar2string(ppar)
        kkeys = pparstr.keys()
        for ikey in kkeys:
#            self.radmc3dPar.ppar[ikey] = ppar[ikey] #not sure if this is necessary   
            self.radmc3dPar.setPar([ikey, pparstr[ikey], '', ''])
        self.radmc3dPar.writeParfile(fdir=self.rundir)

    def getOpacPar(self, specs=None):
        """ opacpar is a bit complicated, so add this method
        Parameters
        ----------
        specs : list, in string
                Names of the refractive index you want. listed in set_dustspec
        """
        if specs is None:
            specs = ['Sil_Draine']
        dustspec = zylutils.set_dustspec.dustspec()
        dustspec.getSpec(specs=specs)
        mixabun = np.array(dustspec.mabun) / sum(dustspec.mabun)

        self.opacpar = {
            'nw':[25, 50, 50], 'wbound':[0.1, 7.0, 25., 1e4],
            'lnk_fname':dustspec.lnkfiles,
            'nscatang':361, 'chopforward':0,
            'logawidth':0.05, 'na':20, 
            'mixgsize':0,
            'gdens':dustspec.swgt, 'mixabun':mixabun,
            'miescat_verbose':False,
            'wfact':3.0, 'extrapolate':True, 'errtol':0.01,
            'kpara0':0.05, 'ksca0':1., 
            'beta':1.0, 'beck0':10. #gsmin, gsmax, etc is already in default parfile
            }
        kkey = self.opacpar.keys()
        for valname in kkey:
            self.radmc3dPar.setPar(parlist=
                [valname, str(self.opacpar[valname]), '', 'Dust opacity'])
        self.radmc3dPar.writeParfile(fdir=self.rundir)

    def getmcOpac(self, waterpar=None):
        if waterpar is None:
            raise ValueError('Need waterpar obj')
        nw = self.opacpar['nw']
        wbound = self.opacpar['wbound']
        radpar = {'nw':nw, 'wbound':wbound}
        grid = reggrid.radmc3dGrid()
        grid.makeWavelengthGrid(ppar=radpar)
        op = dustopac.radmc3dDustOpac()

        if self.opactype is 'Beck':
        # needs beta and beck0
            ppin = {}
            ppin.update(self.opacpar)
            if waterpar.isLabelinPar('beta'):
                beta = waterpar.getParString('beta')
            else:
                beta = ppin['beta']
            if waterpar.isLabelinPar('beck0'):
                beck0 = waterpar.getParString('beck0')
            else:
                beck0 = ppin['beck0']
            op.makeBeckOpac(wav=grid.wav, beck0=beck0, beta=beta, fdir=self.rundir)

        elif self.opactype is 'Mie':
        # needs gsmin, gsmax, ngs, gsdist_powex
            miearg = ['gsmin', 'gsmax', 'ngs', 'gsdist_powex']
            nargs = len(miearg)
            ppin = {}
            ppin.update(self.radmc3dPar.ppar)
            ppin.update(self.opacpar)
            for ii in range(nargs):
                if waterpar.isLabelinPar(miearg[ii]):
                    ival = waterpar.getParString(miearg[ii])
                    ppin.update({miearg[ii]:ival})
            ksca0 = self.opacpar['ksca0']
            op.makeOpac(ppar=ppin, ksca0=ksca0, fdir=self.rundir)

        mopac = op.readMasterOpac(fdir=self.rundir)
        op.readOpac(ext=mopac['ext'], scatmat=mopac['scatmat'][0], fdir=self.rundir)

        return op

    def getmcModel(self, waterpar=None):
        ppin = self.radmc3dPar.ppar
        if waterpar is not None:
            ppin.update(waterpar.ppar)

        model = setup.radmc3dModel(binary=self.binary, model=self.modelname)
        model.readParams(fdir=self.rundir)
        model.writeRadmc3dInp(fdir=self.rundir)
        model.makeGrid(writeToFile=True, fdir=self.rundir)
        model.makeRadSources(writeToFile=True, fdir=self.rundir)
        model.makeVar(ddens=True, writeToFile=True, fdir=self.rundir)
        if self.domctherm is False:
            model.makeVar(dtemp=True, writeToFile=True, fdir=self.rundir)
        return model

    def getmcImage(self, waterpar=None, fname=None):
        if fname is None:
            fname = 'myimage.out'

        inc = waterpar.getParString('inc')
        image.writeCameraWavelength(camwav=self.impar['camwav'], fdir=self.rundir)

        argin = {} #give another dictionary for the actual input, just for easier coding
        argin.update(self.impar)
        os.chdir(self.rundir)
        image.makeImage(npix=argin['npix'], sizeau=argin['sizeau'],
                        loadlambda=1, incl=inc, phi=argin['phi'], 
                        stokes=argin['stokes'], nostar=argin['nostar'],
                        secondorder=argin['secondorder'], fname=fname)
        os.chdir(self.projdir)

        # store the rundir and image name for convenience
        if self.rundir[-1] is '/':
            self.imfname = self.rundir + fname
        else:
            self.imfname = self.rundir + '/' + fname

        # just check if the file exist
        if os.path.isfile(self.imfname) is False:
            raise ValueError('Something wrong with image file name')

def makeModel(waterpar, riverpar):
    """
    get a single image based on parameters in waterpar, riverpar
    """
    cup = radmc3dTool(rundir=riverpar.rundir, projdir=riverpar.projdir)
    cup.modelname = riverpar.modelname
    cup.getDefaultPar()
    cup.updatemcPar(ppar=riverpar.radpar)
    cup.updatemcPar(ppar=waterpar.ppar)
    cup.getOpacPar(specs=riverpar.specs)
    cup.opactype=riverpar.opactype
    op = cup.getmcOpac(waterpar=waterpar)
    model = cup.getmcModel()
    cup.getmcImage(waterpar=waterpar)
    return cup
