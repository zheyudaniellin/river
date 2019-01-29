# emceeTool.py
# emceeTool object as wrapper for calling and using emcee
import emcee
import params
import radmc3dTool
import galarioTool
import time
import numpy as np
import pdb
from radmc3dPy import image

class emceeTool(object):
    """
    Attributes
    ----------
    nwalkers	: int
                  number of walkers
    nthreads	: int. Default 1
                  number of threads to use emcee
    steps	: int
                  number of steps for emcee
    waterpar	: WaterPar object
                  the waterpar that will be continuualy updated by emcee. must be initiated

    """

    def __init__(self):
        self.nwalkers = []
        self.nthreads = 1
        self.steps = []
        self.pball = []
        self.waterpar = []
        self.riverpar = []
        self.sampler = []

    def getPball(self):
        # a whole lot of ways to initialize guess, but pball is a good way
        parCin = np.array(self.waterpar.parC, dtype=np.float64).copy()
        parWin = self.waterpar.parW
        pball = emcee.utils.sample_ball(parCin, parWin, size=self.nwalkers)
        self.pball = pball

    def doSampler(self, obsgrid):
        nparC = len(self.waterpar.parC)
        sampler = emcee.EnsembleSampler(self.nwalkers, nparC, lnprob,
                      args=[self.waterpar, self.riverpar, obsgrid])

        starttime = time.clock()
        sampler.run_mcmc(self.pball, self.steps)
        stoptime = time.clock()

        print('Sampler time: %fs'%(stoptime-starttime))
        self.sampler = sampler 
        
def lnprior(waterpar):
    """
    calculate prior probability function. Most simple form: uniform
    Parameters
    ----------
    waterpar : WaterPar object                   
               must be already initiated
    """
    fail_prior = 0
    nparC = len(waterpar.parC)
    for ipar in range(nparC):
        if waterpar.parC[ipar] < waterpar.parL[ipar]: fail_prior = 1
        if waterpar.parC[ipar] > waterpar.parU[ipar]: fail_prior = 1
    return fail_prior
   
def lnprob(par, waterpar, riverpar, obsgrid):
    """
    Calculates the probability function, called by emcee.EnsembleSampler 
    """
    nparC = len(waterpar.parC)
    for ipar in range(nparC):
        waterpar.updateParC(waterpar.parCL[ipar], par[ipar])

    fail_prior = lnprior(waterpar)
    if fail_prior:
        return -np.inf
    else:
        chi, cup, im, net = getChi(waterpar, riverpar, obsgrid)
        return -chi/2.0

def getChi(waterpar, riverpar, obsgrid):
    """
    calculate the chi value
    """
    # produce model and image
    cup = radmc3dTool.makeModel(waterpar, riverpar)

    # read the image
    im = image.readImage(fname=cup.imfname, dpc=riverpar.dis)

    # produce the visibilities and Chi
    net = galarioTool.makeVis(waterpar, riverpar, im, obsgrid)

    rchi2 = net.rchi2
    chiout = np.sum(rchi2)
    return chiout, cup, im, net

