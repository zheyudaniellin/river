# panelTool.py
""" collection of functions to assist in manipulating all of river
"""

import os
import matplotlib.pyplot as plt
from radmc3dPy import *
import numpy as np
import emcee
import galario
import fntools

class panelTool(object):
    """
    """
    def __init__(self):
        self.rundir = ''
        self.obsfile = ''
        self.riverpar = []
        self.waterpar = []

    @staticmethod
    def panel0(rundir, projdir, obsfile, par=None, riverpar=None):
        """ default panel
        """
        # make your rundir 
        os.system('rm -rf '+rundir)
        os.system('mkdir '+rundir)

        # read in observation
        obsuv = river.UVdata.UVdata()
        obsuv.readHDF5(filename=obsfile)

        obsgrid = [[obsuv]]

        # determine the image pixel size
        nxy, dxy = galario.double.get_image_size(obsuv.uu*1e3, obsuv.vv*1e3, verbose=True)

        # set up par
        if par is None:
            par = {
                'inc': [1, 45, 5, 0, 90, ''],
                'PA' : [1, 0, 5., -150, 150, ''],
                'dRA': [
                    [ 1, 0, 0.1, -0.5, 0.5, '']
                       ], 
                'dDEC': [
                    [ 1, 0, 0.1, -0.5, 0.5, '']
                        ]
                  }

         # set up riverpar
         if riverpar is None:
             riverpar = river.params.RiverPar(rundir=rundir, projdir=projdir)
             riverpar.getDefaultRadPar()
             river.dis = 1.
             riverpar.opactype='Beck'
             riverpar.specs = ['Sil_Draine']
             riverpar.modelname = 'disk_thin'

         # set up waterpar
         waterpar = river.params.WaterPar()
         waterpar.par = par
         waterpar.initiatePar()
         ndim = len(waterpar.parC)

         # do the emcee
         mill = river.emceeTool.emceeTool()
         mill.waterpar = waterpar
         mill.riverpar = riverpar
         mill.setps = 300
         mill.nwalkers = 8
         mill.nthreads = 16
         mill.getPball()
         mill.doSampler(obsgrid)
         sampler = mill.sampler

         return sampler

