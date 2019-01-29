# setup_radmc3d.py
import matplotlib.pyplot as plt
from radmc3dPy import *
import river

rundir = 'testing1/'
projdir = '/scratch/zdl3gk/coding/python/river_code/river_v2/projects/blueprint/'

par = { #on/off, value, sigma, lower, upper, unit
    'rho0': [ 1, 1e-15, 5e-16, 3e-16, 5e-15, ''],
    'T0':   [ 1, 101., 10., 80., 120., ''], 
    'radius': [ 0, 25., 5., 15, 30., 'au'], 
    'inc':  [ 0, 0., 0., 0., 90., '']
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
cup.getDefaultPar()
cup.updatemcPar(ppar=radpar)
cup.updatemcPar(ppar=waterpar.ppar)
cup.getOpacPar(specs=['Sil_Draine'])
cup.opactype='Beck'
op = cup.getmcOpac(waterpar=waterpar)
model = cup.getmcModel()
cup.getmcImage(waterpar=waterpar)

