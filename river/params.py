# params.py
# used to handle parameters
import numpy as np
import pdb

class WaterPar(object):
    """
    Used to split up arguments based on lower, upper limits, parameter width, units, etc
    """
    def __init__(self):
        """
        Initialize object
        Attributes
        ----------
        par 	: dictionary
              	in order of 
                  turn on/off, value, sigma, lower limit, upper limit, unit

        parCL	: list, in string
                label of parameter to mcmc

        parC 	: list, in float
               	value of parameter for mcmc

        par0	: list, in float
                initial value of parameter for mcmc in par

        parCunit : list, in string
                unit of parameter for mcmc, and for sigma and limits

        parW 	: list, in float
               	sigma of parameter

        parL 	: list, in float
               	lower limit of parameter

        parU 	: list, in float
               	upper limit of parameter

        parSL 	: list, in string
               	label of fixed parameter

        parS	: list, in float
               	value of fixed parameter              

        parSunit : list, in string
                units of the fixed parameter

        listpar	: list, in string
               	name of arguments that are a list

        ppar	: dictionary 
               	the actual dictionary that can go into radmc3d
        """
        self.par = {}
        self.parCL = []
        self.parC = []
        self.par0 = []
        self.parCunit = []
        self.parW = []
        self.parL = []
        self.parU = []
        self.parSL = []
        self.parS = []
        self.parSunit = []
        self.listpar = []
        self.ppar = {}

    def initiatePar(self):
        """
        Determine which of the parameters in par are turned on/off
        """
        kkeys = self.par.keys()
        # determine which ones are lists
        listpar = []
        for ik in kkeys:
            if type(self.par[ik][0]) is list:
                listpar.append(ik)
        self.listpar = listpar

        for ik in kkeys:
            # parameters that are in a list
            if ik in self.listpar:
                nelem = len(self.par[ik])
                for ii in range(nelem):
                    if self.par[ik][ii][0] is 0: # turned off
                        self.parSL.append('%s%d'%(ik,ii))
                        self.parS.append(self.par[ik][ii][1])
                        self.parSunit.append(self.par[ik][ii][5])
                    elif self.par[ik][ii][0] is 1: # turned on
                        self.parCL.append('%s%d'%(ik,ii))
                        self.parC.append(self.par[ik][ii][1])
                        self.parW.append(self.par[ik][ii][2])
                        self.parL.append(self.par[ik][ii][3])
                        self.parU.append(self.par[ik][ii][4])
                        self.parCunit.append(self.par[ik][ii][5])
                    else:
                        raise ValueError('error in turn on/off')
            # parameters that are in a single value
            else:
                if self.par[ik][0] is 0: # turned off
                    self.parSL.append(ik)
                    self.parS.append(self.par[ik][1])
                    self.parSunit.append(self.par[ik][5])
                elif self.par[ik][0] is 1: # turned on
                    self.parCL.append(ik)
                    self.parC.append(self.par[ik][1])
                    self.parW.append(self.par[ik][2])
                    self.parL.append(self.par[ik][3])
                    self.parU.append(self.par[ik][4])
                    self.parCunit.append(self.par[ik][5])
                else:
                    raise ValueError('error in turn on/off')
        dumpar0 = np.array(self.parC, dtype=np.float64).copy()
        self.par0 = dumpar0.tolist()
        self.ppar = self.getPpar()

        self.checkPar()

    def checkPar(self):
        for ix in range(len(self.parC)):
            if self.parL[ix] > self.parU[ix]:
                print('Name: %s'%self.parCL[ix])
                print('Lower: %.2e'%self.parL[ix])
                print('Upper: %.2e'%self.parU[ix])
                raise ValueError('Error in upper and lower limits for %s'%self.parCL[ix])

    def updateParC(self, name, val):
        """
        Update value of parC and its correspondence in ppar
        Should not be updating the fixed parameters: parS, parSL
        Parameters
        ----------
        name : string
               the name of parameter
        val  : float
        """
        if name in self.parCL:
            # update parC
            inx = self.parCL.index(name)
            self.parC[inx] = val
        elif name in self.parSL:
            raise ValueError('Should not update the fixed parameters')
        else:
            raise ValueError('Cannot find a name to update Par')

        # update ppar
        self.ppar = self.getPpar()

    def getPpar(self):
        """
        ppar is a dictionary especially for radm3dPy inputs
        where its values are strings, even for lists
        assumes that parC, parS, etc are all set
        """
        ppar = {}
        kkeys = self.par.keys()
        for ikey in kkeys:
            if ikey in self.listpar:
                strval = self.getParStringList(ikey)
            else:
                strval = self.getParString(ikey)
            ppar.update({ikey:strval})
        return ppar

    def getParString(self, name): 
        if name in self.parCL:
            inx = self.parCL.index(name)
            elem = self.parC[inx]
            unit = self.parCunit[inx]
        elif name in self.parSL:
            inx = self.parSL.index(name)
            elem = self.parS[inx]
            unit = self.parSunit[inx]
        else:
            raise ValueError(name+' not in parameters')
        if unit is '':
            spar = '%e'%(elem)
        else:
            spar = '%e*%s'%(elem, unit)
        return spar

    def getParStringList(self, pref):
        pdir = {}
        for ip in range(len(self.parCL)):
            if pref in self.parCL[ip]:
                pdir.update({self.parCL[ip]:[self.parC[ip], self.parCunit[ip]]})
        for ip in range(len(self.parSL)):
            if pref in self.parSL[ip]:
                pdir.update({self.parSL[ip]:[self.parS[ip], self.parSunit[ip]]})
        pkeys = pdir.keys()
        plist = [] # list of parameters in order of names
        for ip in range(len(pkeys)):
            name = '%s%d'%(pref, ip)
            if name in pkeys:
                plist.append(pdir[name])
            else:
                raise ValueError('some problem in transferring list of string for radmc3dpy')
        spar = '['
        for ip in range(len(plist)):
            unit = plist[ip][1]
            if unit is '':
                apex = ''
            else:
                apex = '*'+unit
            spar = spar + '%e'%(plist[ip][0]) + apex
            if ip != (len(plist)-1): 
                spar = spar + ','
        spar = spar + ']'

        return spar

    def getValbyName(self, name):
        # get a non-string value for numbers. usually for single parameter
        if name in self.parCL:
            inx = self.parCL.index(name)
            val = self.parC[inx]
        elif name in self.parSL:
            inx = self.parSL.index(name)
            val = self.parS[inx]
        else:
            raise ValueError('value by name not found')
        return val

    def getValbyIndex(self, pref, inx):
        # all the parC, etc are set. returns the value
        if pref not in self.listpar:
            raise ValueError('Parameter name does not exist')
        pname = pref + '%d'%(inx)
        if pname in self.parCL:
            ii = self.parCL.index(pname)
            val = self.parC[ii]
        elif pname in self.parSL:
            ii = self.parSL.index(pname)
            val = self.parS[ii]
        else:
            raise ValueError('value by index not found')
        return val

    def isLabelinPar(self, name):
        """ test if a label name is in waterpar
        Return: boolean
        """
        sig = False
        if name in self.par.keys():
            sig = True
        return sig

    def writeWaterPar(self, fdir=None):
        """ write to a file
        """
        fname = fdir + 'water.par'

    def readWaterPar(self, fdir=None):
        """ read into a waterpar, then can start to use initiatepar()
        """
        fname = fdir + 'water.par'

class RiverPar(object):
    """
    object to handle parameters for running the whole flow of code
    """
    def __init__(self, rundir=None, projdir=None):
        """ initial attributes
        """
        self.par = [] #input for waterpar
        self.rundir = rundir
        self.projdir = projdir
        self.opactype = []
        self.specs = []
        self.dis = []
        self.radpar = []
        self.modelname= []

    def getDefaultRadpar(self):
        radpar = {
            'scattering_mode_max':0, 'alignment_mode':0,
            'setthreads':16,
            'incl_lines':0, 'incl_dust':1,
            'nphot':1e6, 'nphot_scat':1e6, 'modified_random_walk':1,
            'mc_scat_maxtauabs':15
             }
        self.radpar = radpar

