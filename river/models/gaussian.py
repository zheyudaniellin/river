def getDefaultParams():
    """
    Provide default parameters
    """
    defpar = [
        ['amp', '1.', 'amplitude'],
        ['uv0', '[0,0]', 'phase center'],
        ['sig', '[50, 50]', 'sigma in lambda']
             ]
    return defpar

def getVis(uu=None, vv=None, ppar=None):
    """
    u = 1d array
    v = 1d array
    par = list
        parameters
    """
    nu = len(uu)
    nv = len(vv)
    mesh = np.meshgrid(uu, vv, indexing='ij')
    umesh = mesh[0]
    vmesh = mesh[1]

    usig = ppar['sig'][0]
    vsig = ppar['sig'][1]

    vis = (ppar['amp'] * np.exp(-((umesh-ppar['uv0'][0])/usig)**2.) * 
            np.exp(-((vmesh-ppar['uv0'][1])/vsig)**2.))
    return vis
