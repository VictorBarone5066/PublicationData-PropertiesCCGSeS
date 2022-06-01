import numpy as np

pi = 3.141592654
sin = np.sin
cos = np.cos

#eV per cubic angstrom to GPa
def EvAngstToGPa(v):
    return v*160.2176487


#Gets the direction cosines for a vector defined by u*xhat + v*yhat + w*zhat
def GetDirCosines(u, v, w):
    dist = (u**(2) + v**(2) + w**(2))**(1./2.)
    #cos alpha: between vec(u, v, w) and x direction, cos beta, cos gamma
    return u/dist, v/dist, w/dist

#Make sure sij matrix is indexed such that s[1][1] = s11 instead of s[0][0] = s11.
#ref: Crystal Structure and Mechanical Properties of ThBC2 by Xinchun Zhou and Baobing Zheng
def GetDirYoungsMod(u, v, w, s):
    al, be, ga = GetDirCosines(u, v, w)
    return 1./(s[1][1]*al**(4) + s[2][2]*be**(4) + s[3][3]*ga**(4) + \
               (2*s[1][2] + s[6][6])*al**(2)*be**(2) + \
               (2*s[2][3] + s[4][4])*be**(2)*ga**(2) + \
               (2*s[1][3] + s[5][5])*al**(2)*ga**(2))

#Gives a set of x, y, z points that sample the surface of a unit sphere (i.e. r^2 = x^2 + y^2 + z^2)
def GetCircleSampling(normal, rad=1., ptDen=100):
    x, y, z = [], [], []
    #normal to x
    if(normal == 'x' or normal == 'u' or normal == 'a'):
        for phi in list(np.linspace(0, 2*pi, ptDen)):
            y.append(rad*cos(phi))
            z.append(rad*sin(phi))
        return None, y, z
    #normal to y
    if(normal == 'y' or normal == 'v' or normal == 'b'):
        for phi in list(np.linspace(pi, 3*pi, ptDen)):
            x.append(rad*cos(phi))
            z.append(rad*sin(phi))
        return x, None, z
    #normal to z
    if(normal == 'z' or normal == 'w' or normal == 'c'):
        for phi in list(np.linspace(0, 2*pi, ptDen)):
            x.append(rad*cos(phi))
            y.append(rad*sin(phi))
        return x, y, None
    return None, None, None
