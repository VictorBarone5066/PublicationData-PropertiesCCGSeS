import headerBandPlot as h

from copy import deepcopy
from scipy.optimize import curve_fit as cf
import numpy as np
from scipy.interpolate import CubicSpline

import numpy as np
VBM_MIN_OCC = 0.99 #if VBM_MIN_OCC <= occupancy <= 1.0, occupied by an electron
CBM_MAX_OCC = 0.01 #if 0.0 <= occupancy <= CBM_MAX_OCC, occupied by a hole

ALMOST_EQUAL_TOL = 1E-6

HBAR = 1.054571817E-34 #Js
m0 = 9.1093837015E-31 #kg

def EvToJo(ev):
    return ev*1.602176565e-19

def fmt(s):
    return f'{s:.3f}'

#Get highest energy occupied band info of a k-point.  Returns band instance
def GetHighestElectronBand(bandList):
    #search up from lowest energy
    for i in range(0, len(bandList)):
        if((VBM_MIN_OCC <= bandList[i][2]) and (bandList[i+1][2] < VBM_MIN_OCC)):
            return bandList[i]
    return None ##uh oh

#Get lowest energy unoccupied band info of a k-point.  Returns band instance
def GetLowestHoleBand(bandList):
    #search up from lowest energy since python dosen't optimize while loops
    for i in range(0, len(bandList)):
        if((bandList[i][2] <= CBM_MAX_OCC) and (CBM_MAX_OCC < bandList[i-1][2])):
            return bandList[i]
    return None ##uh oh

#Returns the n'th band from a list of bands
##TODO: Make this not awfully inefficient.  This works for now because I'm scared of 0-index errors
def GetNthBand(n, bandList):
    for band in bandList:
        if(band[0] == n):
            return band

#Get the k-point closest to the top of the valence band
def GetVBMKpoint(kpLis):
    bestKp, bandMaxEnergy = kpLis[0], kpLis[0].bands[0][1] ##initialize maxEnergy to the bottom of the vb
    for kp in kpLis:
        if(GetHighestElectronBand(kp.bands)[1] > bandMaxEnergy):
            bestKp = kp
            bandMaxEnergy = GetHighestElectronBand(kp.bands)[1]
    return bestKp

#Get the k-point closest to the bottom of the conduction band
def GetCBMKpoint(kpLis):
    bestKp, bandMinEnergy = kpLis[-1], kpLis[-1].bands[-1][1] ##initialize minEnergy to the top of the cb
    for kp in kpLis:
        if(GetLowestHoleBand(kp.bands)[1] < bandMinEnergy):
            bestKp = kp
            bandMinEnergy = GetLowestHoleBand(kp.bands)[1]
    return bestKp

#Feed me lists
def VecsAreEqual(a, b):
    return ((abs(b[0]-a[0]) < ALMOST_EQUAL_TOL) and (abs(b[1]-a[1]) < ALMOST_EQUAL_TOL) and \
            (abs(b[2]-a[2]) < ALMOST_EQUAL_TOL))

#Feed me k-point instances
def KPsAreEqual(k1, k2):
    return ((abs(k2.a-k1.a) < ALMOST_EQUAL_TOL) and (abs(k2.b-k1.b) < ALMOST_EQUAL_TOL) and \
            (abs(k2.c-k1.c) < ALMOST_EQUAL_TOL))

#Expects a and b vectors as 3d lists.  Returns 3d list
def UnitVec(a, b):
    scal = ((b[0] - a[0])**(2.) + (b[1] - a[1])**(2.) + (b[2] - a[2])**(2.))**(1./2.)
    return [(b[0] - a[0])/scal, (b[1] - a[1])/scal, (b[2] - a[2])/scal]

def VecDist(a, b):
    return ((b[0] - a[0])**(2.) + (b[1] - a[1])**(2.) + (b[2] - a[2])**(2.))**(1./2.)

def KPDist(a, b):
    return ((b.a - a.a)**(2.) + (b.b - a.b)**(2.) + (b.c - a.c)**(2.))**(1./2.)


#Checks if x is on a line defined by points a and b.  a, b, x are 3d lists.  This only works if it is used
#on a set of points that always include a fixed point
def IsOnLine(a, b, x):
    u1 = UnitVec(a, b)
    u2 = [-u for u in u1]
    xU = UnitVec(x, a)
    return(VecsAreEqual(xU, u1) or VecsAreEqual(xU, u2))

#Returns a list of lines (see format inside function) that include a specific k-point
#Wont return lines that have less points than minListLen (need at least 3 for a parabolic fit)
def GetMaximaLines(fixedKpoint, kpList, minLisLen=3):
    lines = [[fixedKpoint]] #lines = [[point1_1, point1_2, point1_3, ...], [point2_1, point2_2, point2_3, ...], ...]
    for kp in kpList:
        if(not KPsAreEqual(fixedKpoint, kp)): ##don't try to make a line with the fixed point
            for lin in lines:
                if(len(lin) < 2):
                    continue
                if(IsOnLine([lin[0].a, lin[0].b, lin[0].c], [lin[1].a, lin[1].b, lin[1].c],
                            [kp.a, kp.b, kp.c])):
                    lin.append(kp)
                    break
            lines.append([fixedKpoint, kp])
    ret = []
    for lin in lines:
        if(len(lin) >= minLisLen):
            ret.append(lin)
    return ret

#a, b, c are real-space lattice vectors (dim 3) i.e. a = a_x + a_y + a_z
pi = 3.141592654
def KListDirToCart(kList, a1, a2, a3, usc=1.):
    V = np.dot(a1, np.cross(a2, a3))
    b1 = list(2.*pi * (np.cross(a2, a3))/V)
    b2 = list(2.*pi * (np.cross(a3, a1))/V)
    b3 = list(2.*pi * (np.cross(a1, a2))/V)

    for k in kList:
        aOrig, bOrig, cOrig = k.a, k.b, k.c
        k.a = usc * (b1[0]*aOrig + b2[0]*bOrig + b3[0]*cOrig)
        k.b = usc * (b1[1]*aOrig + b2[1]*bOrig + b3[1]*cOrig)
        k.c = usc * (b1[2]*aOrig + b2[2]*bOrig + b3[2]*cOrig)
    return kList

def dispersion(x, a, c):
    return a*x**(2.) + c

#Returns the second derivitive of a cubic spline fit at x = 0.  You must give x as distances from the
#extrema, or this will fail.
def SplineFit(x, y, returnInterps=False, ptsForSpline=1000):
    def ApproxEqual(a, b):
        return (abs(a - b) < ALMOST_EQUAL_TOL)

    ##Make sure coeffs are given in high -> low order!
    def Poly(x, xn, coeffs):
        y = 0.
        for i in range(0, len(coeffs)):
            y += coeffs[i] * (x - xn) ** (float(i))
        return y

    ##For Cubic splines specifically. Make sure coeffs are given in high -> low order!
    def SecondDerivitive(xn, coeffs, evalAt):
        return 2.*coeffs[2] + 6.*coeffs[3]*(evalAt - xn)

    # Get the whole spline interpolation, plot
    spline = CubicSpline(x, y)
    xG = np.linspace(x[0], x[-1], ptsForSpline)
    yG = np.array(spline(xG))

    # Get the index of the distance = 0 point since thats where we're interested in taking the derivitive
    # zeroIndex = 0
    for i in range(0, len(x)):
        if (ApproxEqual(0., x[i])):
            zeroIndex = i
            break

    # The coeffs of interpolation a_i in (a_i * (x - x_n)^i) between datapoint #n and #n+1 are stored as:
    # a_i = spline.c.item(i, n).
    lCoeffs = [spline.c.item(i, zeroIndex - 1) for i in range(3, -1, -1)]
    rCoeffs = [spline.c.item(i, zeroIndex) for i in range(3, -1, -1)]

    # Compute the derivitives.
    dxL = SecondDerivitive(x[zeroIndex - 1], lCoeffs, evalAt=0.)
    dxR = SecondDerivitive(x[zeroIndex], rCoeffs, evalAt=0.)
    dxAvg = (dxL + dxR)/2.  ##just in case the bands are not symmetric about the extrema

    if(returnInterps):
        return dxAvg, xG, yG
    return dxAvg

#Holds information about a line in k-space
class RecipLine:
    fixedKPoint = None #k-point that is taken as the origin
    kPoints = None #non-scaled k-points in k-space
    distances = None #reciprocal-space distance of every kPoint from the fixed one.  parallel with energies
    energies = None #energies of each k-point's cbm or vbm (depending on the paramater hole / electron

    #Expects a fixed kpoint instance, a list of all kpoints, and the string "conduction" or "valence"
    #rsx are real space lattice vectors (like from POSCAR).  Scale by univ scal fact before sending these in!
    def __init__(self, fixed, pointList, bandIndex, rs1, rs2, rs3):
        self.fixedKPoint = None
        self.kPoints = []
        self.distances = []
        self.energies = []
        self.hole, self.electron = False, False

        #Remove identical points to improve fits
        newPoints = [pointList[0]]
        for point in pointList[1:]:
            copied = False
            for nps in newPoints:
                if (KPsAreEqual(point, nps)):
                    copied = True
                    break
            if(not copied):
                newPoints.append(point)

        self.fixedKPoint = deepcopy(fixed) #the fact that I have to do this is completly retarded
        self.kPoints = deepcopy(newPoints)

        #Get distances.  Involves switching from direct recip vects to cartesian recip vects
        fixed = KListDirToCart([deepcopy(fixed)], rs1, rs2, rs3)[0]
        pointsCart = list(KListDirToCart(self.kPoints, rs1, rs2, rs3))
        for point in pointsCart:
            self.distances.append(KPDist(point, fixed))
            self.energies.append(GetNthBand(bandIndex, point.bands[:])[1]) ##first index = band energy

        return

    #Returns the dist from fixed point to the furthest point from it that we care about
    def GetFurthestDist(self, npts=3, useSymm=False):
        x = sorted(self.distances[:])
        if (useSymm and npts % 2 != 0):  ##mirror the points about the smallest distance (presumably zero)
            boundLow, boundHi = int(len(x) - (npts + 1) / 2), int(len(x) + (npts - 3) / 2)
            x = [-x_ for x_ in x[::-1]] + x[1:]
            return max([abs(x_) for x_ in x[boundLow:boundHi + 1]])

        return x[npts - 1] ##if we're interested in n pts, then the n-1 index is the max dist

    #Gives the effective mass by approximating E(k) as a 2nd order polynomial
    #As of now (and probably forever) interpType=spline is only supported for useSymm=True
    def GetEffMass(self, type, npts=3, extraData=False, useSymm=False, interpType="spline"):
        dx = 1E-12
        x, y = zip(*sorted(zip(self.distances[:], self.energies[:])))  ##parallel array sort w/ black magic
        x, y = list(x), list(y)

        if(useSymm and npts%2 != 0): ##mirror the points about the smallest distance (presumably zero)
            boundLow, boundHi = int(len(x) - (npts+1)/2), int(len(x) + (npts-3)/2)
            x = [-x_ for x_ in x[::-1]] + x[1:]
            y = [y_ for y_ in y[::-1]] + y[1:]

            kMax = max([abs(x_) for x_ in x[boundLow:boundHi+1]])
            if(interpType[0] == 's' or interpType[0] == 'S'):
                secDeriv = SplineFit(deepcopy(x), deepcopy(y))
                devs = "No errors for spline interpolations :)"
            else:
                coeffs = cf(dispersion, x[boundLow:boundHi+1], y[boundLow:boundHi+1],
                            p0=[1., y[int((len(y)-1)/2)]],
                            bounds=([-np.inf, y[int((len(y)-1)/2)]], [np.inf, y[int((len(y)-1)/2)] + dx]))[0]
                #devs = [abs(dispersion(x[n],
                #                       coeffs[0], coeffs[1]) - y[n]) for n in range(boundLow, boundHi+1)]
                secDeriv = 2*coeffs[0]
        else:
            kMax = max([abs(x_) for x_ in x[:npts]])
            coeffs = cf(dispersion, x[:npts], y[:npts], p0=[1., y[0]],
                        bounds=([-np.inf, y[0]], [np.inf, y[0] + dx]))[0]
            devs = [abs(dispersion(x[n], coeffs[0], coeffs[1]) - y[n]) for n in range(0, npts)]

        mEff = (HBAR)**(2.)/(EvToJo(secDeriv))*(1E10)**(2.) ##1x10^2 for angstroms to meters

        if(extraData):
            sum = 0
            for d in devs:
                sum += d
            return mEff/m0, sum/npts, kMax
        return mEff/m0

    def GetLineDir(self, acc=3):
        uv = UnitVec([self.kPoints[0].a, self.kPoints[0].b, self.kPoints[0].c],
                       [self.kPoints[-1].a, self.kPoints[-1].b, self.kPoints[-1].c])
        return [f'{s:.3f}' for s in uv]



##TODO: Sort by whether there are poits to either side of the fixed kpoint
##      Sort by how good parabolic approxs will be (the closer to the fixed point, the better the approx)
