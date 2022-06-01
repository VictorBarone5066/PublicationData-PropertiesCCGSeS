from numpy import linalg as la
from copy import deepcopy

M_MIN = 1E-10 ##minimum mass before considered zero

#Read file.  Returns n, p data in that order
def ReadOutfile(infileLoc, n=True, p=True):
    with open(infileLoc, 'r') as infile:
        nDat, pDat = dict(), dict()

        for lin in infile:
            line = lin.split()
            if(line[0] == "type"):
                continue

            type_ = line[0]
            mu = float(line[1])
            temp = float(line[2])
            xx = float(line[3])
            yy = float(line[4])
            zz = float(line[5])
            yz = float(line[6])
            zx = float(line[7])
            xy = float(line[8])

            if(abs(xx) < M_MIN or abs(yy) < M_MIN or abs(zz) < M_MIN or
               abs(yz) < M_MIN or abs(zx) < M_MIN or abs(xy) < M_MIN):
                continue

            if(line[0] == 'n'):
                if mu not in nDat.keys():
                    nDat[mu] = {"xx": {"temp": [], "mass": []},
                                "yy": {"temp": [], "mass": []},
                                "zz": {"temp": [], "mass": []},
                                "yz": {"temp": [], "mass": []},
                                "zx": {"temp": [], "mass": []},
                                "xy": {"temp": [], "mass": []}}
                nDat[mu]["xx"]["temp"].append(temp)
                nDat[mu]["xx"]["mass"].append(xx)
                nDat[mu]["yy"]["temp"].append(temp)
                nDat[mu]["yy"]["mass"].append(yy)
                nDat[mu]["zz"]["temp"].append(temp)
                nDat[mu]["zz"]["mass"].append(zz)
                nDat[mu]["yz"]["temp"].append(temp)
                nDat[mu]["yz"]["mass"].append(yz)
                nDat[mu]["zx"]["temp"].append(temp)
                nDat[mu]["zx"]["mass"].append(zx)
                nDat[mu]["xy"]["temp"].append(temp)
                nDat[mu]["xy"]["mass"].append(xy)

            if(line[0] == 'p'):
                if mu not in pDat.keys():
                    pDat[mu] = {"xx": {"temp": [], "mass": []},
                                "yy": {"temp": [], "mass": []},
                                "zz": {"temp": [], "mass": []},
                                "yz": {"temp": [], "mass": []},
                                "zx": {"temp": [], "mass": []},
                                "xy": {"temp": [], "mass": []}}
                pDat[mu]["xx"]["temp"].append(temp)
                pDat[mu]["xx"]["mass"].append(xx)
                pDat[mu]["yy"]["temp"].append(temp)
                pDat[mu]["yy"]["mass"].append(yy)
                pDat[mu]["zz"]["temp"].append(temp)
                pDat[mu]["zz"]["mass"].append(zz)
                pDat[mu]["yz"]["temp"].append(temp)
                pDat[mu]["yz"]["mass"].append(yz)
                pDat[mu]["zx"]["temp"].append(temp)
                pDat[mu]["zx"]["mass"].append(zx)
                pDat[mu]["xy"]["temp"].append(temp)
                pDat[mu]["xy"]["mass"].append(xy)

    if(n and p):
        return nDat, pDat
    if(n):
        return nDat
    if(p):
        return pDat

#Transform inverse effective mass tensor to normal effective mass tensor
def InvertMassTensor(dat):
    for mu in dat.keys():
        ##xx, yy, zz, ... etc's temp, mass are all of equal length (hopefully!!)
        for i in range(0, len(dat[mu]["xx"]["temp"])):
            thisTensor = [[dat[mu]["xx"]["mass"][i], dat[mu]["xy"]["mass"][i], dat[mu]["zx"]["mass"][i]],
                          [dat[mu]["xy"]["mass"][i], dat[mu]["yy"]["mass"][i], dat[mu]["yz"]["mass"][i]],
                          [dat[mu]["zx"]["mass"][i], dat[mu]["yz"]["mass"][i], dat[mu]["zz"]["mass"][i]]]
            #newTensor = list(list(t) for t in la.inv(thisTensor))
            newTensor = la.inv(thisTensor)
            dat[mu]["xx"]["mass"][i] = newTensor[0][0]
            dat[mu]["yy"]["mass"][i] = newTensor[1][1]
            dat[mu]["zz"]["mass"][i] = newTensor[2][2]
            dat[mu]["yz"]["mass"][i] = newTensor[1][2]
            dat[mu]["zx"]["mass"][i] = newTensor[0][2]
            dat[mu]["xy"]["mass"][i] = newTensor[0][1]

#Sets two extra keys in dat's entries - one for (normalized) eigenvectors, and one for the eigenvalues
#For example, dat[mu] now has three extra entries: "1", "2", "3".
#And dat[mu]["1"] = {"mass": [], "temp": [], "vec": [[v0x, v0y, v0z], ...]}, etc.
#'1', '2', '3' correspond to the lightest, medium, and heaviest mass, respestivly
def SetEigens(dat):
    for mu in dat.keys():
        dat[mu]['1'] = {"temp": dat[mu]["xx"]["temp"][:], "mass": [], "vec": []}
        dat[mu]['2'] = {"temp": dat[mu]["xx"]["temp"][:], "mass": [], "vec": []}
        dat[mu]['3'] = {"temp": dat[mu]["xx"]["temp"][:], "mass": [], "vec": []}

        for i in range(0, len(dat[mu]["xx"]["temp"])):
            thisTensor = [[dat[mu]["xx"]["mass"][i], dat[mu]["xy"]["mass"][i], dat[mu]["zx"]["mass"][i]],
                          [dat[mu]["xy"]["mass"][i], dat[mu]["yy"]["mass"][i], dat[mu]["yz"]["mass"][i]],
                          [dat[mu]["zx"]["mass"][i], dat[mu]["yz"]["mass"][i], dat[mu]["zz"]["mass"][i]]]
            va, ve = la.eig(thisTensor)
            va = list(va)
            ve = list(list(t) for t in ve)
            va, ve = (list(t) for t in zip(*sorted(zip(va, ve))))

            dat[mu]['1']["mass"].append(va[0])
            dat[mu]['1']["vec"].append(ve[0])
            dat[mu]['2']["mass"].append(va[1])
            dat[mu]['2']["vec"].append(ve[1])
            dat[mu]['3']["mass"].append(va[2])
            dat[mu]['3']["vec"].append(ve[2])

#Sets an extra key in dat's entries, for the harmonic average of the '1', '2', '3' masses.
#dat[mu]["avgs"] = {"mass": [], "temp": []}
#This is interpreted as the conductivity effective mass
def SetAvgs(dat):
    for mu in dat.keys():
        dat[mu]["cond"] = {"temp": dat[mu]["xx"]["temp"][:], "mass": []}

        for i in range(0, len(dat[mu]["xx"]["temp"])):
            dat[mu]["cond"]["mass"].append(3*(1./dat[mu]['1']["mass"][i] + \
                                              1./dat[mu]['2']["mass"][i] + \
                                              1./dat[mu]['3']["mass"][i])**(-1.))

#Sets an extra key in dat's entries, for the geometric average of the '1', '2', '3' masses.
#Defined in ashcroft and mermin, ch 28 as Mdos = (m1*m2*m3)^1/3 where m1 thru m3 refer to the princpial
#directions
#dat[mu]["doss"] = {"mass": [], "temp": []}
#This is interpreted as the density of states effective mass
def SetDoss(dat):
    for mu in dat.keys():
        dat[mu]["doss"] = {"temp": dat[mu]["xx"]["temp"][:], "mass": []}

        for i in range(0, len(dat[mu]["xx"]["temp"])):
            dat[mu]["doss"]["mass"].append((dat[mu]['1']["mass"][i] * \
                                            dat[mu]['2']["mass"][i] * \
                                            dat[mu]['3']["mass"][i])**(1./3.))


#Does everything
##Note: for very low temperatures, you may have a different number of reasonable entries for a given mu:
##In this case, averages will be incorrect (I think)!
##Also, you'll get linear algebra errors since the matrices will be written as -0.0, making a singular matrix.
def InfileToEigDat(infileLoc, n=True, p=True):
    nd, pd = None, None
    if(n):
        nd = ReadOutfile(infileLoc, True, False)
        InvertMassTensor(nd)
        SetEigens(nd)
        SetAvgs(nd)
        SetDoss(nd)
    if(p):
        pd = ReadOutfile(infileLoc, False, True)
        InvertMassTensor(pd)
        SetEigens(pd)
        SetAvgs(pd)
        SetDoss(pd)

    if(n and p):
        return nd, pd
    if(n):
        return nd
    if(p):
        return pd


#Functions to approximate the doping carrier concentration
#Requires electron and hole masses to be evaluated for the exact same mu, T...

#Due to numerical reasons, this is often not possible exactly, especially for large band gaps.
#In this case, it is useful to remember that mu has an exponentially decreasing affect on effective masses
#as it increases in distance from the band edge.
#Ex: You have a semiconductor with a 2 eV band gap.  You can sucessfully calculate m_elec at mu = 1.9,
#T = 50K.  m_hole at mu = 1.9, T = 50K isn't physically important (this is heavily n-doped!), but you need
#the numerical values for approximating the carrier concentration - it isn't possible to use the standard
#thermal average process on the valance bands at these conditions, though (e^(e - mu / kt) = e^450)!.
#Instead, use m_hole at mu = 0.5 or higher, T = 50K as an approximation to m_hole at mu = 1.9, T = 50K.
#If you have no idea where to put your approximate mu, set it at half the band gap and then move it closer
#until your calculation is possible.

#Formalism from ashcroft and mermin, ch 28.  Approximations valid for |e - mu| >> kT
#I think that |e - mu| >= 3*kT is a good start.

#all energies in eV, temperatures in kelvin, and masses scaled to the electron rest mass.
#energies and mu should be absolute (not scaled to anything).  The band gap refers to CBM - VBM,
#not the absolute energy of the middle of the gap.
from numpy import exp
from numpy import log
kB = 8.6173E-5 ##eV/K

#Gives a list of intrinsic carrier concentrations in line with the three input arrays, which are the
#dos (relative) effective masses of electrons, holes, and the temperatures they were evaluated at
#carrier concentrations given in units of 1/cm^3
def GiveIntrinsicConcs(dossN, dossP, temps, bandgap):
    ret = []
    for i in range(0, len(temps)):
        ret.append(5./2. * (dossN[i]*dossP[i])**(3./4.) * \
                   (temps[i]/300.)**(3./2.) * \
                   exp(-bandgap/(2*kB*temps[i])) * \
                   1E19)

    return ret

#Gives a list of intrinsic chemical potentials in line with three input arrays.  See above comment.
def GiveIntrinsicMus(dossN, dossP, temps, bandgap, eVBM):
    ret = []
    for i in range(0, len(temps)):
        ret.append(eVBM + bandgap/2. + 3./4.*kB*temps[i]*log(dossP[i]/dossN[i]))

    return ret

#Gives a list of doped carrier concentrations given dos masses, mus, temperatures, and the doped mu
#If the doped mu > bandgap / 2, returns n-doping concentrations.  Otherwise, returns p-doped concentrations
#Concentrations are in 1/cm^3
def GiveDopedConcs(dossN, dossP, temps, bandgap, eVBM, mu):
    iConcs = GiveIntrinsicConcs(dossN=dossN, dossP=dossP, temps=temps, bandgap=bandgap)
    iMus = GiveIntrinsicMus(dossN=dossN, dossP=dossP, temps=temps, bandgap=bandgap, eVBM=eVBM)

    cFac = 1.0 if (mu > (eVBM + bandgap/2.)) else -1.0

    ret = []
    for i in range(0, len(temps)):
        ret.append(iConcs[i]*exp(cFac*(mu - iMus[i])/(kB*temps[i])))

    return ret
