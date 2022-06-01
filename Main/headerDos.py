from matplotlib import pyplot as plt
from matplotlib import gridspec
import re

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.tight_layout = True

#Hope you don't somehow have more than 8 atom types in a POSCAR file
COLORS = [(0, 0, 0), #black
          (255, 0, 0), #red
          (0, 173, 0), #(better) green
          (0, 0, 255), #blue
          (255, 0, 183), #pink
          (15, 243, 255), #cyan
          (235, 109, 0), #orange
          (153, 0, 235) #purple
          ]
cIndex = -1

def GetColor(i):
    return (a / 255 for a in COLORS[i])

#returns a 2D array, where each element is a pair of [atom type, range], where range is in the form [low, high] inclusive
def ScanPoscar(inLoc, atomTypeLineNum = 5, atomNumLineNum = 6):
    atomTypes, atomTypeNums, ret = [], [], []
    with open(inLoc, 'r') as infile:
        for i, line in enumerate(infile):
            if (i == atomTypeLineNum):
                for a in line.split():
                    atomTypes.append(a)
            if (i == atomNumLineNum):
                for a in line.split():
                    atomTypeNums.append(int(a))
                infile.close()
                break    
    
    #Don't worry about how this works
    tot = sum(atomTypeNums)
    ret.append([atomTypes[0], [1, atomTypeNums[0]]])
    tot= tot - atomTypeNums[0]
    i = 1
    while(tot > 0):
        ret.append([atomTypes[i], [ret[i-1][1][1] + 1, ret[i-1][1][1] + atomTypeNums[i]]])
        tot = tot - atomTypeNums[i]
        i = i + 1
        
    return ret

class AtomGroup:
    path = None
    
    atomType = None
    nAtoms = None
    dosData = None
    
    spin = None
    ncl = None
    
    energy = None
    sDosSum = None
    pDosSum = None
    dDosSum = None

    sUpDosSum = None
    sDnDosSum = None    
    pUpDosSum = None
    pDnDosSum = None
    dUpDosSum = None
    dDnDosSum = None
    
    sxDosSum = None
    syDosSum = None
    szDosSum = None
    pxDosSum = None
    pyDosSum = None
    pzDosSum = None
    dxDosSum = None
    dyDosSum = None
    dzDosSum = None
    
    def __init__(self, inLoc, atomRange, nedos, atomVals = None, spin = False, ncl = False, atomType = "???"):
        if(not isinstance(atomRange, list)):
            print("AtomGroup():  atomRange (the second argument to this function) needs to be supplied as a list")
        if(len(atomRange) != 2 and atomVals == None):
            print("AtomGroup Initialize:  atomRange needs two values, in list form: [low, high].  The numbers are inclusive.\n")
            return
        if(atomRange[0] == 0 and atomVals == None):
            print("AtomGroup Initialize:  This class only supports atom data, and atom 0 is designated as the full system DOS.  The first atom is atom 1, and so on.\n")
            return
        if(spin == True and ncl == True):
            print("Spin up/down does not make sense in terms of noncollinear mag.  Only one may be turned on at a time")
            return
        
        self.path = inLoc
        self.atomType = atomType
        self.dosData = []
        self.spin = spin
        self.ncl = ncl

        if (atomVals == None):
            i = atomRange[0]
            self.nAtoms = (atomRange[-1] - atomRange[0] + 1)
            energies, dosS, dosP, dosD, dosSu, dosPu, dosDu, dosSd, dosPd, dosDd = [], [], [], [], [], [], [], [], [], []
            dosSx, dosSy, dosSz, dosPx, dosPy, dosPz, dosDx, dosDy, dosDz = [], [], [], [], [], [], [], [], []
            while(i <= atomRange[-1]):
                thisAtom = GetAtomDosInfo(inLoc, i, self.spin, self.ncl)
                self.dosData.append(thisAtom)
                energies.append(thisAtom["energy"])
    
                if (not self.spin):
                    dosS.append(thisAtom["sDos"])
                    dosP.append(thisAtom["pDos"])
                    dosD.append(thisAtom["dDos"])    
                if(self.spin and not self.ncl):
                    dosSu.append(thisAtom["sDos(up)"])
                    dosPu.append(thisAtom["pDos(up)"])
                    dosDu.append(thisAtom["dDos(up)"])                 
                    dosSd.append(thisAtom["sDos(dn)"])
                    dosPd.append(thisAtom["pDos(dn)"])
                    dosDd.append(thisAtom["dDos(dn)"])             
                if(not self.spin and self.ncl):
                    dosSx.append(thisAtom["sx"])
                    dosSy.append(thisAtom["sy"])
                    dosSz.append(thisAtom["sz"])
                    dosPx.append(thisAtom["px"])
                    dosPy.append(thisAtom["py"])
                    dosPz.append(thisAtom["pz"])
                    dosDx.append(thisAtom["dx"])
                    dosDy.append(thisAtom["dy"])
                    dosDz.append(thisAtom["dz"])                    
                i = i + 1
    
            i = 0
            self.sDosSum, self.pDosSum, self.dDosSum = [0]*nedos, [0]*nedos, [0]*nedos
            self.sUpDosSum, self.pUpDosSum, self.dUpDosSum = [0]*nedos, [0]*nedos, [0]*nedos
            self.sDnDosSum, self.pDnDosSum, self.dDnDosSum = [0]*nedos, [0]*nedos, [0]*nedos  
            self.sxDosSum, self.syDosSum, self.szDosSum = [0]*nedos, [0]*nedos, [0]*nedos
            self.pxDosSum, self.pyDosSum, self.pzDosSum = [0]*nedos, [0]*nedos, [0]*nedos
            self.dxDosSum, self.dyDosSum, self.dzDosSum = [0]*nedos, [0]*nedos, [0]*nedos
            self.energy = energies[0]
    
            while(i < (atomRange[-1] - atomRange[0] + 1)):
                j = 0
                while(j < len(energies[i])):
                    if (not self.spin):
                        self.sDosSum[j] = self.sDosSum[j] + dosS[i][j]
                        self.pDosSum[j] = self.pDosSum[j] + dosP[i][j]   
                        self.dDosSum[j] = self.dDosSum[j] + dosD[i][j]
                    if(self.spin and not self.ncl):
                        self.sUpDosSum[j] = self.sUpDosSum[j] + dosSu[i][j]
                        self.pUpDosSum[j] = self.pUpDosSum[j] + dosPu[i][j]   
                        self.dUpDosSum[j] = self.dUpDosSum[j] + dosDu[i][j]
                        self.sDnDosSum[j] = self.sDnDosSum[j] + dosSd[i][j]
                        self.pDnDosSum[j] = self.pDnDosSum[j] + dosPd[i][j]   
                        self.dDnDosSum[j] = self.dDnDosSum[j] + dosDd[i][j]    
                    if(not self.spin and self.ncl):
                        self.sxDosSum[j] = self.sxDosSum[j] + dosSx[i][j]
                        self.syDosSum[j] = self.syDosSum[j] + dosSy[i][j]
                        self.szDosSum[j] = self.szDosSum[j] + dosSz[i][j]
                        self.pxDosSum[j] = self.pxDosSum[j] + dosPx[i][j]
                        self.pyDosSum[j] = self.pyDosSum[j] + dosPy[i][j]
                        self.pzDosSum[j] = self.pzDosSum[j] + dosPz[i][j]
                        self.dxDosSum[j] = self.dxDosSum[j] + dosDx[i][j]
                        self.dyDosSum[j] = self.dyDosSum[j] + dosDy[i][j]
                        self.dzDosSum[j] = self.dzDosSum[j] + dosDz[i][j]
                    j = j + 1
                i = i + 1        

        if (atomVals != None):
            self.nAtoms = len(atomVals)
            energies, dosS, dosP, dosD, dosSu, dosPu, dosDu, dosSd, dosPd, dosDd = [], [], [], [], [], [], [], [], [], []
            dosSx, dosSy, dosSz, dosPx, dosPy, dosPz, dosDx, dosDy, dosDz = [], [], [], [], [], [], [], [], []
            for i in atomVals:
                thisAtom = GetAtomDosInfo(inLoc, i, self.spin, self.ncl)
                self.dosData.append(thisAtom)
                energies.append(thisAtom["energy"])
    
                if (not self.spin):
                    dosS.append(thisAtom["sDos"])
                    dosP.append(thisAtom["pDos"])
                    dosD.append(thisAtom["dDos"])    
                if(self.spin and not self.ncl):
                    dosSu.append(thisAtom["sDos(up)"])
                    dosPu.append(thisAtom["pDos(up)"])
                    dosDu.append(thisAtom["dDos(up)"])                 
                    dosSd.append(thisAtom["sDos(dn)"])
                    dosPd.append(thisAtom["pDos(dn)"])
                    dosDd.append(thisAtom["dDos(dn)"])             
                if(not self.spin and self.ncl):
                    dosSx.append(thisAtom["sx"])
                    dosSy.append(thisAtom["sy"])
                    dosSz.append(thisAtom["sz"])
                    dosPx.append(thisAtom["px"])
                    dosPy.append(thisAtom["py"])
                    dosPz.append(thisAtom["pz"])
                    dosDx.append(thisAtom["dx"])
                    dosDy.append(thisAtom["dy"])
                    dosDz.append(thisAtom["dz"])   
    
            self.sDosSum, self.pDosSum, self.dDosSum = [0]*nedos, [0]*nedos, [0]*nedos
            self.sUpDosSum, self.pUpDosSum, self.dUpDosSum = [0]*nedos, [0]*nedos, [0]*nedos
            self.sDnDosSum, self.pDnDosSum, self.dDnDosSum = [0]*nedos, [0]*nedos, [0]*nedos        
            self.sxDosSum, self.syDosSum, self.szDosSum = [0]*nedos, [0]*nedos, [0]*nedos
            self.pxDosSum, self.pyDosSum, self.pzDosSum = [0]*nedos, [0]*nedos, [0]*nedos
            self.dxDosSum, self.dyDosSum, self.dzDosSum = [0]*nedos, [0]*nedos, [0]*nedos
            self.energy = energies[0]

            i = 0    
            while(i < len(atomVals)):
                j = 0
                while(j < len(energies[i])):
                    if (not self.spin):
                        self.sDosSum[j] = self.sDosSum[j] + dosS[i][j]
                        self.pDosSum[j] = self.pDosSum[j] + dosP[i][j]   
                        self.dDosSum[j] = self.dDosSum[j] + dosD[i][j]
                    if(self.spin):
                        self.sUpDosSum[j] = self.sUpDosSum[j] + dosSu[i][j]
                        self.pUpDosSum[j] = self.pUpDosSum[j] + dosPu[i][j]   
                        self.dUpDosSum[j] = self.dUpDosSum[j] + dosDu[i][j]
                        self.sDnDosSum[j] = self.sDnDosSum[j] + dosSd[i][j]
                        self.pDnDosSum[j] = self.pDnDosSum[j] + dosPd[i][j]   
                        self.dDnDosSum[j] = self.dDnDosSum[j] + dosDd[i][j]      
                    if(not self.spin and self.ncl):
                        self.sxDosSum[j] = self.sxDosSum[j] + dosSx[i][j]
                        self.syDosSum[j] = self.syDosSum[j] + dosSy[i][j]
                        self.szDosSum[j] = self.szDosSum[j] + dosSz[i][j]
                        self.pxDosSum[j] = self.pxDosSum[j] + dosPx[i][j]
                        self.pyDosSum[j] = self.pyDosSum[j] + dosPy[i][j]
                        self.pzDosSum[j] = self.pzDosSum[j] + dosPz[i][j]
                        self.dxDosSum[j] = self.dxDosSum[j] + dosDx[i][j]
                        self.dyDosSum[j] = self.dyDosSum[j] + dosDy[i][j]
                        self.dzDosSum[j] = self.dzDosSum[j] + dosDz[i][j]
                    j = j + 1
                i = i + 1

    def PlotThisAtom(self, suppressS = False, suppressP = False, suppressD = False):
        global cIndex
        cIndex = (cIndex + 1)%8
        if (not self.spin and not self.ncl):
            if(not suppressS):
                plt.plot([e - self.dosData[0]["eFermi"] for e in self.energy], [s/self.nAtoms for s in self.sDosSum], color=tuple(GetColor(cIndex)), label=self.atomType+" "+'(s)', linewidth=.8)
            if(not suppressP):
                plt.plot([e - self.dosData[0]["eFermi"] for e in self.energy], [p/self.nAtoms for p in self.pDosSum], color=tuple(GetColor(cIndex)), label=self.atomType+" "+'(p)', linewidth=.8, linestyle='--')
            if (not suppressD):
                plt.plot([e - self.dosData[0]["eFermi"] for e in self.energy], [d/self.nAtoms for d in self.dDosSum], color=tuple(GetColor(cIndex)), label=self.atomType+" "+'(d)', linewidth=.8, linestyle=':')
        if (self.spin and not self.ncl):
            if(not suppressS):
                plt.plot([e - self.dosData[0]["eFermi"] for e in self.energy], [s/self.nAtoms for s in self.sUpDosSum], color=tuple(GetColor(cIndex)), label=self.atomType+" "+'(s)', linewidth=.8)
            if(not suppressP):
                plt.plot([e - self.dosData[0]["eFermi"] for e in self.energy], [p/self.nAtoms for p in self.pUpDosSum], color=tuple(GetColor(cIndex)), label=self.atomType+" "+'(p)', linewidth=.8, linestyle='--')
            if (not suppressD):
                plt.plot([e - self.dosData[0]["eFermi"] for e in self.energy], [d/self.nAtoms for d in self.dUpDosSum], color=tuple(GetColor(cIndex)), label=self.atomType+" "+'(d)', linewidth=.8, linestyle=':')
            if(not suppressS):
                plt.plot([e - self.dosData[0]["eFermi"] for e in self.energy], [-s/self.nAtoms for s in self.sDnDosSum], color=tuple(GetColor(cIndex)), label=self.atomType+" "+'(s)', linewidth=.8)
            if(not suppressP):
                plt.plot([e - self.dosData[0]["eFermi"] for e in self.energy], [-p/self.nAtoms for p in self.pDnDosSum], color=tuple(GetColor(cIndex)), label=self.atomType+" "+'(p)', linewidth=.8, linestyle='--')
            if (not suppressD):
                plt.plot([e - self.dosData[0]["eFermi"] for e in self.energy], [-d/self.nAtoms for d in self.dDnDosSum], color=tuple(GetColor(cIndex)), label=self.atomType+" "+'(d)', linewidth=.8, linestyle=':')                


def GetNEDOS(inLoc):
    with open(inLoc, 'r') as infile:
        for i, line in enumerate(infile):
            if (i == 5):
                infile.close()
                break
    return int((line.split())[2]) 

def GetLineRange(atomNum, NEDOS):    
    if(atomNum < 1):
        return [5, 5+NEDOS]
    #[low, high]
    return [5+atomNum*(1+NEDOS), 5+atomNum+(atomNum+1)*NEDOS]
    

def ReadDoscarHead(inLoc):
    ret = {"nAtoms": None,
           "PDOS": None}
    
    lines = []
    with open(inLoc, 'r') as infile:
        for line in infile:
            lines.append(line)
            if(len(lines) > 4):
                infile.close()
                break    

    ret["nAtoms"] = int((lines[0].split())[0])
    ret["PDOS"] = bool((lines[0].split())[2])
    return ret
    
def GetEnergyInfo(inLoc, atomNum):
    ret = {"eMax": None,
           "eMin": None,
           "NEDOS": None,
           "eFermi": None}
    
    lineRange = GetLineRange(atomNum, GetNEDOS(inLoc))

    with open(inLoc, 'r') as infile:
        for i, line in enumerate(infile):
            if (i == lineRange[0]):
                infile.close()
                break
            
    ret["eMax"] = float((line.split())[0])
    ret["eMin"] = float((line.split())[1])
    ret["NEDOS"] = int((line.split())[2])        
    ret["eFermi"] = float((line.split())[3])
    return ret

def GetAtomDosInfo(inLoc, atomNum, spin = False, ncl = False):
    #Need cases for spin-pol and cases for getting specific atoms or the entire system.
    #General:
    ret = {"type" : None, #all cases.  For determining what is inside ret.  see below
           "path" : None, #all cases
           "eFermi": None, #all cases
           "energy": None, #all cases
           "dos": None, #entire system, spin pol off
           "dos(up)" : None, #entire system, spin pol on
           "dos(dn)" : None, #entire system, spin pol on
           "intDos": None, #entire system, spin pol on
           "intDos(up)" : None, #entire system, spin pol on
           "intDos(dn)" : None, #entire system, spin pol on      
           "sDos" : None, #specific atom, spin pol off
           "pDos" : None, #specific atom, spin pol off
           "dDos" : None, #specific atom, spin pol off
           "sDos(up)" : None,#<-\
           "sDos(dn)" : None,#<- \
           "pDos(up)" : None,#<-  \specific atom, spin pol on
           "pDos(dn)" : None,#<-  /
           "dDos(up)" : None,#<- /
           "dDos(dn)" : None,#<-/
           "stot" : None,   #\
           "sx" : None,     # \
           "sy" : None,     #  \
           "sz" : None,     #   \
           "ptot" : None,   #    \
           "px" : None,     #     \ -Noncollinear l-decomposed (NOT lm decomposed)
           "py" : None,     #     / -The (x, y, z) are magnetization densities
           "pz" : None,     #    /  -stot, ptot, dtot are actually the same as sDos, pDos, dDos but have
           "dtot" : None,   #   /   different variables because i'm dumb
           "dx" : None,     #  /
           "dy" : None,     # /
           "dz" : None      #/
           }  
    
    eInfo = GetEnergyInfo(inLoc, atomNum)
    NEDOS_ = eInfo["NEDOS"]
    eFermi_ = eInfo["eFermi"]   
    lineRange = GetLineRange(atomNum, NEDOS_)
    
    ret["path"] = inLoc
    ret["eFermi"] = eFermi_ 
    
    #TYPE 0: entire system, spin pol off.  Also works with ncl:
    if(atomNum < 1 and not spin):
        e, d, iD = [], [], [] #energy, dos, integrated dos
        infile = open(inLoc, 'r')
        for i, line in enumerate(infile):
            if (lineRange[0] < i <= lineRange[1]):
                e.append(float(line.split()[0]))
                d.append(float(line.split()[1]))    
                iD.append(float(line.split()[2])) 
        infile.close()
       
        ret["type"] = 0
        ret["energy"] = e    
        ret["dos"] = d
        ret["intDos"] = iD
    
    #TYPE 1: entire system, spin pol on:
    if(atomNum < 1 and spin and not ncl):
        e, du, dd, iDu, iDd = [], [], [], [], [] #energy, dos up & down, int dos up & down
        infile = open(inLoc, 'r')
        for i, line in enumerate(infile):
            if (lineRange[0] < i <= lineRange[1]):
                e.append(float(line.split()[0]))
                du.append(float(line.split()[1]))    
                dd.append(float(line.split()[2])) 
                iDu.append(float(line.split()[3]))    
                iDd.append(float(line.split()[4]))             
        infile.close()
       
        ret["type"] = 1
        ret["energy"] = e    
        ret["dos(up)"] = du
        ret["dos(dn)"] = dd
        ret["intDos(up)"] = iDu   
        ret["intDos(dn)"] = iDd     

    #TYPE 2: specific atom, spin pol off:
    if(atomNum >= 1 and not spin and not ncl):    
        e, s, p, d = [], [], [], [] #energy, s-dos, p-dos, d-dos
        infile = open(inLoc, 'r')
        for i, line in enumerate(infile):
            if (lineRange[0] < i <= lineRange[1]):
                e.append(float(line.split()[0]))
                s.append(float(line.split()[1]))    
                p.append(float(line.split()[2])) 
                d.append(float(line.split()[3]))                
        infile.close()
       
        ret["type"] = 2
        ret["energy"] = e    
        ret["sDos"] = s
        ret["pDos"] = p
        ret["dDos"] = d   
    
    #TYPE 3: specific atom, spin pol on: (yikes)
    if(atomNum >= 1 and spin and not ncl):
        e, su, sd, pu, pd, du, dd = [], [], [], [], [], [], [] #energy, s-dos up/dn, p-dos up/dn, d-dos up/dn
        infile = open(inLoc, 'r')
        for i, line in enumerate(infile):
            if (lineRange[0] < i <= lineRange[1]):
                e.append(float(line.split()[0]))
                su.append(float(line.split()[1]))    
                sd.append(float(line.split()[2])) 
                pu.append(float(line.split()[3]))  
                pd.append(float(line.split()[4]))    
                du.append(float(line.split()[5])) 
                dd.append(float(line.split()[6]))               
        infile.close()
       
        ret["type"] = 3
        ret["energy"] = e    
        ret["sDos(up)"] = su
        ret["sDos(dn)"] = sd
        ret["pDos(up)"] = pu
        ret["pDos(dn)"] = pd 
        ret["dDos(up)"] = du
        ret["dDos(dn)"] = dd    
    
    #TYPE 4: specific atom, spin pol off, ncl on (no f orbital yet): (super yikes)
    if(atomNum >= 1 and spin and not ncl):
        e, s, sx, sy, sz, p, px, py, pz, d, dx, dy, dz = [], [], [], [], [], [], [], [], [], [], [], [] #energy, sdos(x, y, z), pdos(x, y, z), ddos(x, y, z)
        infile = open(inLoc, 'r')
        for i, line in enumerate(infile):
            if (lineRange[0] < i <= lineRange[1]):
                e.append(float(line.split()[0]))
                s.append(float(line.split()[1]))    
                sx.append(float(line.split()[2])) 
                sy.append(float(line.split()[3]))  
                sz.append(float(line.split()[4]))    
                p.append(float(line.split()[5])) 
                px.append(float(line.split()[6]))         
                py.append(float(line.split()[7]))  
                pz.append(float(line.split()[8]))    
                d.append(float(line.split()[9])) 
                dx.append(float(line.split()[10]))
                dy.append(float(line.split()[11]))  
                dz.append(float(line.split()[12]))    
        infile.close()
       
        ret["type"] = 4
        ret["energy"] = e    
        ret["stot"] = s
        ret["sx"] = sx
        ret["sy"] = sy
        ret["sz"] = sz 
        ret["ptot"] = p
        ret["px"] = px
        ret["py"] = py
        ret["pz"] = pz  
        ret["ptot"] = d
        ret["dx"] = dx
        ret["dy"] = dy
        ret["dz"] = dz  
        
        ret["sDos"] = s
        ret["pDos"] = p
        ret["dDos"] = d
        
    return ret

def GetDisplayName(ret):
    if(re.search("pure", ret)):
        return "Pure"
    if(re.search("sv", ret)):
        return "SV"
    if(re.search("dv", ret) and not re.search("Skew", ret)):
            return r"DV$^{a}$"
    if(re.search("dv", ret) and re.search("Skew", ret)):
            return r"DV$^{b}$"
    if(re.search("555777", ret)):
        return "555777"
    if(re.search("sws", ret) and re.search("Skew", ret)):
            return r"STW$^{a}$"
    if(re.search("sws", ret) and not re.search("Skew", ret)):
            return r"STW$^{b}$"
    return "???"

def GetDirection(key):
    if(re.search("xStrain", key)):
        return 'X'
    if(re.search("yStrain", key)):
        return 'Y'    
    if(re.search("minEnergy", key)):
        return 'Min'
    
def GetAxisName(path):
    if(GetDirection(path) == 'X'):
        return "2% Zig-Zag Strain"
    if(GetDirection(path) == 'Y'):
        return "2% Armchair Strain"
    if(GetDirection(path) == 'Min'):
        return "Zero Strain"    
    
def DOSGridGraphAtom(pathToIn, pathToOut, outType, fileList, suppressS = False, suppressP = False, suppressD = False):
    fig = plt.figure(figsize=(8.5, 11))
    gs = gridspec.GridSpec(7, 3, fig, hspace=0.1)
    
    def GetRowCol(path):
        
        #Determine the correct row
        if(GetDisplayName(path) == "Pure"):
            row = 0
        elif(GetDisplayName(path) == "SV"):
            row = 1
        elif(GetDisplayName(path) == r"DV$^{a}$"):
            row = 2
        elif(GetDisplayName(path) == r"DV$^{b}$"):
            row = 3
        elif(GetDisplayName(path) == "555777"):
            row = 4
        elif(GetDisplayName(path) == r"STW$^{a}$"):
            row = 5
        elif(GetDisplayName(path) == r"STW$^{b}$"):
            row = 6

        
        #Determine the correct column
        if(GetDirection(path) == 'Min'):
            col = 0
        if(GetDirection(path) == 'X'):
            col = 1            
        if(GetDirection(path) == 'Y'):
            col = 2      
            
        return row, col

    for f in fileList:
        row, col = GetRowCol(f.path)
        ax = fig.add_subplot(gs[row, col])
        ax.set_facecolor((0,0,0,0))
        ax.set_xlim(-5, 5)
        if(not f.spin):
            ax.set_ylim(0, 25)
        else:
            ax.set_ylim(-25, 25)
        
        if(col == 2):
            ax.set_title(GetDisplayName(f.path), rotation = 90, x=1.1, y=0.4)       
        if(row == 6):
            ax.set_xlabel(GetAxisName(f.path))
        else:                    
            ax.tick_params('x', labelbottom=False)
        if(row == 3 and col == 0):
            ax.set_ylabel('Electronic DOS (|states| / eV)', labelpad=10)  
        if(row == 0 and col == 1):
            ax.set_title(r"$E-E_\mathrm{F}$ (eV)")              
        if(col == 1 or col == 2):
            ax.tick_params('y', labelleft=False)
      
        if (not f.spin):
            if(not suppressS):
                ax.plot([e - f.dosData[0]["eFermi"] for e in f.energy], f.sDosSum, color='k', label='2s', linewidth=.8)
            if(not suppressP):
                ax.plot([e - f.dosData[0]["eFermi"] for e in f.energy], f.pDosSum, color='k', label='2p', linewidth=.8, linestyle='--')
            if (not suppressD):
                ax.plot([e - f.dosData[0]["eFermi"] for e in f.energy], f.dDosSum, color='k', label='2d', linewidth=.8, linestyle=':')
        if (f.spin):
            if(not suppressS):
                ax.plot([e - f.dosData[0]["eFermi"] for e in f.energy], f.sUpDosSum, color='k', label='2s', linewidth=.8)
            if(not suppressP):
                ax.plot([e - f.dosData[0]["eFermi"] for e in f.energy], f.pUpDosSum, color='k', label='2p', linewidth=.8, linestyle='--')
            if (not suppressD):
                ax.plot([e - f.dosData[0]["eFermi"] for e in f.energy], f.dUpDosSum, color='k', label='2d', linewidth=.8, linestyle=':')
            if(not suppressS):
                ax.plot([e - f.dosData[0]["eFermi"] for e in f.energy], [-s for s in f.sDnDosSum], color='k', linewidth=.8)
            if(not suppressP):
                ax.plot([e - f.dosData[0]["eFermi"] for e in f.energy], [-p for p in f.pDnDosSum], color='k', linewidth=.8, linestyle='--')
            if (not suppressD):
                ax.plot([e - f.dosData[0]["eFermi"] for e in f.energy], [-d for d in f.dDnDosSum], color='k', linewidth=.8, linestyle=':')                
    
    fig.savefig(str(pathToOut), dpi = 300)

def DOSGridGraphFull(pathToIn, pathToOut, outType, fileList):
    fig = plt.figure(figsize=(8.5, 11))
    gs = gridspec.GridSpec(7, 3, fig, hspace=0.1)
    
    def GetRowCol(path):
        
        #Determine the correct row
        if(GetDisplayName(path) == "Pure"):
            row = 0
        elif(GetDisplayName(path) == "SV"):
            row = 1
        elif(GetDisplayName(path) == r"DV$^{a}$"):
            row = 2
        elif(GetDisplayName(path) == r"DV$^{b}$"):
            row = 3
        elif(GetDisplayName(path) == "555777"):
            row = 4
        elif(GetDisplayName(path) == r"STW$^{a}$"):
            row = 5
        elif(GetDisplayName(path) == r"STW$^{b}$"):
            row = 6

        
        #Determine the correct column
        if(GetDirection(path) == 'Min'):
            col = 0
        if(GetDirection(path) == 'X'):
            col = 1            
        if(GetDirection(path) == 'Y'):
            col = 2      
            
        return row, col

    for f in fileList:
        row, col = GetRowCol(f["path"])
        ax = fig.add_subplot(gs[row, col])
        ax.set_facecolor((0,0,0,0))
        ax.set_xlim(-10, 5)
        if(not f["type"] == 1):
            ax.set_ylim(0, 45)
        else:
            ax.set_ylim(-25, 25)
        
        if(col == 2):
            ax.set_title(GetDisplayName(f["path"]), rotation = 90, x=1.1, y=0.4)       
        if(row == 6):
            ax.set_xlabel(GetAxisName(f["path"]))
        else:                    
            ax.tick_params('x', labelbottom=False)
        if(row == 3 and col == 0):
            ax.set_ylabel('Electronic DOS (states / eV)', labelpad=10)  
        if(row == 0 and col == 1):
            ax.set_title(r"$E-E_\mathrm{F}$ (eV)")              
        if(col == 1 or col == 2):
            ax.tick_params('y', labelleft=False)
      
        if (not f["type"] == 1):
            ax.plot([e - f["eFermi"] for e in f["energy"]], f["dos"], color='k', linewidth=.5)
        if (f["type"] == 1):
            ax.plot([e - f["eFermi"] for e in f["energy"]], f["dos(up)"], color='k', linewidth=.5)
            ax.plot([e - f["eFermi"] for e in f["energy"]], [-d for d in f["dos(dn)"]], color='k', 
                    linewidth=.5)   
  
    fig.savefig(str(pathToOut), dpi = 300)
    
def Difference(datA, datB):
    ret = []
    i = 0
    while(i < len(datA)):
        ret.append(datA[i] - datB[i])
        i = i + 1
    return ret
    
    
    
    
    
    
    
    