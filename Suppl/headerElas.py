from numpy import polyfit

#Makes stuff less annoying to type, and python has no #define equiv.
def S(s, i):
    return s.split(',')[i]

def Slope(x, y):
    return polyfit(x, y, deg=1)[0]

def GetQuadAParam(x, y):
    return polyfit(x, y, deg=2)[0]

#eV per cubic angstrom to GPa
def EvAngstToGPa(v):
    return v*160.2176487

#For parsing output directly fro linux csv
class Datapoint:
    directory = None
    wDirectory = None
    energy = None
    dEnergy = None
    a, b, c = None, None, None
    al, be, ga = None, None, None
    vol= None
    delta = None #equivalent to strain (sometimes)
    epsilons = None #equivalent to strain (all the time)
    x = None
    defID = None
    stresses = None

    #Converts kilobar (vasp default) to gpa (a normal person's default).  Also switches the sign
    #to conform to normal people standards
    def KbToGPa(self):
        for k, v in self.stresses.items():
            self.stresses[k] = -0.1*v

    #Gets delta / strain from the current directory
    def GetDelta(self, s):
        ##If it is not a zero strain model, it will contain an identifier giving the delta
        if("initRel" not in list(s.split('/'))[1:]):
            wDir = list(s.split('/'))[1:][-1]
            if(wDir[0] == 'n'):
                return -1.0 * float(wDir[1:])
            elif(wDir[0] == 'p'):
                return float(wDir[1:])
        return 0.0

    #Gets the concentration X from the directory
    def GetX(self, s):
        for piece in list(s.split('/'))[1:]:
            if(piece[0] == 'x' and len(piece) == 5): ##reasonable that this is the identifier
                return float(int(piece[1:]) / 1000)

    #Get the deform matrix ID from the full directory.
    def GetDefID(self, s):
        if("initRel" not in list(s.split('/'))[1:]):
            return int(s.split('/')[-2])
        return "Init"

    #Get the epsilon equivalents to delta
    #See https://doi.org/10.1063/1.368733 for where this awfulness comes from
    def GetEpsilons(self, delta):
        self.epsilons = {"XX": 0.0,
                         "YY": 0.0,
                         "ZZ": 0.0,
                         "YZ": 0.0,
                         "ZX": 0.0,
                         "XY": 0.0}

        fac = 1.0 / ((1.0 - delta**(2.0))**(1.0/3.0))
        ##Symmetric Deformations
        if(self.defID == 1):
            self.epsilons["XX"] = delta
        if(self.defID == 2):
            self.epsilons["YY"] = delta
        if(self.defID == 3):
            self.epsilons["ZZ"] = delta
        ##Shear Deformations (monoclinic)
        if(self.defID == 4):
            self.epsilons["XX"], self.epsilons["YY"], self.epsilons["ZZ"] = fac - 1.0, fac - 1.0, \
                                                                            fac - 1.0
            self.epsilons["YZ"] = 2.0*delta*fac
        if(self.defID == 5):
            self.epsilons["XX"], self.epsilons["YY"], self.epsilons["ZZ"] = fac - 1.0, fac - 1.0, \
                                                                fac - 1.0
            self.epsilons["ZX"] = 2.0*delta*fac
        if(self.defID == 6):
            self.epsilons["XX"], self.epsilons["YY"], self.epsilons["ZZ"] = fac - 1.0, fac - 1.0, \
                                                                            fac - 1.0
            self.epsilons["XY"] = 2.0*delta*fac
        ##Other Deformations (orthorhombic)
        inc, con, dec = (1.0 + delta)*fac - 1.0, fac - 1.0, (1.0 - delta)*fac - 1.0
        if(self.defID == 7):
            self.epsilons["XX"] = inc
            self.epsilons["YY"] = dec
            self.epsilons["ZZ"] = con
        if(self.defID == 8):
            self.epsilons["XX"] = inc
            self.epsilons["YY"] = con
            self.epsilons["ZZ"] = dec
        if(self.defID == 9):
            self.epsilons["XX"] = con
            self.epsilons["YY"] = inc
            self.epsilons["ZZ"] = dec

        return

    def __init__(self, line):
        self.directory = None
        self.wDirectory = None
        self.energy = None
        self.dEnergy = None
        self.a, self.b, self.c = None, None, None
        self.al, self.be, self.ga = None, None, None
        self.vol= None
        self.delta = None
        self.epsilons = {}
        self.x = None
        self.defID = None
        self.stresses = {}

        #Initialization from CSV
        self.directory = S(line, 0)
        self.wDirectory = S(line, 1)
        self.energy = float(S(line, 2))
        self.dEnergy = float(S(line, 3))
        self.a, self.b, self.c = float(S(line, 6)), float(S(line, 7)), float(S(line, 8))
        self.al, self.be, self.ga = float(S(line, 9)), float(S(line, 10)), float(S(line, 11))
        self.vol = float(S(line, 12))
        self.stresses["XX"], self.stresses["YY"], self.stresses["ZZ"], self.stresses["XY"], \
        self.stresses["YZ"], self.stresses["ZX"] = float(S(line, 13)), float(S(line, 14)), \
                                                   float(S(line, 15)), float(S(line, 16)), \
                                                   float(S(line, 17)), float(S(line, 18))

        #Useful information not directly from csv
        self.KbToGPa()
        self.delta = self.GetDelta(self.directory)
        self.x = self.GetX(self.directory)
        self.defID = self.GetDefID(self.directory)
        self.GetEpsilons(delta=self.delta)

#Helper function to fill Cij / Sij matrices.  Returns {index: ['c' or 's', i, j]}
def IDontKnowWhatToNameThis(s):
    ret = {}
    st = s.split(',')
    for i in range(1, len(st)):
        ret[i] = [st[i][0], int(st[i][1]), int(st[i][2])]
    return ret

#For parsing clean output Cij / Sij
class ElasSet:
    x = None
    cij = None
    sij = None

    #Line is in the form x, c11,...,c66, c12, c13, c23, s11,...,s66, s12, s13, s23
    def __init__(self, header, line):
        self.x = None
        self.cij, self.sij = [], []
        for i in range(0, 7):
            self.cij.append([])
            self.sij.append([])
            for j in range(0, 7):
                self.cij[i].append(0.0)
                self.sij[i].append(0.0)

        self.x = float(line.split(',')[0])
        #Fill Cij and Sij matrices
        matHelper = IDontKnowWhatToNameThis(header)
        for index, elem in enumerate(line.split(',')[1:], start=1):
            val = matHelper[index]
            if(val[0] == 'c'):
                self.cij[val[1]][val[2]] = float(elem)
                self.cij[val[2]][val[1]] = float(elem) #doing this since the matrices are symmetric
            elif(val[0] == 's'):
                self.sij[val[1]][val[2]] = float(elem)
                self.sij[val[2]][val[1]] = float(elem)


#Returns a multidimensional dictionary of the form:
#return = {x value 0: {deform id 0: [Datapoint 0, Datapoint 1, ...],
#                      deform id 1: [Datapoint 0, Datapoint 1, ...],
#                      ...
#                      },
#          x value 1: {deform id 0: [Datapoint 0, Datapoint 1, ...],
#                      deform id 1: [Datapoint 0, Datapoint 1, ...],
#                      ...
#                      },
#          ...
#         }
def ParseOutfile(infileLoc):
    #Get all the lines into easy to work with format
    allList = []
    with open(infileLoc, 'r') as infile:
        for num, line in enumerate(infile):
            if(num > 0): ##skip header
                allList.append(Datapoint(line))
        infile.close()

    ret = {}
    #Yes, this can be done more efficiently.  But it would be gross to look at
    #First: initialize the x value keys
    for dataPoint in allList:
        if(dataPoint.x not in ret.keys()):
            ret[dataPoint.x] = {}
    #Second: for each x value, initialize the deform matrix id dicts
    for x, xDict in ret.items():
        for dataPoint in allList:
            if(x == dataPoint.x and dataPoint.defID not in xDict.keys()):
                xDict[dataPoint.defID] = [dataPoint]
            elif(x == dataPoint.x and dataPoint.defID in xDict.keys()):
                xDict[dataPoint.defID].append(dataPoint)

    #We want the initial strain point in each deform group, so do that
    for x in ret.keys():
        for defID in ret[x].keys():
            ret[x][defID].append(ret[x]["Init"][0]) ##[0] bc there is only be one "Init"

    #Finially remove the init keys from each x value - theyre already taken care of
    for x in ret.keys():
        del ret[x]["Init"]

    return ret

#Returns a list of ElasSet instances
def ParseCijSijFile(infileLoc):
    allList = []
    with open(infileLoc, 'r') as infile:
        header = ""
        for num, line in enumerate(infile):
            if(num == 0):
                header = line
            else:
                allList.append(ElasSet(header, line))
        infile.close()

    return allList

"""
=====================================
Begin Stuff for Directional Moduli...
=====================================
"""
