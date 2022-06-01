import numpy as np

TOL = 0.000001
def Approx(a, b):
    return abs(a - b) < TOL

def AlmostEqual(u, v):
    return ((Approx(u[0], v[0])) and (Approx(u[1], v[1])) and (Approx(u[2], v[2])))

class kPoint:
    idChar = None
    a, b, c = None, None, None
    weight = None
    bands = None

    def __init__(self, header, bandInfo):
        self.idChar = None #set this later - no way to obtain from EIGENVAL file alone
        self.a, self.b, self.c = None, None, None
        self.weight = None
        self.bands = []

        self.a, self.b, self.c = float(header.split()[0]), float(header.split()[1]), \
                                 float(header.split()[2])
        self.weight = float(header.split()[3])
        for i in range(0, len(bandInfo)):
            self.bands.append([int(bandInfo[i].split()[0]),     #Band ID
                               float(bandInfo[i].split()[1]),   #Band Energy
                               float(bandInfo[i].split()[2])])  #Occupancy (not int bc of part. occ)

    def SetId(self, char):
        self.idChar = char

    def RemoveLowOcc(self, rmBelow=0.5):
        newBands = []
        for b in self.bands:
            if(b[2] > rmBelow):
                newBands.append(b)
        self.bands = newBands[:]

    #As of now, only shifts the c.b. up by specified amount
    def ApplyScissor(self, scissor = 0, tol = 0.5):
        for b in self.bands:
            if(b[2] < tol): #if occupancy less than tolerance
                b[1] = b[1] + scissor #increase energy by scissor in eV

class LineSegment:
    startID, endID = None, None
    start, end = None, None
    kPoints = None
    length = None
    scaledLength = None

    scaledKPoints = None

    def __init__(self, kplis, st = None, ed = None):
        self.startID, self.endID = None, None
        self.start, self.end = None, None
        self.kPoints = []
        self.length = None
        self.scaledLength = None
        self.scaledKPoints = {}


        self.startID, self.endID = kplis[0].idChar, kplis[-1].idChar
        self.kPoints = kplis
        self.length = ((kplis[-1].a - kplis[0].a)**2 + (kplis[-1].b - kplis[0].b)**2
                       + (kplis[-1].c - kplis[0].c)**2)**(1./2.)

        #Safer way to initialize the length - guarentees that unordered kpts wont be an issue
        if(st != None and ed != None):
            self.start, self.end = st, ed
            self.startID, self.endID = st.idChar, ed.idChar
            self.length = ((ed.a - st.a)**2 + (ed.b - st.b)**2 + (ed.c - st.c)**2)**(1./2.)

    def GetScaledKpts(self, totalScaleFactor = 1, xOffset = 0):
        for k in self.kPoints:
            x = ((k.a - self.start.a)**2 + (k.b - self.start.b)**2 +
                 (k.c - self.start.c)**2)**(1./2.) / self.length
            x = xOffset + x*totalScaleFactor
            self.scaledKPoints[x] = k

#Returns a list of k-point instances from the EIGENVAL file
def ParseEigenval(infileLoc, occupanciesListed=True, eFermi=None):
    if((not occupanciesListed) and (eFermi == None)):
        print("ParseEigenval: ERR: if the occupancies are not explicitly written to EIGENVAL, you need to")
        print("provide the fermi energy to this function!")
        return [None]

    kpts = []

    with open(infileLoc, 'r') as infile:
        readingBands = False
        head, info = None, []
        for num, line in enumerate(infile):
            #Obtain important params before reading in the k-points
            if(num == 5):
                ##nElectrons = int(line.split()[0])
                nKpts = int(line.split()[1])
                nBands = int(line.split()[2])
                bandCounter = int(line.split()[2])

            #Check for reading k-point heads
            if(num > 5 and not readingBands and len(line.split()) == 4):
                head = line
                readingBands = True

            #Check for reading in that k-point's band info
            if(occupanciesListed):
                if(num > 5 and readingBands and len(line.split()) == 3):
                    info.append(line)
                    bandCounter = bandCounter - 1
            else:
                if(num > 5 and readingBands and len(line.split()) == 2):
                    if(float(line.split()[1]) < eFermi): ##if e < eFermi, assume fully occupied
                        line += " 1.000000"
                    elif(float(line.split()[1]) > eFermi): ##if e > eFermi, assume completly unoccupied
                        line += " 0.000000"
                    else:
                        print("ParseEigenval: Uh oh: State located at eFermi\n")

                    info.append(line)
                    bandCounter = bandCounter - 1

            #Check for ending that k-point's info
            if(num > 5 and readingBands and (len(line.split()) == 0 or bandCounter == 0)):
                kpts.append(kPoint(head, info))
                bandCounter = nBands
                head, info = None, []
                readingBands = False

        infile.close()

        if(len(kpts) == nKpts):
            return kpts
        else:
            print("I couldn't read in all the k-points.  SAD!!!")
        return []


#Reads in a line-mode style KPOINTS file (presumably from pymatgen) and edits an existing kpoints
#list to add band structure labels to them
def LabelKpoints(infileLoc, kpList):
    with open(infileLoc, 'r') as infile:
        go = False
        for line in infile:
            ##Skip header junk
            if(not go and (line.split()[0][0] == 'R' or line.split()[0][0] == 'r')):
                go = True
                continue

            ##If there are actually coordinates on this line...
            if(go and '!' in line.split()):
                #go through the k point list and set the correct labels
                label = line.split()[4]
                for i in range(0, len(kpList)):
                    if(AlmostEqual([float(k) for k in line.split()[0:3]],
                                   [kpList[i].a, kpList[i].b, kpList[i].c])):
                        kpList[i].SetId(label)
                        break

        infile.close()

#If you ever find yourself having to debug this function, it'd be more time efficient to start
#over from scratch.  Sorry - I was in a rush
from copy import deepcopy
def GetLineSegments(linemodeLoc, kpLis):
    #Get initial grid of kpoints
    segs = []
    with open(linemodeLoc, 'r') as infile:
        go = False
        counter = 0
        start, end = [], []
        intersects = 0
        for num, line in enumerate(infile):
            ##Get number of intersections per line segment
            if(num == 1):
                intersects = int(line)
            ##Skip the rest of the header
            if(not go and (line.split()[0][0] == 'R' or line.split()[0][0] == 'r')):
                go = True
                continue

            ##If there are actually coordinates on this line...
            if(go and '!' in line.split()):
                ###Read the beginning of the line segment
                if(counter%2 == 0):
                    start = [float(k) for k in line.split()[0:3]]
                ###Read the end of the line segment, add the line segment to a list...
                if(counter%2 == 1):
                    end = [float(k) for k in line.split()[0:3]]
                    thisKpLis = []
                    for a, b, c in zip(np.linspace(start[0], end[0], intersects),
                                       np.linspace(start[1], end[1], intersects),
                                       np.linspace(start[2], end[2], intersects)):
                        thisKpLis.append([a, b, c])
                    ####...by searching through every k-point and matching it.  Yikes
                    thisSegLis = []
                    start_, end_ = None, None
                    for tkl in thisKpLis:
                        for k in kpLis:
                            if(AlmostEqual(start, [k.a, k.b, k.c])):
                                start_ = k
                            if(AlmostEqual(end, [k.a, k.b, k.c])):
                                end_ = k
                            if(AlmostEqual(tkl, [k.a, k.b, k.c])):
                                thisSegLis.append(k)
                                break
                    segs.append(LineSegment(thisSegLis[:], start_, end_))

                counter = counter + 1

        infile.close()
    return segs
