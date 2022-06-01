#For calculating charge transfer values
import headerPoscar as P

class model:
    poscar = None
    valenceCounts = None ##values defined for specific POTCARs
    totChgLis = None ##in parallel with poscar's atom list.  directly from ACF.dat
    actChgLis = None ##in parallel with poscar's atom list.  total charge - valence

    def __init__(self, pos, valenDic, totChgs):
        self.poscar = None
        self.valenceCounts = {}
        self.totChgLis = []
        self.actChgLis = []

        self.poscar = P.deepcopy(pos)
        self.valenceCounts = P.deepcopy(valenDic)
        self.totChgLis = totChgs[:]

        #Initialize the actual charges.  Assumes POSCAR's atom order has not been modified -
        #make model instances immediatly after reading in POSCAR, and then only use model's poscar
        for i, atom in enumerate(self.poscar.atoms):
            self.actChgLis.append(self.totChgLis[i] - float(self.valenceCounts[atom.atomType]))

def ParseValence(infileLoc):
    ret = {}
    with open(infileLoc, 'r') as infile:
        for line in infile:
            #Capitalize() forces elements to be read in in a way that Poscar.h can recognize
            ret[line.split()[0].capitalize()] = int(line.split()[1])
        infile.close()
    return ret

def ParseAcf(infileLoc):
    chgs = []

    with open(infileLoc, 'r') as infile:
        hyphCount = 0 ##after first hyphen, begin reading.  After second, end reading
        for line in infile:
            if("-----" in line):
                hyphCount = hyphCount + 1
                continue
            if(hyphCount == 1):
                #Format of line: atom#  X  Y  Z  totChg  minDist  atomicVol
                chgs.append(float(line.split()[4]))
            if(hyphCount == 2):
                break
        infile.close()

    return chgs[:]
