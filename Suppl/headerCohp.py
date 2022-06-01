#ONLY WORKS FOR NON SPIN POL. RIGHT NOW
#READ THIS MESSAGE ^^^

#Takes a line that looks like this:   No.<ID>:<ELEM1><ELEM_ID1>-><ELEM2><ELEM_ID2>(<DIST>)
#And returns this:                    <ID>, <ELEM1>, <ELEM_ID1>, <ELEM2>, <ELEM_ID2>, <DIST>
def ParseHead(s):
    def Rip(s_):
        lets = s_.rstrip("0123456789") ##numbers are always after letters, so this is fine
        nums = s_[len(lets):]
        return lets, nums
    
    pairId = s.split(':')[0].split('.')[1]
    el1, elId1 = Rip(s.split("->")[0].split(':')[1])
    el2, elId2 = Rip(s.split('(')[0].split("->")[1])
    dist = s.split(')')[0].split('(')[1]

    return int(pairId), el1, int(elId1), el2, int(elId2), float(dist)

class Pair:
    pairId = None
    elem1, elem2 = None, None
    elemId1, elemId2 = None, None ##IDs are the absolute atom number, not relative.  
    dist = None                   ##i.e. Si15 is the 15th element, not necessairily the 15th Si  
    
    energies = None ##energy is ALREADY SCALED s.t. E = E - Efermi
    pCohps = None
    ipCohps = None ##integrated COHPs.  i.e. ipCohps[6] is integrated from energy[0] to energy[6]  
    
    #Initialize with the pair header, then have to read in energies and *cohps later.  Not my fault
    def __init__(self, line, avgAll = False):
        self.pairId = None
        self.elem1, self.elem2 = None, None
        self.elemId1, self.elemId2 = None, None
        self.dist = None
        self.energies = []
        self.pCohps = []
        self.ipCohps = []
        
        if(avgAll): ##for the average over all atom pairs
            self.pairId = 0
            self.elem1, self.elem2 = "all", "all"
            self.elemId1, self.elemId2 = 0, 0
            self.dist = -1
        else:
            self.pairId, self.elem1, self.elemId1, self.elem2, self.elemId2, self.dist = \
            ParseHead(line)
            
        
    #Adds an energy and this pair's pCohp and ipCohp from the full COHPCAR line
    #line = <nrg> <pCohp all pairs> <ipCohp all pairs> <pCohp pair 1> <ipCohp pair 1> <pcohp pair 2>  
    def AddVals(self, line):
        self.energies.append(float(line.split()[0]))
        self.pCohps.append(float(line.split()[1 + 2*self.pairId]))
        self.ipCohps.append(float(line.split()[2 + 2*self.pairId]))
        
#Returns a list of every interaction pair's COHP data
def ParseCohpFile(infileLoc, retEFermi = False):
    allPairs = []
    eFermi = -999
    with open(infileLoc, 'r') as infile:
        allPairs.append(Pair("hi, I don't do anything in this context", avgAll=True))
        for num, line in enumerate(infile):
            #Skip file header except to read in the fermi energy
            if(num < 3): 
                if(num == 1):
                    eFermi = float(line.split()[5])
            
            #Read in the pair header info
            elif(num >= 3 and line[0] == 'N'): ##every pair header begins with "No.N:..."
                allPairs.append(Pair(line))
            
            #Fill in each pair's COHP info
            else:
                for pair in allPairs:
                    pair.AddVals(line)  
        infile.close()
        
    if(retEFermi):
        return allPairs, eFermi
    return allPairs        















