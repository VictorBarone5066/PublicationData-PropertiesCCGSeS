class line:
    h = None
    k = None
    l = None
    d = None
    fr = None
    fi = None
    fAbs = None
    ang = None
    intsty = None
    m = None
    wavID = None
    phase = None
    
    def __init__(self, lin):
        self.h = None
        self.k = None
        self.l = None
        self.d = None
        self.fr = None
        self.fi = None
        self.fAbs = None
        self.ang = None
        self.intsty = None
        self.m = None
        self.wavID = None
        self.phase = None
        
        self.h, self.k, self.l = int(lin.split()[0]), int(lin.split()[1]), int(lin.split()[2])
        self.d = float(lin.split()[3])
        self.fr = float(lin.split()[4])
        self.fi = float(lin.split()[5])
        self.fAbs = float(lin.split()[6])
        self.ang = float(lin.split()[7])
        self.intsty = float(lin.split()[8])
        self.m = int(lin.split()[9])
        self.wavID = int(lin.split()[10])
        self.phase = int(lin.split()[11])
        
    def Print(self):
        print("(h, k, l) = ("+ str(self.h) + ", " + str(self.k) + ", " + str(self.l) + ')')
        print("2-Theta: " + f'{self.ang:.3f}')
        print("Intensity: " + f'{self.intsty:.3f}')
    
#Reads in a file with formatting from VESTA's XRD output
#Returns a list of 'line' objects created from the outfile
def ParseOutfile(infileLoc):
    ret = []
    with open(str(infileLoc), 'r', errors='ignore') as infile:     
        for lineNum, line_ in enumerate(infile):
            if(lineNum > 0): ##don't read in the header
                ret.append(line(line_))
    infile.close()
        
    return ret           
        
        
        
