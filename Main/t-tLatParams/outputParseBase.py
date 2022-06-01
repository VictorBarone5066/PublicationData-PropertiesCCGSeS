class line:
    directory = None
    wDirectory = None
    nrg = None
    nConv = None
    irredKpts = None
    aV, bV, cV = None, None, None
    al, be, ga = None, None, None
    vol = None
    
    def __init__(self, lin):
        #The fact that I have to do this is completly insane.
        self.directory = None
        self.wDirectory = None
        self.nrg = None
        self.nConv = None
        self.irredKpts = None
        self.aV, self.bV, self.cV = None, None, None
        self.al, self.be, self.ga = None, None, None
        self.vol = None
        
        #Now to the actual initialization from csv file
        self.directory = lin.split(',')[0]
        self.wDirectory = lin.split(',')[1]
        self.nrg = float(lin.split(',')[2])
        self.nConv = int(lin.split(',')[3])
        self.irredKpts = int(lin.split(',')[4])
        self.aV, self.bV, self.cV = float(lin.split(',')[5]), float(lin.split(',')[6]), \
                                    float(lin.split(',')[7])
        self.al, self.be, self.ga = float(lin.split(',')[8]), float(lin.split(',')[9]), \
                                    float(lin.split(',')[10])
        self.vol = float(lin.split(',')[11])
        
    def Print(self, directory=True, wDirectory=True, nrg=True, nConv=True, irredKpts=True,
              abc=True, albega=True, vol=True):
        if(directory):
            print("Dir: " + self.directory)
        if(wDirectory):
            print("Work Dir: " + self.wDirectory)
        if(nrg):
            print("Energy: " + f'{self.nrg:.3f}')
        if(nConv):
            print("Num Converged Steps: " + str(self.nConv))
        if(irredKpts):
            print("Irreducible K-points in 1st BZ: " + str(self.irredKpts))
        if(abc):
            print("(a, b, c) = (" + f'{self.aV:.3f}' + ", " + f'{self.bV:.3f}' + ", " + 
              f'{self.cV:.3f}' + ")")
        if(albega):
            print("(alpha, beta, gamma) = (" + f'{self.al:.3f}' + ", " + f'{self.be:.3f}' + ", " + 
              f'{self.ga:.3f}' + ")")
        if(vol):
            print("Volume: " + f'{self.vol:.3f}')
        

#Reads in a very specific file format
#Returns a list of 'line' objects created from the outfile
def ParseOutfile(infileLoc):
    ret = []
    with open(str(infileLoc), 'r') as infile:     
        for lineNum, line_ in enumerate(infile):
            if(lineNum > 0): ##don't read in the header
                ret.append(line(line_))
    infile.close()
        
    return ret    
        
        
        
        
        