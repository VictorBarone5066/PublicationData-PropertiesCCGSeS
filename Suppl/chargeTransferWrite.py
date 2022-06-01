import numpy as np

import headerPoscar as P
import headerChargeTransfer as C

tetraFolder = "t-tChargeTransfer"
orthoFolder = "o-oChargeTransfer"

VALENCE_LOC = {"tetra": tetraFolder + "//" + "nElect.txt",
               "ortho": orthoFolder + "//" + "nElect.txt"}
POSCAR_LOC = {"tetra": tetraFolder + "//" + "CONTCARx",
              "ortho": orthoFolder + "//" + "CONTCARx"}
ACF_LOC = {"tetra": tetraFolder + "//" + "acfx",
           "ortho": orthoFolder + "//" + "acfx"}
OUTFILE_LOC = "allChgTransfer.csv"

with open(OUTFILE_LOC, 'w') as outfile:    
    #Do the following for each set of tetra and ortho xs
    for phaseId in ["tetra", "ortho"]:
        include = list(np.arange(0.000, 1.125, 0.125))
        include = [(f'{i:.3f}'[0] + f'{i:.3f}'[2:]) for i in include] ##remove decim pt
        #Read in each ACF.dat, CONTCAR pair
        pairs = []
        for conc in include:
            pairs.append([conc, C.model(P.Poscar(POSCAR_LOC[phaseId] + conc), 
                                        C.ParseValence(VALENCE_LOC[phaseId]), 
                                        C.ParseAcf(ACF_LOC[phaseId] + conc))])
        
        #Write outfile.  Format:  conc, site 1 elem, site 1 actChg, site 2 elem, site 2 actChg, ...
        #Header
        outfile.write(phaseId + "\n")
        head = "conc*1000"
        for i in range(0, len(pairs[0][1].poscar.atoms)):
            head = head + ",elem" + str(i + 1) + ",chg" + str(i + 1) 
        outfile.write(head + "\n")
        #The rest
        for pair in pairs:
            line = pair[0] ##the concentration
            for i in range(0, len(pair[1].actChgLis)):
                line += ',' + pair[1].poscar.atoms[i].atomType + ',' + str(pair[1].actChgLis[i])    
            outfile.write(line + "\n")
        
    outfile.close()







    