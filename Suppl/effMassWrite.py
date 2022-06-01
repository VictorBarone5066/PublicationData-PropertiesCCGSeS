import headerPoscar as POH
import headerBandPlot as PBH
import headerEffMass as EMH
import numpy as np

tetraFolder = "t-tEffMass"
orthoFolder = "o-oEffMass"

EIGENVAL_LOC = {"tetra": tetraFolder + "//" + "EIGENVALx",
              "ortho": orthoFolder + "//" + "EIGENVALx"}
POSCAR_LOC = {"tetra": tetraFolder + "//" + "POSCARx",
               "ortho": orthoFolder + "//" + "POSCARx"}
OUTFILE_LOC = "allEffMassses.csv"

#Fitting params
INTERP_TYPE = "quad" ##spline or quad
PTS_FOR_FIT = 3 ##total number of points including extrema.  Spline ignores this

#Finds the max number of effective mass entries for the multi-layered dictionary defined later
def GetMaxEffMassEntries(dic):
    max = -1
    for a in dic.values():
        for b in a.values():
            for c in b.values():
                ##Now, we're at the lists
                if(len(c) > max):
                    max = len(c)

    return max - 2 ##subtract 2 because there are two lists that are not effective masses

#Initialize a disgusting dictionary that looks something like:
#massData =  {tetraDict = {}, orthoDict = {}}       keys = "tetra", "ortho"
#...Dict =   {holeDict_ = {}, elecDict_ = {}}       keys = "hole", "elec"
#...Dict_ =  {concDict0 = [], concDict1 = [], ...}  keys = 0.000, 0.125, 0.250, 0.375, ...
#concDictx = [[extrema k-point locs in direct recip coords (length 3)],
#             [1-indexed band number of extrema (length 1)],
#             followed by n lists (one for each line), length 4 formatted as [nth effective mass,
#                                                                             *three unit vector components*]
#             ]

#Do the following for each set of tetra and ortho xs
massData = {}
for phaseId in ["tetra", "ortho"]:
    massData[phaseId] = {}
    massData[phaseId]["hole"] = {}
    massData[phaseId]["elec"] = {}

    include = list(np.arange(0.000, 1.125, 0.125))
    include = [(f'{i:.3f}'[0] + f'{i:.3f}'[2:]) for i in include] ##remove decim pt
    for conc in include:
        #Get all k-points, get important band stats
        kps = PBH.ParseEigenval(EIGENVAL_LOC[phaseId] + conc, occupanciesListed=True)

        #Get real-space lattice vectors
        pos = POH.Poscar(POSCAR_LOC[phaseId]+ conc)
        pos.ConvertToCartesian() ##apply univ scal fact
        a1, a2, a3 = pos.superCellVecA, pos.superCellVecB, pos.superCellVecC

        #---Hole data---
        vbm = EMH.GetVBMKpoint(kps)
        bandStats = EMH.GetLowestHoleBand(vbm.bands)
        massData[phaseId]["hole"][float(conc)/1000.] = [[vbm.a, vbm.b, vbm.c], [bandStats[0] - 1]]

        lines = EMH.GetMaximaLines(vbm, kps)
        for num, line in enumerate(lines):
            ##Get actual data
            thisLine = EMH.RecipLine(vbm, line, bandStats[0] - 1, a1, a2, a3) ## - 1 to sample holes at VBM
            effMass = thisLine.GetEffMass("hole", npts=PTS_FOR_FIT, useSymm=True, interpType=INTERP_TYPE)
            uv = thisLine.GetLineDir(acc=9) ##direct unit vector along line starting at extrema
            massData[phaseId]["hole"][float(conc)/1000.].append([effMass, uv[0], uv[1], uv[2]])


        #---Elec data---
        cbm = EMH.GetCBMKpoint(kps)
        bandStats = EMH.GetHighestElectronBand(cbm.bands)
        massData[phaseId]["elec"][float(conc)/1000] = [[cbm.a, cbm.b, cbm.c], [bandStats[0] + 1]]

        lines = EMH.GetMaximaLines(cbm, kps)
        for line in lines:
            #Get actual data
            thisLine = EMH.RecipLine(cbm, line, bandStats[0] + 1, a1, a2, a3) ## + 1 to sample elecs at CBM
            effMass = thisLine.GetEffMass("elec", npts=PTS_FOR_FIT, useSymm=True, interpType=INTERP_TYPE)
            uv = thisLine.GetLineDir(acc=9) ##direct unit vector along line starting at extrema
            massData[phaseId]["elec"][float(conc) / 1000].append([effMass, uv[0], uv[1], uv[2]])

#Now, write everything
with open(OUTFILE_LOC, 'w') as outfile:
    ##Prepare header
    maxEntries = GetMaxEffMassEntries(massData)
    colLabs = "x,extrA,extrB,extrC,extBandNum"
    for i in range(0, maxEntries):
        colLabs += ',' + str(i) + "_m*/m0"
        for j in ['A', 'B', 'C']:
            colLabs += ',' + str(i) + "_uv" + j

    for phaseKey, phaseDict in massData.items():
        outfile.write(str(phaseKey) + "\n")

        for massKey, massDict in phaseDict.items():
            outfile.write(str(massKey) + "\n" + colLabs + "\n")

            for concKey, concList in massDict.items():
                outfile.write(str(concKey))
                for kpInd in concList[0]:
                    outfile.write(',' + str(kpInd))
                for bandNum in concList[1]:
                    outfile.write(',' + str(bandNum))
                for massListItr in range(2, len(concList)):
                    for thisMassData in concList[massListItr]:
                        outfile.write(',' + str(thisMassData))
                outfile.write("\n")

    outfile.close()
