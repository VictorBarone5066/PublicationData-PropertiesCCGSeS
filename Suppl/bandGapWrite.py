import numpy as np
import headerDensityOfStates as dos

tetraFolder = "t-tBandgap"
orthoFolder = "o-oBandgap"
DOSCAR_LOC =  {"tetra": tetraFolder + "//" + "DOSCARx",
               "ortho": orthoFolder + "//" + "DOSCARx"}
BG_DATA_SAVE_LOC = "allBandgap.csv"

TOL = 0.100 #in whatever units the DOS is written as

def FermiLoc(nrgs, ef):
    bst, loc = abs(nrgs[0] - ef), 0
    for i in range(0, len(nrgs)):
        if(abs(nrgs[i] - ef) < bst):
            bst = abs(nrgs[i] - ef)
            loc = i
    return loc

data = {} #data[phase] = [[concs], [bandgaps], [line details]]
with open(BG_DATA_SAVE_LOC, 'w') as outfile:
    for phaseId in ["tetra", "ortho"]:
        outfile.write(phaseId + "\n" + "conc,bg,eF\n")

        include = list(np.arange(0.000, 1.125, 0.125))
        include = [(f'{i:.3f}'[0] + f'{i:.3f}'[2:]) for i in include] ##gets rid of the decimal point
        for num, conc in enumerate(include):
            #Find fermi level / its location, get energies and dos values in a more workable format
            eFermi = dos.GetEnergyInfo(DOSCAR_LOC[phaseId] + conc, 0)["eFermi"]
            vals = dos.GetAtomDosInfo(DOSCAR_LOC[phaseId] + conc, 0, spin = False)
            energies = vals["energy"]
            eFermiLoc = FermiLoc(energies, eFermi)
            doss = vals["dos"]

            #Verify that there are no states at eFermi
            if(doss[eFermiLoc] >= TOL):
                print("State(s) at the fermi level (" + str(eFermi) + ") eV!")
                exit

            #Get bandgap value
            ##Search backwards from eFermi
            i = eFermiLoc
            while(i > 0):
                if(doss[i] >= TOL):
                    vbm = 1./2. * (energies[i + 1] + energies[i])
                    break
                i = i - 1
            ##Serch forwards from eFermi
            i = eFermiLoc
            while(i < len(energies)):
                if(doss[i] >= TOL):
                    cbm = 1./2. * (energies[i] + energies[i - 1])
                    break
                i = i + 1

            #Write data
            outfile.write(str(float(conc)/1000) + ',' + str(cbm - vbm) + ',' + str(eFermi) + "\n")
    outfile.close()
