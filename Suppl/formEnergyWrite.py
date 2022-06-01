tetraFolder = "t-tFormationEnergies"
orthoFolder = "o-oFormationEnergies"

BULK_ENERGY_LOC = {"tetra": tetraFolder + "//" + "bulkEnergies.csv",
                   "ortho": orthoFolder + "//" + "bulkEnergies.csv"}
MATERIAL_ENERGY_LOC = {"tetra": tetraFolder + "//" + "convOutput.csv",
                       "ortho": orthoFolder + "//" + "convOutput.csv"}
OUTFILE_LOC = "allFormEnergies.csv"

def AtomCountsFromConc(x):
    ret = {"Cd": None,
           "Cu": None,
           "Ge": None,
           "Se": None,
           'S': None}

    #Hardcoded because I need this done fast.  Later, read in CONTCARs.  Or not bc I'll only need
    #this once
    ret["Cd"] = 2
    ret["Cu"] = 4
    ret["Ge"] = 2
    ret['S'] = int(8.*x)
    ret["Se"] = int(8.*(1. - x))

    return ret

with open(OUTFILE_LOC, 'w') as outfile:
    for phaseId in ["tetra", "ortho"]:
        bulkNrgPerAtom = {"Cd": None,
                          "Cu": None,
                          "Ge": None,
                          "Se": None,
                          'S': None}

        #Get the energy per atom from the bulk calculations
        with open(BULK_ENERGY_LOC[phaseId], 'r') as infile:
            #infile format: elem, total energy, number of atoms in cell
            for num, line in enumerate(infile):
                if(num > 0): ##skip header
                    cont = line.split(',')
                    bulkNrgPerAtom[cont[0]] = float(cont[1])/float(cont[2])
            infile.close()

        #Read total energies from converged output
        with open(MATERIAL_ENERGY_LOC[phaseId], 'r') as infile:
            #At the same time, write results to output
            outfile.write(phaseId + "\n")
            outfile.write("conc,eForm\n")
            for num, line in enumerate(infile):
                ##skip header/useless concs
                if(num > 0 and line.split(',')[1] != "x04375" and
                   line.split(',')[1] != "x05625"):
                    thisConc = float(line.split(',')[1][1:])/1000.
                    thisNrg = float(line.split(',')[2])
                    ##Form energy is (energy of material) - (energy of parts)
                    bulkNrg = 0.
                    for atom, count in AtomCountsFromConc(thisConc).items():
                        bulkNrg += bulkNrgPerAtom[atom]*count
                    outfile.write(str(thisConc) + ',' + str(thisNrg - bulkNrg) + "\n")
            infile.close()
    outfile.close()
