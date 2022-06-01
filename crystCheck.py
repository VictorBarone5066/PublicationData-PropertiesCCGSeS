import pymatgen as mg
import pymatgen.core.structure as struc
import pymatgen.symmetry.analyzer as sa

def Exit():
    exit(1)

START = "C://Users//baron//Desktop//Cu2CdGe(SxSe1-x)4//"

thisSymprec = 1.
for i in range(0, 10):
    thisSymprec /= 10.
    for sym in ["tetra-tetra", "ortho-ortho"]:
        posLoc = START + "//old//" + sym + "//POSCARx"
        conLoc = START + sym + "//convergedModels//CONTCARx"

        #print("\n===***===" + sym + "===***===")
        fileDict = {}
        for fi in [posLoc, conLoc]:
            #print("===" + fi.split('/')[-1][:-1] + "===")
            for nu, x in enumerate(["0000", "0125", "0250", "0375", "0500", "0625", "0750", "0875", "1000"]):
                s = struc.Structure.from_file(fi + x)
                a = sa.SpacegroupAnalyzer(s, symprec=thisSymprec, angle_tolerance=5)
                n = a.get_space_group_number()
                k = a.get_crystal_system()
                r = a.get_space_group_symbol()

                fileDict[str(nu) + fi.split('/')[-1][:-1]] = [n, k, r]
                #print(str(nu) + ": " + str(n) + ' ' + str(k) + ' ' + str(r) + "\n")

        #Check to see if n, k, r for POSCAR matches n, k, r for the corresponding CONTCAR
        for i in range(0, len(fileDict.keys())):
            for j in range(i+1, len(fileDict.keys())):
                if(list(fileDict.keys())[i][0] == list(fileDict.keys())[j][0]):
                    #print(list(fileDict.keys())[i], list(fileDict.keys())[j])
                    #print(list(fileDict.values())[i], list(fileDict.values())[j])
                    if(list(fileDict.values())[i] == list(fileDict.values())[j]):
                        #print("Match")
                        continue
                    else:
                        print("Failed at symprec= " + str(thisSymprec))
                        Exit()
