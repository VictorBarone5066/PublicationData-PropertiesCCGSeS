import numpy as np
import headerElas as h

tetraFolder = "t-tElas"
orthoFolder = "o-oElas"

STRESSES_LOC = {"tetra": tetraFolder + "//" + "outputstresses.csv",
                "ortho": orthoFolder + "//" + "outputstresses.csv"}
CIJSIJ_LOC = {"tetra": tetraFolder + "//" + "tetraCijsSijs.csv",
               "ortho" : orthoFolder + "//" + "orthoCijsSijs.csv"}

def GetCij(dataList, i, j, c11=None, c22=None, c33=None):
    #Initial paramaters needed for (1/v0)d^2E/de^2      
    deltas, energies, v0 = [], [], -1 #-1 for error checking
    for dat in dataList:
        if(dat.delta == 0):
            v0 = dat.vol
        deltas.append(dat.delta)
        energies.append(dat.energy)

    #Error checking
    if(v0 == -1):
        print("Volume not initialized.  Could not find delta=0")
        return -1
    ##Need c11, c22, c33 for some elastic constants
    if((i == 1 and j == 2) or (i == 2 and j == 1)):
        if(c11 == None and c22 == None):
            print("C" + str(i) + str(j) + ": c11 and c22 not passed")
    if((i == 1 and j == 3) or (i == 3 and j == 1)):
        if(c11 == None and c33 == None):
            print("C" + str(i) + str(j) + ": c11 and c33 not passed")
    if((i == 2 and j == 3) or (i == 3 and j == 2)):
        if(c22 == None and c33 == None):
            print("C" + str(i) + str(j) + ": c22 and c33 not passed")      
    
    #Elastic constant calculations
    secDeriv = 2.0*h.GetQuadAParam(deltas, energies)
    ##c11, c22, c33: non-volume conserving
    if(i == j and i <= 3):
        return (1.0/v0)*secDeriv
    ##c44, c55, c66: shear
    if(i == j and i <=6 and i>= 4):
        return (1.0/(4.0*v0)*secDeriv)
    ##c12, c13, c23: remaining
    if(i == 1 and j == 2):
        return (1.0/2.0)*(c11 + c22 - (1.0/v0)*secDeriv)
    if(i == 1 and j == 3):
        return (1.0/2.0)*(c11 + c33 - (1.0/v0)*secDeriv)
    if(i == 2 and j == 3):
        return (1.0/2.0)*(c22 + c33 - (1.0/v0)*secDeriv)    

#dict[id] = [i, j]
deformIdMap = {1: [1, 1],
               2: [2, 2],
               3: [3, 3],
               4: [4, 4],
               5: [5, 5],
               6: [6, 6],
               7: [1, 2],
               8: [1, 3],
               9: [2, 3]}


for phaseId in ["tetra", "ortho"]:
    data = h.ParseOutfile(STRESSES_LOC[phaseId])
    
    #Write outfile
    with open(CIJSIJ_LOC[phaseId], 'w') as outfile:
        outfile.write("x,c11,c22,c33,c44,c55,c66,c12,c13,c23,s11,s22,s33,s44,s55,s66,s12,s13,s23\n")
    
        for x, xDefIdDict in data.items():
            outline = str(x)
            
            cijs = [] ###makes it easier to ref cijs. size 7 instead of 6 bc of 0-indexing
            for i in range(0, 7):
                cijs.append([])
                for j in range(0, 7):
                    cijs[i].append(0.0)
            #note cijs=[[None]*7]*7 SHOULD WORK, but this creates a list of 7 arrays that are all 
            #referencing one another - changing one affects the others.  This is absolutely retarded
            
            for iden, dataList in xDefIdDict.items():
                i, j = deformIdMap[iden][0], deformIdMap[iden][1]
                if(i == j):
                    cijs[i][i] = GetCij(dataList[:], i, j)            
                if(i == 1 and j == 2):
                    cijs[i][j] = GetCij(dataList[:], i, j, c11=cijs[1][1], c22=cijs[2][2])
                    cijs[j][i] = GetCij(dataList[:], i, j, c11=cijs[1][1], c22=cijs[2][2])
                if(i == 1 and j == 3):
                    cijs[i][j] = GetCij(dataList[:], i, j, c11=cijs[1][1], c33=cijs[3][3])   
                    cijs[j][i] = GetCij(dataList[:], i, j, c11=cijs[1][1], c33=cijs[3][3])
                if(i == 2 and j == 3):
                    cijs[i][j] = GetCij(dataList[:], i, j, c22=cijs[2][2], c33=cijs[3][3])
                    cijs[j][i] = GetCij(dataList[:], i, j, c22=cijs[2][2], c33=cijs[3][3])
            
            ##Output stability according to Born criteria:
            ##From Phys. Rev. B90, 224104 (2014)
            c = cijs
            print("\n\n----- CONC " + str(x) + " -----")
            #Orthorhombic
            print("***ORTHORHOMBIC***")
            print("c11 > 0?\t\t\t" + str(c[1][1] > 0.0))
            print("c44 > 0?\t\t\t" + str(c[4][4] > 0.0))
            print("c55 > 0?\t\t\t" + str(c[5][5] > 0.0))
            print("c66 > 0?\t\t\t" + str(c[6][6] > 0.0))
            print("c11*c22 > c12^2?\t\t" + str(c[1][1]*c[2][2] > (c[1][2])**(2.))) 
            un = c[1][1]*c[2][2]*c[3][3]
            do = 2*c[1][2]*c[1][3]*c[2][3]
            tr = c[1][1]*(c[2][3])**(2.)
            qu = c[2][2]*(c[1][3])**(2.)
            si = c[3][3]*(c[1][2])**(2.)
            print("cii + 2*cij - cii*cij^2?\t" + str(un + do - tr - qu - si > 0.0))
            
            ##Get into GPa units
            for i_ in range(0, len(cijs)):
                for j_ in range(0, len(cijs[i_])):
                    cijs[i_][j_] = h.EvAngstToGPa(cijs[i_][j_])
    
            ##Get Sij by taking the inverse of the Cij matrix
            Cij = []
            for i_ in range(0, 6):
                Cij.append([])
                for j_ in range(0, 6):
                    Cij[i_].append(0.0)
            for i_ in range(0, len(cijs) - 1):
                for j_ in range(0, len(cijs[i_]) - 1):
                    Cij[i_][j_] = cijs[i_+1][j_+1]
            Sij = np.linalg.inv(Cij)
            Sij = np.ndarray.tolist(Sij) ##bc numpy arrays are annoying
            ##To stay consistent with my other formatting, i want sijs[1][1] to equal S11 (not S22)
            sijs = [[0.0]*6] + Sij[:] 
            for i_ in range(0, len(sijs)):
                sijs[i_] = [0.0] + sijs[i_]
            
            ##Add cijs to the output string
            for i_ in range(1, 7):
                outline += ',' + str(cijs[i_][i_])
            outline += ',' + str(cijs[1][2])
            outline += ',' + str(cijs[1][3])
            outline += ',' + str(cijs[2][3])
            ##Add sijs to the output string
            for i_ in range(1, 7):
                outline += ',' + str(sijs[i_][i_])
            outline += ',' + str(sijs[1][2])
            outline += ',' + str(sijs[1][3])
            outline += ',' + str(sijs[2][3])    
            outfile.write(outline + "\n")
        outfile.close()
    

    
    
    
    
    