import headerPoscar as hp
import headerEigenval as he
import headerKTransform as hk

which = "x"

PATH_TO_POSCAR = "POSCAR" + which
PATH_TO_EIGENVAL = "EIGENVAL" + which
PATH_TO_OUTFILE = "eigOrtho-" + which + ".aem"

#Get reciprocal lattice
poscar = hp.Poscar(PATH_TO_POSCAR)
A = [poscar.superCellVecA,
     poscar.superCellVecB,
     poscar.superCellVecC]
B = hk.RealToRecip(A)

m = max([(B[i][0]*B[i][0] + B[i][1]*B[i][1] + B[i][2]*B[i][2])**(1./2.) for i in range(0, 3)])
m = 1.05*m
print(B[0], "\n", B[1], "\n", B[2], "\n\n", m)


#Transform k from direct to cartesian
origKps = he.ParseEigenval(PATH_TO_EIGENVAL, occupanciesListed=True)#[:400]
cartKps = hk.TransformKps(origKps=origKps, B=B, isym=0, cubeLen=m)

#Write new formatting
with open(PATH_TO_OUTFILE, 'w') as outfile:
     outfile.write(str(len(cartKps)) + ' ' + str(len(cartKps[0].bands)) + "\n")
     for i in range(0, len(cartKps)):
          outfile.write(str(cartKps[i].a) + ' ' + str(cartKps[i].b) + ' ' + str(cartKps[i].c) + "\n")
          for j in range(0, len(cartKps[i].bands)):
               outfile.write(str(cartKps[i].bands[j][1]) + "\n")
     outfile.write("STOP")
     outfile.close()
