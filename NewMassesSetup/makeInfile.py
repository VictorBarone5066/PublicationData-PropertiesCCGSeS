import headerPoscar as hp
import headerEigenval as he
import headerKTransform as hk

def ff(d, dec=3):
    return f'{d:.{dec}f}'

which = "x"

PATH_TO_POSCAR = "POSCAR" + which
PATH_TO_EIGENVAL = "EIGENVAL" + which
PATH_TO_OUTFILE = "inOrtho-" + which + ".aem"


KPTS = [6, 6, 5] ##a, b, c



#Get reciprocal lattice
poscar = hp.Poscar(PATH_TO_POSCAR)
A = [poscar.superCellVecA,
     poscar.superCellVecB,
     poscar.superCellVecC]
B = hk.RealToRecip(A)

b1 = (B[0][0]*B[0][0] + B[0][1]*B[0][1] + B[0][2]*B[0][2])**(1./2.) * (1.05)
b2 = (B[1][0]*B[1][0] + B[1][1]*B[1][1] + B[1][2]*B[1][2])**(1./2.) * (1.05)
b3 = (B[2][0]*B[2][0] + B[2][1]*B[2][1] + B[2][2]*B[2][2])**(1./2.) * (1.05)

kpLis = he.ParseEigenval(PATH_TO_EIGENVAL, occupanciesListed=True)
vbmkp = he.GetVBMKpoint(kpLis, spin=False)
cbmkp = he.GetCBMKpoint(kpLis, spin=False)
vbmBand = he.GetHighestOccupiedBand(vbmkp, spin=False)
cbmBand = he.GetLowestUnoccupiedBand(cbmkp, spin=False)
vbmInd = vbmBand[0] - 1
vbm = vbmBand[1]
cbmInd = cbmBand[0] - 1
cbm = cbmBand[1]


with open(PATH_TO_OUTFILE, 'w') as outfile:
     outfile.write("BEGIN DOMAIN\n")
     outfile.write(ff(-b1/2.) + " " + ff(+b1/2.) + "\n")
     outfile.write(ff(-b2/2.) + " " + ff(+b2/2.) + "\n")
     outfile.write(ff(-b3/2.) + " " + ff(+b3/2.) + "\n")
     outfile.write("END DOMAIN\n\n")

     outfile.write("BEGIN GRIDPOINTS\n91\n91\n91\nEND GRIDPOINTS\n\n")

     outfile.write("BEGIN LATTICE\n")
     outfile.write(ff(B[0][0], dec=8) + " " + ff(B[0][1], dec=8) + " " + ff(B[0][2], dec=8) + "\n")
     outfile.write(ff(B[1][0], dec=8) + " " + ff(B[1][1], dec=8) + " " + ff(B[1][2], dec=8) + "\n")
     outfile.write(ff(B[2][0], dec=8) + " " + ff(B[2][1], dec=8) + " " + ff(B[2][2], dec=8) + "\n")
     outfile.write("END LATTICE\n\n")

     outfile.write("BEGIN IDWPARAMS\n")
     outfile.write(ff((b1/KPTS[0]*b1/KPTS[0] + b2/KPTS[1]*b2/KPTS[1] + b3/KPTS[2]*b3/KPTS[2])**(1./2.) * 2,
                      dec=3) + "\n")
     outfile.write("3\n")
     outfile.write("END IDWPARAMS\n\n")

     outfile.write("BEGIN SMEARPARAMS\n")
     outfile.write(ff(91/KPTS[0]) + "\n")
     outfile.write(ff(91/KPTS[1]) + "\n")
     outfile.write(ff(91/KPTS[2]) + "\n")
     outfile.write("END SMEARPARAMS\n\n")

     outfile.write("BEGIN VTEMPS\nexplicit\n1\n300\nEND VTEMPS\n\n")
     outfile.write("BEGIN CTEMPS\nexplicit\n1\n300\nEND CTEMPS\n\n")

     outfile.write("BEGIN VMUS\nexplicit\n1\n")
     outfile.write(ff(vbm + (cbm - vbm)/100*2, dec=6) + "\n") ##2% above the VBM
     outfile.write("END VMUS\n\n")

     outfile.write("BEGIN CMUS\nexplicit\n1\n")
     outfile.write(ff(cbm - (cbm - vbm)/100*2, dec=6) + "\n") ##2% below the CBM
     outfile.write("END CMUS\n\n")

     outfile.write("BEGIN VBANDS\nrange\n")
     outfile.write(str(0) + " " + str(vbmInd) + "\n")
     outfile.write("END VBANDS\n\n")

     outfile.write("BEGIN CBANDS\nrange\n")
     outfile.write(str(cbmInd) + " " + str(len(kpLis[0].bands) - 1) + "\n")
     outfile.write("END CBANDS\n\n")

     outfile.write("STOP\n")
