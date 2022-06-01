import Poscar as h

POSCARS_LOC = "old//tetra-tetra//"
CONTCARS_LOC = "tetra-tetra//convergedModels//"

OUT_LOC = "C://Users//baron//Desktop//devFrac.csv"

DIFF = 0.75 #fractional

def Dist(st, fi):
    return ((fi.a-st.a)**2. + (fi.b-st.b)**(2.) + (fi.c-st.c)**(2.))**(1./2.)

with open(OUT_LOC, 'w') as outfile:
    outfile.write("Deviation from POSCAR in Angstroms\n")
    outfile.write("conc*1000,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,avg\n")

    for x in ["0000", "0125", "0250", "0375", "0500", "0625", "0750", "0875", "1000"]:
        pos = h.Poscar(POSCARS_LOC + "POSCARx" + x)
        con = h.Poscar(CONTCARS_LOC + "CONTCARx" + x)

        pos.ConvertToDirect()
        con.ConvertToDirect()

        outfile.write(str(x))

        avg = 0.
        for p, c in zip(pos.atoms, con.atoms):
            #a:
            if(c.a > p.a + DIFF):
                c.a = c.a - 1.0
            elif( p.a > c.a + DIFF):
                p.a = p.a - 1.0
            #b
            if(c.b > p.b + DIFF):
                c.b = c.b - 1.0
            elif( p.b > c.b + DIFF):
                p.b = p.b - 1.0
            #c
            if(c.c > p.c + DIFF):
                c.c = c.c - 1.0
            elif( p.c > c.c + DIFF):
                p.c = p.c - 1.0

#            pos.ConvertToCartesian()
#            con.ConvertToCartesian()

            avg += Dist(p, c)
            outfile.write(',' + str(Dist(p, c)))

#            pos.ConvertToDirect()
#            con.ConvertToDirect()
        outfile.write(',' + str(avg/16.) + "\n")

    outfile.close()
