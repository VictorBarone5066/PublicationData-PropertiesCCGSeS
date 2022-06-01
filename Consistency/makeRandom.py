import os
import random
import headerPoscar as hp

home=os.getcwd()

def numSe(x):
    x_ = float(x[0] + '.' + x[1:])
    if(x_ == 0.0):
        return 8
    if(x_ == 0.125):
        return 7
    if(x_ == 0.25):
        return 6
    if(x_ == 0.375):
        return 5
    if(x_ == 0.5):
        return 4
    if(x_ == 0.625):
        return 3
    if(x_ == 0.750):
        return 2
    if(x_ == 0.875):
        return 1
    if(x_ == 1.0):
        return 0


for phase in ['tetra', 'ortho']:
    os.chdir(phase)
    for x in ['0000', '0125', '0250', '0375', '0500', '0625', '0750', '0875', '1000']:
        basePoscar = hp.Poscar("POSCARx" + x)
        os.makedirs(x)
        os.chdir(x)
        for randId in ['1', '2', '3']:
            os.mkdir(randId)
            os.chdir(randId)
            theseInds = [8, 9, 10, 11, 12, 13, 14, 15]
            random.shuffle(theseInds)
            theseInds = theseInds[:numSe(x)]
            newPoscar = hp.deepcopy(basePoscar)
            for n in range(0, len(newPoscar.atoms)):
                if(n in theseInds):
                    newPoscar.atoms[n].atomType = "Se"
                    continue
                elif(n >= 8):
                    newPoscar.atoms[n].atomType = 'S'
                    continue
            newPoscar.Refresh()
            if(x == '0000'):
                newPoscar.ChangeAtomOrder(['Cd', 'Cu', 'Ge', 'Se'])
            elif(x == '1000'):
                newPoscar.ChangeAtomOrder(['Cd', 'Cu', 'Ge', 'S'])
            else:
                newPoscar.ChangeAtomOrder(['Cd', 'Cu', 'Ge', 'Se', 'S'])               
            newPoscar.Write("POSCAR")

            os.chdir("..//")
        os.chdir("..//")
    os.chdir("..//")
