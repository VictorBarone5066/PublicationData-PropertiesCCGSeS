OUTFILE_LOC = "baderConvTest.csv"

#default [NGXF, NGYF, NGZF]
DEFAULTS = [80., 80., 96.]

#Gives a list of multiplicitave factors that transform current NG*F to default NG*F.
#i.e. if the name is acf_160_240_144, and defaults are [80, 80, 96], this will return [2, 3, 1.5]
def ParseName(s):
    ret = []
    for i in range(0, len(DEFAULTS)):
        ret.append(str(float(s.split('_')[i + 1]) / DEFAULTS[i]))
    return ret

toOpen = [[80., 80., 96.], 
          [120., 120., 144.],
          [160., 160., 192.],
          [200., 200., 240.],
          [240., 240., 288.],
          [280., 280., 336.]]

with open(OUTFILE_LOC, 'w') as outfile:
    #Setup outfile loc
    head = ""
    for i in range(1, 17):
        head += ',' + str(i)
    head += "\n"
    outfile.write("mx,my,mz" + str(head))
    
    #Read each infile
    for nameNums in toOpen:
        name = "acf_" + str(int(nameNums[0])) + '_' + str(int(nameNums[1])) + '_' + \
                str(int(nameNums[2]))
        with open(name, 'r') as infile:
            hyphSum = 0 #at the first hyphen, begin reading.  At the second, end reading
            dat = [None]*16          
            for line in infile:
                if(line[5] == '-' and line[6] == '-'):
                    hyphSum = hyphSum + 1
                    continue
                if(hyphSum == 1):
                    dat[int(line.split()[0]) - 1] = line.split()[4] 
                if(hyphSum == 2):
                    break
            infile.close()
                
        #Add each to the outfile
        p = ParseName(name)
        outline = p[0] + ','+ p[1] + ',' + p[2]
        for d in dat:
            outline += ',' + d
        outline += "\n"
        outfile.write(outline)
    
    outfile.close()
                
                
                
                
                
                
    outfile.close()