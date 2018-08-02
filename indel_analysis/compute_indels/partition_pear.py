import io, os, sys, shutil
from Bio import SeqIO

if len(sys.argv) != 3:
    print('Usage: parition_pear.py <pearfile> <nump>')
else:
    pearfile = sys.argv[1]
    nump = eval(sys.argv[2])
    dirname = '/'.join(pearfile.split('/')[:-1])
    pearsuffix = pearfile.split('/')[-1]

    for filename in os.listdir(dirname):
        if filename[:len(pearfile)+1] == ('%s_' % pearsuffix):
            os.remove(dirname + '/' + filename)

    fouts = []
    for i in range(nump):
        fouts.append(io.open(pearfile[:-5] + '_%d.' % i + pearfile[-5:],'w'))

    i, j = 0,0
    f = io.open(pearfile)
    for line in f:
        fouts[i%nump].write(line)
        j += 1
        i += (j%4 == 0)
    f.close()
    for i in range(nump):
        fouts[i].close()		

    
