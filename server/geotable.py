import os, sys
from os.path import dirname, basename

csvName=sys.argv[3]
jobName=basename(sys.argv[1])

ircFile = open(sys.argv[1], 'r')
ircList  = []
for line in ircFile:
   llist = line.split()
   lllen = len(llist)
   if lllen != 0:
     ircList.append(llist)
ircFile.close()


conFile = open(sys.argv[2], 'r')
conList  = []
for line in conFile:
   llist = line.split()
   lllen = len(llist)
   if lllen != 0:
     conList.append(llist)
conFile.close()


coordFile = open(csvName, "a")
coordFile.write('"","' + jobName + ',' + '   -' + conList[0][0] + '"')
coordFile.write('\n')
for Index, Molecule in enumerate(ircList[0:]):
    coordFile.write('"","' + "  ".join(str(bit) for bit in Molecule) + '"')
    coordFile.write('\n')
coordFile.close()
