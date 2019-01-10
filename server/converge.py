import os, string, re, sys
from os.path import dirname, basename
from shutil import copyfile
ircFile = open(sys.argv[1], 'r')
ircList  = []
for line in ircFile:
   llist = line.split()
   lllen = len(llist)
   if lllen != 0:
     ircList.append(llist)
ircFile.close()

baseName=basename(sys.argv[1])
dirName=dirname(sys.argv[1])


for index, item in enumerate(['ener', 'enech', 'gmaxs', 'grmss', 'smaxs', 'srmss']):
  coordFile = open(dirName + "/table_" + baseName + item + '.csv', "w+")
  for Index, Molecule in enumerate(ircList[0:-1]):
      coordFile.write(str(Index) + ',0,0,0,0,' + Molecule[index-1] + ',0')
      coordFile.write('\n')
  coordFile.write(str(len(ircList)-1) + ',0,0,0,0,' + ircList[-1][0] + ',0')
  coordFile.close()
