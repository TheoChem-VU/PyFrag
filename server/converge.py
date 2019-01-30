import os, string, re, sys
from os.path import dirname, basename

ircFile = open(sys.argv[1], 'r')
ircList  = []
for line in ircFile:
   llist = line.split()
   lllen = len(llist)
   if lllen != 0 and lllen != 11:
      ircList.append(llist)
ircFile.close()


coordFile = open(sys.argv[1] + '.csv', "w+")
for Index, Molecule in enumerate(ircList[0:]):
   coordFile.write(str(Index)+ ',' + ",".join(str(bit) for bit in Molecule))
   coordFile.write('\n')
   if Molecule[3]=="T" and Molecule[6]=="T" and Molecule[9]=="T" and Molecule[12]=="T" and Molecule[15]=="T":
      break
coordFile.close()


# coordFile = open(sys.argv[1] + '.csv', "w+")
# for Index, Molecule in enumerate(ircList[0:]):
#     coordFile.write(str(Index)+ ',' + ",".join(str(bit) for bit in Molecule))
#     coordFile.write('\n')
# coordFile.close()

# baseName=basename(sys.argv[1])
# dirName=dirname(sys.argv[1])

# for index, item in enumerate(['ener', 'enech', 'gmaxs', 'grmss', 'smaxs', 'srmss']):
#   coordFile = open(dirName + "/table_" + baseName + item + '.csv', "w+")
#   for Index, Molecule in enumerate(ircList[0:-1]):
#       coordFile.write(str(Index) + ',0,0,0,0,' + Molecule[index-1] + ',0')
#       coordFile.write('\n')
#   coordFile.write(str(len(ircList)-1) + ',0,0,0,0,' + ircList[-1][0] + ',0')
#   coordFile.close()
