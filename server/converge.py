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


if len(ircList) == 1:
   ircList1=ircList
else:
   ircList1=ircList
   ircList1[0] = ircList[1]
   ircList1[0][0] = float(ircList[1][0]) + 2
   ircList1[0][1] = float(ircList[1][2]) + 0.001
   ircList1[0][3] = 'F'
   coordFile = open(sys.argv[1] + '.csv', "w+")
   for Index, Molecule in enumerate(ircList1[1:]):
      coordFile.write(str(Index)+ ',' + ",".join(str(bit) for bit in Molecule))
      coordFile.write('\n')
      if Molecule[3]=="T" and Molecule[6]=="T" and Molecule[9]=="T" and Molecule[12]=="T" and Molecule[15]=="T":
         break
   coordFile.close()



# coordFile = open(sys.argv[1] + '.csv', "w+")
# for Index, Molecule in enumerate(ircList[0:]):
#    coordFile.write(str(Index)+ ',' + ",".join(str(bit) for bit in Molecule))
#    coordFile.write('\n')
#    if Molecule[3]=="T" and Molecule[6]=="T" and Molecule[9]=="T" and Molecule[12]=="T" and Molecule[15]=="T":
#       break
# coordFile.close()
