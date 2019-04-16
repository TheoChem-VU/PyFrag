import os, string, re, sys

ircFile = open(sys.argv[1], 'r')
ircRaw  = []
for line in ircFile:
   llist = line.split()
   ircRaw.append(llist)
ircRawList = [_f for _f in ircRaw if _f]
ircFile.close()

currentPath = os.getcwd()
fileName  = os.path.join(currentPath, str('PYFRAG.csv'))
coordFile = open(str(fileName), "w")
for Index, Molecule in enumerate(ircRawList):
   for atom in Molecule[0:-1]:
      coordFile.write(atom + ',')
   coordFile.write(Molecule[-1])
   coordFile.write('\n')
coordFile.close()

# csvName=sys.argv[2]

# tableFile = open(csvName, "a")
# for Index, Molecule in enumerate(ircRawList):
#     tableFile.write('"","' + "  ".join(str(bit) for bit in Molecule) + '"')
#     tableFile.write('\n')
# tableFile.close()

