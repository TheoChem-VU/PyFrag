import os, string, re, sys

ircFile = open(sys.argv[1], 'r')
ircRaw  = []
for line in ircFile:
   llist = line.split()
   ircRaw.append(llist)
ircRawList = [_f for _f in ircRaw if _f]
print (ircRawList)
ircFile.close()

currentPath = os.getcwd()
fileName  = os.path.join(currentPath, str('PYFRAG.csv'))
# fileName  = '/Users/xiaobo/Sites/node/bokeh/examples/app/stocks/PYFRAG.csv'
coordFile = open(str(fileName), "w")
for Index, Molecule in enumerate(ircRawList):
   for atom in Molecule[0:-1]:
      coordFile.write(atom + ',')
   coordFile.write(Molecule[-1])
   coordFile.write('\n')
coordFile.close()
