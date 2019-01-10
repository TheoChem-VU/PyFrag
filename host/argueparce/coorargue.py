import os, string, re, sys
from shutil import copyfile

def defragment(file_path):
   with open(file_path, 'r') as f:
        lines = len(list(filter(lambda x: x.strip(), f)))
   atomNums = list(range(1, int(lines)+1))
   fragment = ' '.join(str(atomNum) for atomNum in atomNums)
   return fragment

# path = "/Users/xiaobo/Desktop/1.xyz"

ircFile = open(sys.argv[1], 'r')
ircRaw  = [[]]
for line in ircFile:
   llist = line.split()
   lllen = len(llist)
   # if lllen == 4:
   if lllen != 0:
      # append coordinate
      ircRaw[-1].append(llist)
   else:
      # initiate new coordinate block
      ircRaw.append([])
ircRawList = [_f for _f in ircRaw if _f]
ircFile.close()

currentPath = os.getcwd()
for Index, Molecule in enumerate(ircRawList):
   fileDir  = os.path.join(currentPath, str(Index+1))
   filepath = os.path.join(fileDir, str(Index+1) + '.xyz')
   os.mkdir(fileDir)
   os.chdir(fileDir)

   coordFile = open(str(Index+1) + '.xyz', "w")
   for atom in Molecule[1:]:
      for coordinate in atom:
         coordFile.write(coordinate + '   ')
      coordFile.write('\n')
   coordFile.close()

   infoFile  = open(str(Index+1)  + '.txt', "w")
   for term in Molecule[0]:
      infoFile.write(term + '   ')
   infoFile.close()

   if Molecule[0][0] == 'TS:':
      tsFile  = open(str(Index+1)  + '.ts', "w")
   #change to -8 [0, 4] or -4 for 2 or 1 constrain
      tsLT = [Molecule[0][-12:][i:i+4] for i in [0,4,8]]
      for constaint in tsLT:
         for step in constaint:
            tsFile.write(step + '   ')
         tsFile.write('\n')
      tsFile.close()
      print ('--tspath ', filepath, ' \\')
      print ('--tsfragment ',defragment(filepath),' \\')

   elif Molecule[0][0] == 'R1:':
      print ('--r1path ', filepath, ' \\')
      print ('--r1fragment ',defragment(filepath),' \\')

   elif Molecule[0][0] == 'R2:':
      print ('--r2path ', filepath, ' \\')
      print ('--r2fragment ',defragment(filepath),' \\')

   elif Molecule[0][0] == 'RC:':
      print ('--rcpath ', filepath, ' \\')
      print ('--rcfragment ',defragment(filepath),' \\')

   elif Molecule[0][0] == 'P:':
      print ('--ppath ', filepath, ' \\')
      print ('--pfragment ',defragment(filepath),' \\')

   os.chdir(currentPath)
