import numpy as np

# Calcualte bond energy of He dimers for a series of bond 
# distances using ADF

# type of atoms
atom1 = 'He'
atom2 = 'He'

# interatomic distance values
dmin = 2.2
dmax = 4.2
step = 0.2

# create a list with interatomic distances
distances = np.arange(dmin,dmax,step)

# calculation parameters (single point, TZP/PBE+GrimmeD3)
sett = Settings()
sett.input.ams.task = 'SinglePoint'
sett.input.adf.basis.type = 'TZP'
sett.input.adf.xc.gga = 'PBE'
sett.input.adf.xc.dispersion = 'Grimme3'

energies = []
for d in distances:
    mol = Molecule()
    mol.add_atom(Atom(symbol=atom1, coords=(0.0, 0.0, 0.0)))
    mol.add_atom(Atom(symbol=atom2, coords=(  d, 0.0, 0.0)))
    job = AMSJob(molecule=mol, settings=sett, name=f'dist_{d:.2f}')
    job.run()
    energies.append(job.results.get_energy(unit='kcal/mol'))

# print
print('== Results ==')
print('d[A]    E[kcal/mol]')
for d,e in zip(distances, energies):
    print(f'{d:.2f}    {e:.3f}')
