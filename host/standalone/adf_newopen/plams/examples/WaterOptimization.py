# Perform a geometry optimization of a water molecule and compute
# the vibrational normal modes using GFN1-xTB.

# You could also load the geometry from an xyz file: 
# molecule = Molecule('path/my_molecule.xyz')
molecule = Molecule()
molecule.add_atom(Atom(symbol='O', coords=(0,0,0)))
molecule.add_atom(Atom(symbol='H', coords=(1,0,0)))
molecule.add_atom(Atom(symbol='H', coords=(0,1,0)))

settings = Settings()
settings.input.ams.Task = 'GeometryOptimization'
settings.input.ams.Properties.NormalModes = 'Yes'
settings.input.dftb.Model = 'GFN1-xTB'

job = AMSJob(molecule=molecule, settings=settings, name='water_optimization')
result = job.run()

energy = result.get_energy(unit='kcal/mol')
frequencies = result.get_frequencies(unit='cm^-1')
optimized_molecule = result.get_main_molecule()

# Unlike python lists, where the index of the first element is 0, 
# the index of the first atom in the molecule object is 1
bond_angle = optimized_molecule[1].angle(optimized_molecule[2], optimized_molecule[3])

print('== Results ==')
print('Optimized geometry:')
print(optimized_molecule)
print('Energy      : {:.3f} kcal/mol'.format(energy))
print('Bond angle  : {:.1f} degrees'.format(Units.convert(bond_angle, 'rad', 'degree')))
print('Frequencies : {} cm^-1'.format(frequencies))