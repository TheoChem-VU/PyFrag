Main Specifications
====================

You can print more information in the final activation strain analysis by adding the following specifications between the PyFrag and PyFrag END
in the job input.

These following statements allow you to define the fragments in the analysis. For each of the fragments you supply a list of the numbers for the atoms as they exist on the LT or IRC file. The program will check if they match with the order as in the IRC or LT calculation file. It will print a statement in the error log if it thinks the order is wrong. ::

  fragment atomnrs

  for example:

  fragment 2
  fragment 1 3 4 5 6

The following possibilities are optional. Users can choose what information to print after a PyFrag calculation. For instance, one can specify if the strain energy is to be printed for the fragments. Or one can specify the equilibrium energies (in kcal/mol) for the fragments. Beware that the order of this specification should correspond to the order of the fragment definition. This value will then simply be subtracted from the energy of the fragment in question. The program will print the individual strain values for each fragment, plus the total strain and total energy. ::

  strain fragenergy

  for example:

  strain -301.01
  strain -19.02

To specify the bond length to be printed for each geometry step, just indicate the atom numbers as they appear in the input order of the total molecule. Specifying bond_diff as well subtracts this value from the actual bond length. ::

  bondlength atomnr1 atomnr2 bond-diff

To specify the angle between atoms 1, 2 and 3 to be printed for each geometry step, just indicate the atom numbers as they appear in the input order of the total molecule. Specifying angle_diff as well subtracts this value from the actual angle. ::

  angle atomnr1 atomnr2 atomnr3 angle-diff

The following will print the Hirshfeld charges for the fragment. Note that the order of the fragments as used internally by ADF may differ from what you would expect. Note also that Hirshfeld charges as computed in a fragment analysis differ from those obtained in a ‘normal’ calculation from basic (spherical average-of configurations) ADF atoms like simple single point calculations. ::

  hirshfeld frag1

The following prints the VDD charges on the atoms with numbers as given by atomnrs. Note also that VDD charges as computed in a fragment analysis differ from those obtained in a ‘normal’ calculation from basic (spherical average-of configurations) atoms like single point calculation. ::

  VDD atomnrs

  for example:

  VDD 1 2

Through in the input below, the value of the orbital interaction energy will be printed for certain available irrep. The irrep symbol relates to the symmetry of the whole molecule. ::

  irrepOI oi irrep

  for example:

  irrepOI AA

The statement below will print the overlap between orbital numbers orb1 on frag1 and orb2 on frag2 in irrep as they appear in the fragment analysis calculation. It can also print the overlap between the HOMO's of fragment one (frag1) and the LUMOs of fragment two (orbitals as found on the fragment calculations). If the orbitals found in this way differ in symmetry, an overlap value of zero is returned. It should be noted that irrep and orbital number refers to each fragment, rather than the whole molecule. Especially, when it comes to frozen core situation, the count of orbital number should not include core orbitals. This rule also applies to the situation of printing orbital energy and population for orbitals of certain fragment. ::

  overlap frag1 HOMO/LUMO frag2 HOMO/LUMO
  overlap irrep1 frag1 orb1 irrp2 frag2 orb2

  for example:

  overlap frag1 HOMO frag2 LUMO
  overlap S frag1 5 AA frag2 4

Through the input below, the orbital energy for a fragment orbital will be printed in available irrep. ::

  orbitalenergy frag HOMO/LUMO
  orbitalenergy irrep frag orb

  for example:

  orbitalenergy frag1 HOMO
  orbitalenergy AA frag2 5

Through the input below, the Mulliken population for a fragment orbital will be printed in available irrep. ::

  opulation frag HOMO/LUMO
  population irrep frag orb

  for example:

  population frag1 HOMO
  population AA frag2 5
