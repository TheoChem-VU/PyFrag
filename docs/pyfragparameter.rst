Main Specifications
====================

The user can print additional information in the final activation strain analysis by adding the following specifications between the PyFrag and PyFrag END in the job input.

These following statements allow the user to define the fragments in the analysis. For each of the fragments you supply a list of the numbers for the atoms as they exist in the supplied XYZ coordinate file of the reaction path from an IRC or LT calculation. The program will check if they match with the order in the supplied XYZ coordinate file of the reaction path from an IRC or LT calculation. If the atom ordering is incorrect, a statement will be printed in the error log. ::

  fragment atomnrs

  for example:

  fragment 2
  fragment 1 3 4 5 6

The following possibilities are optional. The user can choose what information to print during the Activation Strain Analysis (ASA). For instance, the user can specify to print the strain energy for the fragments. Or one can specify the equilibrium energies (in kcal/mol) for the fragments. Beware that the order of this specification should correspond to the order of the fragment definition. This value will then simply be subtracted from the energy of the fragment in question. The program will print the individual strain values for each fragment, plus the total strain and total energy. ::

  strain fragenergy

  for example:

  strain -301.01
  strain -19.02

To specify the bond length to be printed for each geometry step, the user needs to indicate the atom numbers as they appear in the input order of the total molecule. Specifying bond_diff as well subtracts this value from the actual bond length. ::

  bondlength atomnr1 atomnr2 bond-diff

To specify the angle between atoms 1, 2, and 3 to be printed for each geometry step, just indicate the atom numbers as they appear in the input order of the total molecule. Specifying angle_diff as well subtracts this value from the actual angle. ::

  angle atomnr1 atomnr2 atomnr3 angle-diff

The following statement will print the Hirshfeld charges for the fragment. Note that the order of the fragments as used internally by ADF may differ from what you would expect. Note also that Hirshfeld charges as computed in a fragment analysis differ from those obtained in a ‘normal’ calculation from basic (spherical average-of configurations) ADF atoms like simple single point calculations. ::

  hirshfeld frag1

The following statement will print the VDD charges on the atoms with numbers as given by atomnrs. Note also that VDD charges as computed in a fragment analysis differ from those obtained in a ‘normal’ calculation from basic (spherical average-of configurations) atoms like single point calculation. ::

  VDD atomnrs

  for example:

  VDD 1 2

The following statement will print the orbital interaction energy per available irrep. The irrep symbol relates to the symmetry of the whole molecule. ::

  irrepOI oi irrep

  for example:

  irrepOI AA

The following statement will print the overlap between orbital numbers orb1 on frag1 and orb2 on frag2 in irrep as they appear in the fragment analysis calculation. It can also print the overlap between the HOMOs of fragment 1 (frag1) and the LUMOs of fragment 2 (orbitals as found on the fragment calculations). If the orbitals found in this way differ in symmetry, an orbital overlap value of zero is returned. It should be noted that irrep and orbital number refers to each fragment, rather than the whole molecule. Especially, when it comes to frozen core situation, the count of orbital number should not include core orbitals. This rule also applies to the situation of printing orbital energy and population for orbitals of certain fragment. ::

  overlap frag1 HOMO/LUMO frag2 HOMO/LUMO
  overlap irrep1 frag1 orb1 irrp2 frag2 orb2

  for example:

  overlap frag1 HOMO frag2 LUMO
  overlap frag1 HOMO-1 frag2 LUMO+3
  overlap S frag1 5 AA frag2 4

The following statement will print the orbital energy for a fragment orbital per available irrep. The irrep symbol relates to the symmetry of the fragment. ::

  orbitalenergy frag HOMO/LUMO
  orbitalenergy irrep frag orb

  for example:

  orbitalenergy frag1 HOMO
  orbitalenergy frag1 HOMO-2
  orbitalenergy AA frag2 5

The following statement will print the gross Mulliken population for a fragment orbital. ::

  population frag HOMO/LUMO
  population irrep frag orb

  for example:

  population frag1 HOMO
  population frag2 HOMO-1
  population AA frag2 5
