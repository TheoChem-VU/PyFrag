# Changelog
All notable changes to this project will be documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),


## PyFrag 2025 - unreleased

### Added
 * PyFrag jobs now automically restart when a folder with the same name already exists. The old folder will get the “.res” extension and will be deleted when the new calculation finished succesfully. This is to prevent quick reprinting of the results table upon changing PyFrag specifications such as "overlap", "bondlength", etc. without having to redo all computations.
 * A logger which level is determined by the “log_level”, resulting in a “[calculation_name].log” file. Try to specify “log_level debug” in you input within the PyFrag section to see all the details of the calculation. The log file is saved in the same folder as the input file.
 * A top line header in the output table which contains information about the units of the keys. For example, the key "bondlength" now contains the unit "Angstrom" in the header.
 * The keys "frag1_indices", "frag2_indices", "frag1_energy", "frag2_energy" as these are more informative than the previous keys "fragment" and "strain" which were not clear about the meaning of the numbers. The old keys are still available for backward compatibility.

### Changed
 * Inputfiles are now compatible with AMS/ADF input files before and after 2019. We recommended using the >AMS2019 inputformat since this receives frequent updates and keeps being supported by plams.
 * The All ADF keys are now recognized and correctly parsed including UnrestrictedFragments, FragOccupations, and IrrepOccupations.
 * The calculation folders (e.g., [name].001]) are cleaned up after the calculation has finished successfully to reduce disk space.
 * The header of the output table contains more information about the relevant key. For example, the key "bondlength" now includes the atom numbers and symbols of the atoms involved in the bondlength calculation. The same applies to the keys "overlap" -> "overlap_[orbital1_index]-[orbital1_index]_[orbital2_index]-[orbital2_index]", "population", "orbitalenergy", and "angle". [#15](https://github.com/TheoChem-VU/PyFrag/issues/15)
 * The unrestricted and unrestricted modules have merged as there was considerable overlap which made the code difficult to maintain. Both are now in the adf_new folder in the host/standalone folder. The adf_new folder is now the default folder for all ADF calculations.
 * The input blocks such as "AMS", "PyFrag", "JobSub", "Complex EXTRA" etc. are now case-insensitive. This means that the input keys can be written in any case, e.g., "pyfrag", "PyFrag", "PYFRAG" etc. will all be recognized as the same key.

### Removed
 * Removed the IrrepOI input key as it prints now automatically the orbital interactions decomposed into symmetry irreps (if symmetry is used) [#10](https://github.com/TheoChem-VU/PyFrag/issues/10)
 * The input keys "ircpath", "lt", "irctwo" are now deprecated and replaced by the general "coordfile" key. The program decides now automatically if the file is a coordinate file or an irc path based on the file extension. The program raises a warning if the deprecated keys are used.
 * The "newopen" executable (`pyfrag -x newopen`) is now deprecated and replaced by the "adf" executable. The program raises a warning if the deprecated key is used.

### Fixed
 * Ordering of the header in the output table (pyfrag_[jobname].txt) with correct spacing. Keys such as "bondlength10" will not be printed before "bondlength2" anymore, but after "bondlength9". So, applying a natural order scheme instead of the default ordering of the sorted python function [#14](https://github.com/TheoChem-VU/PyFrag/issues/14).
 * Reading in charges in the Systems block goes correctly now [#11](https://github.com/TheoChem-VU/PyFrag/issues/11)
 * Getting an error when multiple bondlengths were specified in the input file.
 * Reading in .amv files which would have (by accident) an header line that splits in exactly four parts, which would correspond to a "atom line" such as "C 0.0000 0.0000 0.0000" [#6](https://github.com/TheoChem-VU/PyFrag/issues/6)

## PyFrag2019 - 19/10/2018

### Added
 * Quantum Dots builder functionality

### Changed

 * Use [noodles==0.3.0](https://github.com/NLeSC/noodles/releases)
 * Replace [nose](https://nose.readthedocs.io/en/latest/) with [pytest](https://docs.pytest.org/en/latest/)
 * Imported only the core functionality of the library
 * Used [intra-package-references](https://docs.python.org/3/tutorial/modules.html#intra-package-references)
 * Used `__all__` to limit exposed fuctionality

### Removed

 * Dead code related to all noodles API
 * All the `import *`
 * Dead code from components including PES

### Fixed

 * Job manager issue when removing a SCM job


## PyFrag2019 - 19/10/2018

### Changed
divide all code into three level: pyfrag module; easy job submit and restart and summary(without the requirement of open local server and applescript); full functionalities.

