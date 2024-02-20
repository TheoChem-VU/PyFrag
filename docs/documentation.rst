Further Information
===================

PyFrag Log & Bugfixes
---------------------

**PyFrag 2019.0.1**

Bugfixes: 
 - Inputfiles are now compatible with AMS/ADF input files before and after 2019. We recommended using the >AMS2019 inputformat since this receives frequent updates and keeps being supported by plams. The difference in PyFrag inputfiles is that the ADF / ADF END has been changed to AMS / AMS END, and that the content within the AMS / AMS END should be the same as specified in the ADF GUI for >AMS2019 with an "engine", "system", and "task" block
 - All ADF keys are now recognized such as UNRESTRICTEDFRAGMENTS.
 - The calculation folders (e.g., [name].001]) are cleaned up after the calculation has finished successfully to reduce disk space. 
 - PyFrag jobs now automically restart when a folder with the same name already exists. The old folder will get the ".res" extension and will be deleted when the new calculation finished succesfully. This is to prevent multiple copies and reducing disk space.
 - A logger has been added through the key "log_level" which makes a "[calculation_name].log" file. Try to specify "log_level debug" in you input within the "PyFrag" section
 - Printing the relativistically-scaled orbital energies instead of unscaled energies if scalar ZORA is specified
 - Automatically optimizing fragments when strain keys are not specified

History of PyFrag
-----------------

**PyFrag 2008**

The original version of PyFrag 2008 was developed by Willem-Jan van Zeist, Lando P. Wolters, F. Matthias Bickelhaupt, and Célia Fonseca Guerra at the Theoretical Chemistry Department at the Vrije Universiteit Amsterdam. PyFrag is mainly used to enable a user-friendly analysis of reaction paths in terms of the Extended Activation Strain model of chemical reactivity(ASM). The explanation and application can be found in the references below. For users still using this version, additional useful information can be found on  here_ or here for a more `concise version`_. This version is no longer maintained.

**PyFrag 2016**

The PyFrag 2016 program was rewritten by Xiaobo Sun and Thomas Soini using the PLAMS library in the ADF package and has been included in the script collection in ADF_ 2017 and later version. Compared to the old PyFrag, the new version is more compact and easy to be maintained, expanded and upgraded. Also, due to its high compatibility with other python library tools developed by SCM, such as, PLAMS and QMworks, it can be used as a module in line with these computational chemistry job management tools to streamline a large flow of job. For this version, description can be found using `this site`_.

**PyFrag 2019**

The PyFrag 2019 program was specially designed to facilitate the analysis of reaction mechanism in a more efficient and user-friendly way. PyFrag 2019 has automated and reduced the time-consuming and laborious task of setting up, running, analyzing, and visualizing computational data from reaction mechanism studies to a single job. PyFrag 2019 resolves three main challenges associated with the automatized computational exploration of reaction mechanisms: 1) the management of multiple parallel calculations to automatically find a reaction path; 2) the monitoring of the entire computational process along with the extraction and plotting of relevant information from large amounts of data; and 3) the analysis and presentation of these data in a clear and informative way. The activation strain and canonical energy decomposition results that are generated, relate the characteristics of the reaction profile in terms of intrinsic properties (strain, interaction, orbital overlaps, orbital energies, populations) of the reactant species.

Activation Strain Model
-----------------

For more information on the Activation Strain Model (ASM) of chemical reactivity, the user is directed to the references provided below. An easy exercise_ for activation strain analysis of reaction mechanism using ADF is also included.

**Literature** ::

  1 W.-J. van Zeist, C. Fonseca Guerra, F. M. Bickelhaupt, J. Comput. Chem. 2008, 29, 312-315.
  2 I. Fernandez, F. M. Bickelhaupt, Chem. Soc. Rev. 2014, 43, 4953-4967.
  3 L. P. Wolters, F. M. Bickelhaupt, WIRES Comput. Mol. Sci. 2015, 5, 324-343.
  4 F. M. Bickelhaupt, K. N. Houk Angew. Chem. 2017, 129, 10204-10221; Angew. Chem. Int. Ed. 2017, 56, 10070-10086.


.. _here : http://www.few.vu.nl/~xsn800/Home.html
.. _concise version: https://sunxb05.github.io/pyfragold/
.. _ADF: https://www.scm.com/doc/ADF/Input/PyFrag.html
.. _this site: http://www.few.vu.nl/~bickel/page-2/pyfrag.html
.. _exercise: https://github.com/sunxb05/PyFrag/blob/master/docs/exerciseforPyFrag.docx
