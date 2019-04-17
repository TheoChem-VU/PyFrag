Further Information
===================

History of PyFrag
-----------------

**PyFrag 2008**

The original version of PyFrag was developed by Willem-Jan van Zeist, Lando P. Wolters, F. Matthias Bickelhaupt, and Célia Fonseca Guerra at the Theoretical Chemistry department at Vrije Universiteit Amsterdam. PyFrag is mainly used to enables a user-friendly analysis of reaction paths in terms of the Extended Activation Strain model of chemical reactivity(ASM). The explanation and application can be found in a list of liteartures as below. For old users more useful information can be found on here_ or more `concise version`_. Noted this version is not maintained anymore.

**PyFrag 2016**

The new PyFrag program was rewritten by Xiaobo Sun and Thomas Soini using the PLAMS library in the ADF package and has been included in the script collection in ADF_ 2017 and later version. Compared to the old PyFrag, the new version is more compact and easy to be maintained, expanded and upgraded. Also, due to its high compatibility with other python library tools developed by SCM, such as, PLAMS and QMworks, it can be used as a module in line with these computational chemistry job management tools to streamline a large flow of job. For this version, description can be found in `this site`_.

**PyFrag 2019**

Current version. The PyFrag program is specially designed to facilitates the study of reaction mechanism in a more efficient and user-friendly way. It automates the process of finding transition states, potential energy surface by using one simple input file. It follows by an activation strain analysis on the energy profile to characterize the feature of the reaction mechanism and gain insights into the overall reaction energies. Moreover, users can have an real-time monitoring of the running process via a webpage which vividly displays the updated data in the form of videos and figures and, if necessary, user can rerun the job immediately from where it stops. In this way, the three respects of computational chemistry–job management, data management and analysis management can all be contained in a framework and thus allow chemists to focus on the interpretation and creation work rather than waste time and energy on the finding and processing of massive data.

Activation Strain Model
-----------------------

For information of Activation Strain model of chemical reactivity(ASM), reader can refre to literatures. An easy exercise_ for activation strain analysis of reaction mechanism using ADF and PyFrag is also included.

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
