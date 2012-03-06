Changes
=======

Release 2.0
^^^^^^^^^^^
Bopix now working correctly.
New CAMspec version.
New smica version.
New low l likelihood lowlike.

Better support of lapack mkl. It will now work on more infrastructure. Add of Apple lapack support. 

Simplified (again) configuration procedure. Better (more streamlined) config lines.


Release 1.6
^^^^^^^^^^^
Addition of the egfs foreground model from `arXiv:1102.5195 <http://arxiv.org/abs/1102.5195>`_. This model can be added to inverse gamma, gaussian and smica likelihoods with parametrizable list of default and free parameters. Added (undocumented) tools to include the model in clik files.

Numerous bug fixes. 

Simplified (again) install procedure. More robust test for the availability of external dependencies.
Recommandation is now to let clik install all the external dependencies.
Added ``--install_all_deps`` to implement this recommandation.
Removed the ``--local`` option, install is now local by default. 

Release 1.5
^^^^^^^^^^^

Simplified (hum, yes, it's much simpler) install procedure. Do not require pmclib any more. Better test for the requirements.

Add (fast) low-l pixel based likelihood T+P, BOPIX likelihood. Can call the WMAP7 likelihood code along with the WMAP7 data.

Add tools to manipulate the likelihood files.


Release 1.0
^^^^^^^^^^^

Initial release.

Smica, inverse gamma, inverse gamma copula and gaussian likelihood approximations.

Tool to generate smica likelihoods.