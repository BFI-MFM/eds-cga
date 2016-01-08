Dynare
======

For information about how to use Dynare, you should have a look at the
documentation located in the 'doc' subdirectory of your Dynare installation (you
should have a shortcut in your Start Menu to access it directly).

Beginners should start with the Dynare user guide (under PDF format in
'UserGuide.pdf'). There is also a complete reference manual documenting all
Dynare functions (under HTML format in 'dynare.html\index.html', under PDF
format in 'dynare.pdf').

You can also get more information on the web, on Dynare homepage:

  http://www.dynare.org

Or on Dynare Wiki:

  http://www.dynare.org/DynareWiki

NOTE: Dynare comes with an automated uninstaller, which you can run from the
"Add/Remove Programs" menu of the Control Panel.


Using Dynare with MATLAB (R)
----------------------------

Dynare requires MATLAB (R) version 7.0 or above. With older versions of MATLAB (R),
it may fail or give unexpected results.

To use Dynare, you just have to add the 'matlab' subdirectory of your Dynare
installation to MATLAB (R) path. You have two options for doing that:

* Use the addpath command, by typing the following (assuming that you have
  installed Dynare at the standard location, and replacing '4.x.y' by correct
  version number):

    addpath c:\dynare\4.x.y\matlab

  MATLAB (R) will not remember this setting next time you run it, and you will
  have to do it again.

* Select the "Set Path" entry in the "File" menu, then click on "Add
  Folder...", and select the 'matlab' subdirectory of your Dynare
  installation. Note that you SHOULD NOT use "Add with Subfolders...". Apply
  the settings by clicking on "Save". Note that MATLAB (R) will remember this
  setting next time you run it.

You can test your installation by typing 'dynare' at the MATLAB (R)
prompt. This should give you an error message complaining that you did not
specify a MOD file.


Using Dynare with Octave
------------------------

Dynare is now available for GNU Octave, a free clone of MATLAB (R) (see
<http://www.octave.org>).

The recommended Octave distribution is the Octave/MinGW 3.2.4
precompiled binaries from Octave Forge, available at:

  http://sourceforge.net/projects/octave/files/Octave_Windows%20-%20MinGW/

Every time you run Octave, you should type the two following commands (assuming
that you have installed Dynare at the standard location, and replacing '4.x.y'
by correct version number):

  addpath c:\dynare\4.x.y\matlab

NOTE: if you don't want to type this command every time you run
Octave, you can put it in a file called '.octaverc' in your home
directory (generally 'c:\Documents and Settings\USERNAME\'). This file
is run by Octave at every startup.

You can test your installation by typing 'dynare' at the Octave prompt. This
should give you an error message complaining that you did not specify a MOD
file.

For more information about Dynare for Octave, go to:

  http://www.dynare.org/DynareWiki/DynareOctave


Dynamic Loadable Libraries
--------------------------

For better performance, some parts of Dynare are written in the C++ language,
which is faster than standard M-files. These parts are compiled and distributed
as dynamic loadable libraries (DLL), located in the 'mex' subdirectory of your
Dynare installation.

If the DLL are correctly detected by MATLAB (R) or Octave, the following should
be displayed when you launch Dynare:

  Configuring Dynare ...
  [mex] Generalized QZ.
  [mex] Sylvester equation solution.
  [mex] Kronecker products.
  [mex] Sparse kronecker products.
  [mex] Bytecode evaluation.
  [mex] k-order perturbation solver.
  [mex] k-order solution simulation.

On the contrary, if DLL are not detected, Dynare will fallback on
slower alternatives written in M-files (only for some of the DLLs),
and display the following:

  Configuring Dynare ...
  [m]   Generalized QZ.
  [m]   Sylvester equation solution.
  [m]   Kronecker products.
  [m]   Sparse kronecker products.
  [no]  Bytecode evaluation.
  [no]  k-order perturbation solver.
  [no]  k-order solution simulation.

In this last case, Dynare will run correctly on the basic features,
but with suboptimal speed, and some features will be missing. There
could be several reasons for MATLAB (R) or Octave failing to detect
the DLL:

* Your path settings may be wrong. Make sure that the 'matlab' subdirectory of
  your Dynare installation is the only Dynare directory present in the path
  variable.

* Your MATLAB (R) or Octave version may be incompatible with the provided
  binaries.

* You may have a custom M-file in your search path with the same name than a
  DLL, therefore overriding it.


Credits
-------

MATLAB (R) is a registered trademark of The Mathworks, Inc.
