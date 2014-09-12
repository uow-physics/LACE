README of the LACE package
==========================
Contact: John Back (J.J.Back@warwick.ac.uk)

Introduction
------------

This package (Lpc Algorithm in a C++ Environment) is a C++ 
implementation of the local principal curve (Lpc) algorithm 
described in Eur. Phys. J. C (2014) 74:2832, also available 
on the arXiv: http://arxiv.org/abs/1312.6059. 
It is based on the R-code package "LPCM"
http://cran.r-project.org/web/packages/LPCM/index.html,
written by Jochen Einbeck and Ludger Evers.

This software takes a collection of data points and tries 
to find their local principal curves. Each point is represented 
as an "LpcHit" with N-dimensional co-ordinates (x,y,z,...) and a 
weight (e.g. energy or charge). Each collection of hits 
(data point cloud) is stored in a so-called "LpcEvent" pointer. 
The curves, feature points and any clustering information can 
then be retrieved from each LpcEvent.

Further online documentation is available at
http://universityofwarwick.github.io/LACE


Building the code
-----------------

The package only depends on Eigen, the C++ template library 
for linear algebra: http://eigen.tuxfamily.org. Since Eigen
only contains header files, there is no need to "install" it;
only the Eigen header-file directory needs to be present.

This software can also use ROOT, a general purpose C++ data 
analysis framework (http://root.cern.ch), for input and 
output files, and for finding specific features in the 
principal curves. However, using ROOT is not required nor 
mandatory.

Use the configure script to set various compiler options
("configure --help" lists the available options), then
build using make:

```sh
$ ./configure
$ make
```

If successful, this should create a shared library 
`lib/libLACE.so` as well as the binary program `bin/LACEMain`,
which is based on the example [LACEMain.cc](LACEMain.cc) 
file in the base LACE directory.

The file [pathlib.sh](pathlib.sh) gives an example of setting the
`LD_LIBRARY_PATH` environment variable to include the LACE
shared library.


Running the code
----------------

The example program `bin/LACEMain` can be used to run the 
LACE code. It takes one extra argument on the command line, 
which specifies the name of the Lpc parameter file

```sh
$ cd example
$ ../bin/LACEMain mupPars.txt
```

An example of ten 770 MeV neutrino to muon-proton 
events is provided in the [mup770MeV.txt](example/mup770MeV.txt) 
file, which is text based and contains the following lines:

```
NumberOfEvents NumberOfDimensions
EventNumber1 NumberOfHits1
HitX1 HitX2 HitX3 HitWeight
Other hits...
EventNumber2 NumberOfHits2
HitX1 HitX2 HitX3 HitWeight
Other hits...
etc...
```

A description of the lpc parameters are given below, with suggested 
default values given in parenthesis:

```
infile          Filename of the input file
informat        Format for the input file (text or root)
outfile         Filename of the output file
outformat       Format for the output file (text or root)
 
kernelwidth     Value of the lpc scaled kernel width (0.05)
stepsize        Scaled step size for the lpc (0.05)
npoints         Required number of lpc points (250)
penalisation    Angle penalisation factor (2.0)
eigenratio      Minimum 2nd/1st eigenvalue ratio for possible branching points (0.4)
boundary        Boundary condition for delaying convergence in the tails (0.005)
convergence     Convergence level for the Lpc (1e-6)
branchlevel     Number of possible branching generation levels (0 = none, default)
gapsize         Scaled size of the gap (gapsize*stepsize) to start finding a new branch (1.5)

peakfinder      Select the peak finder method: 0 = none, 1 = simple (default if not using ROOT), 
                                               2 = TSpectrum from ROOT (default if using ROOT)

cosanglecut     Minimum value for finding features in the 1-|cosAngle| distribution (0.01)
minpeakfrac     Minimum fractional height for the next feature peak w.r.t the previous peak (0.01)
peakdiffsq      Squared-distance limit for checking if lpc feature peaks overlap (3.0)

clustering      Set the clustering method: 0 = none, 1 = method in paper (default), 
		                           2 = all hits put into 1 cluster

dobranchvtx     Enable (1) or disable (0 = default) vertexing and clustering for branches
minvtxrescut    Minimum selection cut value for hit-to-lpc residuals for vertexing (20.0)

convexhull      Minimum value of the convex hull ratio d_transverse/d_longitudinal for showers (0.12)
showerres       Minimum threshold for hit-to-lpc residuals for showers (20.0)
showerresratio  Minimum hit-to-lpc residual ratio for showers (0.3)
showerresfrac   Minimum fraction of residuals that need to be above the shower threshold (0.9)

firstevent      Integer specifying the first event (usually this is set to 0)
lastevent       Integer specifying the last event (-1 means process all available events)
```

If "firstevent" and "lastevent" are not set, then all events found in the input 
file are processed.


This software supports finding branches within lpc curves ("branchlevel > 0"),
but this is still in the experimental stage. Feedback regarding this
feature is desired and indeed welcome.

Note that all integer values used and returned by functions and variables
follow the C++ convention and start at zero, _except for branches_.
The integer index number of the main curve is nominally set to zero.
This implies that the indices of branches start at one (branch = 0 means
the main curve that has index 0). Also, the branch generation number 
starts at 1 (generation of zero means no branches, or just the main curve).


Code documentation
------------------

Automatic code documentation can be generated using the
Doxygen system (http://www.doxygen.org). Simply run doxygen

```
$ doxygen
```

in the main LACE directory to create a doxygen sub-directory 
containing a html-based interface to the description of classes 
and all of their functions and variables (doxygen/html/index.html).


License
-------

This software is distributed under the Boost Software License, Version 1,
(Aug 17 2003). See [LICENSE_1_0.txt](LICENSE_1_0.txt) 
(optionally the original at http://www.boost.org/LICENSE_1_0.txt)
for details.


Authors and contributors
------------------------

The main author of the C++ version is John Back.

Other contributors to this work are:

* Gary Barker
* Steve Boyd
* Daniel Brunt
* Jochen Einbeck
* Ludger Evers
* Martin Haigh
* Harmanjeet Khera
* Ben Morgan
* Ben Oakley
* Yorck Ramachers
* Dan Roythorne
* Jamie Wynn

