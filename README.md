The DESIGN Package for GAP
==========================

This README describes the DESIGN package for GAP. The DESIGN package is
for constructing, classifying, partitioning and studying block designs.

All DESIGN functions are written entirely in the GAP language.
However, DESIGN requires the GRAPE package (version
at least 4.4) to be installed, and makes use of certain GRAPE
functions, some of which make use of B.D.~McKay's nauty package.
These GRAPE functions can only be used on a fully installed version
of GRAPE.

The DESIGN package is Copyright (C) Leonard H. Soicher
2003--2011.  DESIGN is part of a wider project, which received EPSRC
funding under grant GR/R29659/01, to provide a web-based resource for
design theory; see <http://designtheory.org>.

DESIGN is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version. For details, see 
<http://www.gnu.org/licenses/gpl.html>.

If you use DESIGN to solve a problem then please send a short email
about it to \Mailto{L.H.Soicher@qmul.ac.uk}, and reference the DESIGN 
package as follows:

L.H. Soicher, The DESIGN package for GAP, Version 1.6, 2011,
http://designtheory.org/software/gap_design/.

Any comments or bug reports should also go to
<L.H.Soicher@qmul.ac.uk>.

Installing the DESIGN Package
-----------------------------

The DESIGN package is included in the standard GAP
distribution. You only need to download and install DESIGN if you need
to install the package locally or are installing an upgrade of DESIGN
to an existing installation of GAP (see the main GAP reference
section "Reference: Installing a GAP Package").  If you do need to download
DESIGN, you can find archive files for the package in various formats
at <https://www.gap-system.org/Packages/design.html>, and then your
archive file of choice should be downloaded and unpacked in the `pkg'
subdirectory of an appropriate GAP root directory (see the main GAP
reference section "Reference: GAP Root Directories").

The DESIGN package is written entirely in GAP code, and requires
no further installation.  However, the DESIGN package has complete
functionality only in an environment in which the GRAPE package is
fully installed.  Usually this will be on a Unix system.

-------------------------------------------------------------------------
Leonard Soicher
Queen Mary, University of London
November 2011
