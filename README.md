## Licensing

FRIDA (fast reliable inelastic data analysis) is a program (<https://sourceforge.net/projects/frida/>, <https://jugit.fz-juelich.de/mlz/frida>) for generic spectral analysis, with many specialized routines for inelastic neutron scattering. The FORTRAN version Frida-1 is an updated version of Joachim Wuttke's IDA, with contributions from the community. The maintainer is Florian Kargl <f_kargl@users.sourceforge.net>.
Frida-1 is released under the GNU General Public License (GPL).
The modifications made by Artem Panchenko <artem.panchenko@fau.de> are released under the GNU General Public License version 3 (GPLv3).
These modifications add support for the 2025 version of PSI FOCUS instrument data files and to enable compilation with modern GNU Fortran, including GNU Fortran 11.4.0. A detailed list of the modifications is provided in `ChangeLog_icsp.txt`.

(C) Joachim Wuttke 1990-2001
(C) Florian Kargl 2006

Third-party numerical software included in this repository retains its original licensing or public-domain status. Original source comments, copyright notices, and attribution information have been preserved.

### MINPACK

This project contains selected routines from MINPACK distributed by Netlib:

https://www.netlib.org/minpack/

These routines are distributed under the Minpack License (`SPDX-License-Identifier: Minpack`).

The original source comments, copyright notices, and license terms have been retained. See `minpack/disclaimer`.

### SLATEC

This project contains selected routines from the SLATEC Common Mathematical Library, Version 4.1 (July 1993), distributed by Netlib:

https://www.netlib.org/slatec/

The SLATEC documentation identifies the library as public-domain software. Original source comments, authorship information, references, and attribution have been retained.

The SLATEC routines used by this project originate from several components of the Netlib SLATEC distribution:

- SLATEC core routines:
  https://www.netlib.org/slatec/src/

- FNLIB special-function routines:
  https://www.netlib.org/slatec/fnlib/

- FISHFFT / FFT routines:
  https://www.netlib.org/slatec/fishfft/

The FISHFFT routines include Fast Fourier Transform routines by P. N. Swarztrauber (NCAR), as identified in the Netlib distribution.

See `slatec/readme`.

### LINPACK

This project contains selected routines originating from the LINPACK library distributed by Netlib:

https://www.netlib.org/linpack/

LINPACK is historically distributed as public-domain software. The original source comments, authorship information, and attribution have been retained.

The principal authors of LINPACK are Jack J. Dongarra, Cleve B. Moler, J. R. Bunch, and G. W. Stewart.

See `linpack/readme`.

### NeXus API

This project includes source and interface files from the NeXus API.

Copyright © NeXus International Advisory Committee.

The NeXus API files are distributed under the GNU Lesser General
Public License, version 2 or, any later version
(SPDX-License-Identifier: LGPL-2.0-or-later).

Original copyright and license notices have been retained.

https://www.nexusformat.org/
