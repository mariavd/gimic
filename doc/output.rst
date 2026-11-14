

Description of generated files
==============================

Note "file.txt" means text file format and "file.vti" stands for VTK file format, which is compatible with ParaView. The text files start with a given grid point (x,y,z) followed by either information about the modulus of the current density f(x,y,z) or the components of the current density vector Jx,Jy,Jz.

The ".vti" files are binary, which keeps them small and quick to load in
ParaView. They cannot be read as text; use ``vti2csv`` (see below) to get the
numbers as a table, or set ``vti_format=ascii`` in the Advanced section of the
input to have GIMIC write them as text.

* ``MOL``: contains molecular coordinates and basis functions, needed as input file for a gimic calculation
* ``XDENS``: contains information about the AO density and first order perturbed AO density, needed as input for a gimic calculation
* ``gimic.inp``: gimic input file, here either the cubic grid or a grid across a bond is specified, needed as input information for a gimic calculation
* ``gimic.out``: gimic output file, contains information about the integrated current strength susceptibility
* ``grid.xyz``: contains molecular coordinates in Ångstrøm and cube points (8) or integration plane corner points (4)
* ``current_profile.dat``: with ``slice_thickness`` in the grid section: the current through each slice of the integration plane, one line per slice from the start of the integration: distance (bohr), net, positive and negative contribution (nA/T)
* ``vertical_profile.dat``: the same for the slices along the height of the plane, when the height is sliced
* ``profile_points.dat``: maxima and minima of the net, positive and negative current per bohr along the profile, the zero crossings of the net current and the current integrated between them, and the blobs of the positive and the negative contribution with their boundaries, peaks and integrated currents
* ``profile.dat``: the same with the geometry of the plane, the slice centres in space, the slices along the height if any, and the current per unit length at every quadrature node
* ``squares_profile.dat``: when both axes are sliced, the net current through each rectangle as a matrix and as a list with coordinates
* ``mol.xyz``: contains molecular coordinates in Ångstrøm

closed-shell calculation:

* ``acid.vti``: contains information about ACID function 
* ``jvec.vti``: contains information about the current density vector function 
* ``jmod.txt``: contains a grid point and the modulus of the current density vector function, x,y,z, f_mod(x,y,z)
* ``jmod.vti``: contains the signed modulus of the current density vector function 

divergence calculation (``calc=divj``):

* ``divj.txt``: one line per grid point, in the order of the grid points: x,y,z in bohr and the divergence of the current density in atomic units
* ``divj.vti``: the divergence of the current density, written for evenly spaced grids
* ``divj.vtu``: the divergence of the current density on the tetrahedral mesh ``grid.1.ele``, written for grids read from a file and for Gauss grids if that file is present

For open-shell calculations the files ``divjalpha``, ``divjbeta`` and ``divjspindens`` are written in addition.

only for open-shell calculations:

* ``jmodalpha.txt``: contains the signed modulus of the alpha current density vector function 
* ``jmodalpha.vti``: contains the signed modulus of the alpha current density vector function 
* ``jmodbeta.txt``: contains the signed modulus of the beta current density vector function 
* ``jmodbeta.vti``: contains the signed modulus of the beta current density vector function 
* ``jmodspindens.txt``: contains the signed modulus of the difference (alpha - beta) current density vector function 
* ``jmodspindens.vti``: contains the signed modulus of the difference (alpha - beta) current density vector function 

* ``jvecalpha.vti``: contains information about the alpha current density vector function 
* ``jvecbeta.vti``: contains information about the beta current density vector function 
* ``jvecspindens.vti``: contains information about the difference (alpha - beta) current density vector function 

Getting the numbers out of a .vti file: vti2csv
------------------------------------------------

``vti2csv`` writes the data of a ".vti" file as a CSV table with one row per
grid point, which can be read by a spreadsheet, gnuplot or pandas::

   vti2csv jvec.vti              # writes jvec.csv
   vti2csv jvec.vti -o out.csv
   vti2csv jvec.vti -o -         # writes to standard output

The columns are ``x,y,z,jx,jy,jz`` for a file with vectors (``jvec.vti``) and
``x,y,z,value`` for a file with one number per point (``jmod.vti``,
``acid.vti``, ``divj.vti``, ``tropicity.vti``). The coordinates are those of
the ".vti" file, in bohr, with x running fastest; the data are in the units of
the file. A ".vti" file describes a grid along the Cartesian axes, so the
coordinates are those of the grid points only if the basis vectors of the grid
point along x, y and z; for a rotated grid take the positions from the grid
definition instead (see the note in :doc:`visualization`). Both binary and text ".vti" files are read. The command needs only
Python and is installed next to ``gimic``.

A CSV file is several times larger than the binary ".vti" file, so convert
only the files you need as a table.
