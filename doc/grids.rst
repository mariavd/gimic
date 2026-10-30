.. role:: math(raw)
   :format: html latex


Grids
=====

There are three possible grid topologies in gimic: 'std (or base)' which is
3-dimensional (although the extend in one direction may be zero), 'bond' which
is 2-dimensional, and 'file' which is an arbitrary set of points read from a
file. 'std' grids are especially useful for visualisations of an entire current
vector field. 'bond' grid are used to integrate the current that flows through
the grid.

For an 'std' or 'bond' grid one can specify the spacing (type) of points in
each dimension. The choices are 'even' for an equidistant grid, 'gauss' for a
Gauss grid and 'lobatto' for a Gauss-Lobatto grid.  The recommended choices are
'even' for 'std' and 'gauss' for 'bond'.  When a quadrature grid is specified
the order of the quadrature must also be specified with the ’gauss\_order’
keyword. The number of grid points in each direction is specified either
explicitly using either of the array keywords ’grid\_points’ or ’spacing’. If
the chosen grid is not a simple even spaced grid, the actual number of grid
points will be adjusted upwards to fit the requirements of the chosen
quadrature.

The shape of the grid can also be modified by the ’radius’ key, which
specifies a cutoff radius. This can be useful for integration. Sometimes
it’s practical to be able to specify a grid relative to a well know
starting point. The ’rotation’ keyword specifies Euler angles for
rotation according to the x->y->z convention. Note that the magnetic
field is not rotated, unless it is specified with ’magnet\_axis=i,j or
k’.

GIMIC automatically outputs a number of .xyz files containing dummy
points to show how the grids defined actually are laid out in space.

Basic grids
-----------

The ’std’ grid is defined by giving an ’origin’ and two orthogonal basis
vectors ’ivec’ and ’jvec’ which define a plane. The third axis is
determined from :math:`\vec k=\vec i\times\vec j`. The array ’lengths’
specifies the grid dimensions in each direction.

grid(std)
~~~~~~~~~

::

  type=even
  origin=-8.0, -8.0, 0.0
      Origin of grid
  
  ivec=1.0, 0.0, 0.0
      Basis vector i
  
  jvec=0.0, 1.0, 0.0
      Basis vector j ( k = i x j )
  
  lengths=16.0, 16.0, 0.0
      Lenthts of (i,j,k)
  
  spacing=0.5, 0.5, 0.5
      Spacing of points on grid (i,j,k)
  
  grid\_points=50, 50, 0
      Number of gridpoints on grid (i,j,k)
  
  rotation=0.0, 0.0, 0.0
      Rotation of (i,j,k) -> (i’,j’,k’) in degrees. Given as Euler angles
      in the x->y->z convention.

Bond grids
----------

The ’bond’ type grids define a plane through a bond, or any other
defined vector. The plane is orthogonal to the vector defining the bond.
The bond can be specified either by giving two atom indices,
’bond=[1,2]’, or by specifying a pair of coordinates, ’coord1’ and
’coord2’. The position of the grid between two atoms is determined by
the ’distance’ key, which specifies the distance from atom 1 towards
atom 2. For analysing dia- and paramagnetic contributions, the positive
direction of the bond is taken to be from atom 1 towards atom 2. Since
one vector is not enough to uniquely defining the coordinate system
(rotations around the bond are arbitrary), a fixpoint must be specified
using either the ’fixpoint’ atom index or the ’fixcoord’ keyword. This
triple of coordinates is also used to fix the direction of the magnetic
field when the ’magnet\_axis=T’ is used.

The shape and size of the bond grid is specified by the keywords 'width' and
'height'.  Each of the four components is relative to a point on the line
connecting the two reference atoms/coordinates.  Typically (but not
necessarily) the first component of both ranges is negative and the second
component is positive.

::

    width=[-1.5, 5.0]
    height=[-5.0, 5.0]

grid(bond):
~~~~~~~~~~

type=gauss\|lobatto
    Use uneven distribution of grid points for quadrature

bond=1,2
    Atom indices for bond specification

fixpoint=5
    Atom index to use for fixing the magnetic field and grid orientation

coord1=0.0, 0.0, 3.14
    Coordinate of atom 1

coord2=0.0, 0.0, -3.14
    Coordinate of atom 2

fixcoord=0.0, 0.0, 0.0
    Fixation coordinate

distance=1.5
    Place grid ’distance’ between atoms 1 and atom 2

gauss\_order=9
    Order for Gauss quadrature

spacing=0.5, 0.5, 0.5
    Spacing of points on grid (i,j,k) (approximate)

grid\_points=50, 50, 0
    Number of grid points on grid (i,j,k) (approximate)

height=-4.0, 4.0
    Grid size relative to grid center

width=-1.0, 6.0
    Grid size relative to grid center

radius=3.0
    Create a round grid by cutting off at radius

rotation=0.0, 0.0, 0.0
    Rotation of (i,j,k) -> (i’,j’,k’) in degrees. Given as Euler angles
    in the x->y->z convention.

towards=1, 2, 3, 4, 5, 6
    Reference centre: the centroid of these atoms (typically the ring the
    integration starts from). Two effects: for bond grids it is the default
    fixed point, so that the integration starts on its side of the bond, and
    the sign of the integrated current is chosen such that a current flowing
    in the diatropic sense around it is positive, whatever the order of the
    bond atoms or the sign of the field. See `Current profiles`_.

towards_coord=0.0, 0.0, 0.0
    The reference centre as a point instead of atoms.

slices=0.02, 0.0
    Current profile: slice thickness along the width and along the height
    of the plane in bohr (0 = not sliced). See `Current profiles`_.

Current profiles
----------------

The current through a plane can be resolved into slices in one run. The
quadrature grid of GIMIC is a composite Gauss rule, i.e. a sequence of blocks
of ``gauss_order`` points each, so with ``slices`` the number of grid points
along a sliced axis is set to one block per slice (``spacing`` and
``grid_points`` are ignored for that axis) and the current through each
slice is the sum over its block of the column integrals, with no extra cost.
A profile with 500 slices of 0.02 bohr is therefore one calculation over the
whole plane instead of 500 calculations, and it is parallel with OpenMP over
the columns of the grid (set ``OMP_NUM_THREADS``; the ``gimic`` launcher
keeps BLAS single-threaded).

::

    magnet=[0.0, 0.0, -1.0]
    Grid(bond) {
        type=gauss
        bond=[1,2]
        towards=[1,2,3,4,5,6]     # the ring: fixed point and sign convention
        distance=1.3147
        gauss_order=9
        spacing=[0.025, 0.0, 0.0] # along the height; the width is set by the slices
        height=[-8.0, 8.0]
        width=[-2.6, 8.0]         # from the ring centre to 8 bohr outside
        slices=[0.02, 0.0]        # 530 slices along the width
    }

The sliced axis must use ``type=gauss`` or ``lobatto``. ``slices=[0.0, d]``
gives a vertical profile, ``slices=[dw, dh]`` a map of rectangles. Files
written: ``current_profile.dat`` (one line per slice, from the start of the
integration: distance, net, positive and negative current in nA/T; the format
of the current-profile scripts), ``vertical_profile.dat`` (the same for the
slices along the height), ``profile.dat`` (the same with the geometry,
the slice centres in space and the current per unit length at every Gauss
node), and ``squares_profile.dat`` when both axes are sliced (a matrix, rows
along the height and columns along the width, followed by a list with
coordinates). For open-shell calculations the files carry the suffixes
``_alpha``, ``_beta`` and ``_spindens``.

Sign convention. The integral is the flux of the current density through the
plane along its normal, which for a bond grid is the direction from atom 1 to
atom 2: with an explicit ``magnet`` vector the sign therefore flips with the
order of the two atoms, and the fixed point only decides on which side of the
bond the integration starts. With ``towards`` the normal is oriented such that
a current circulating around the reference centre in the diatropic sense is
positive, evaluated at the end of the profile (the "out" edge of a bond grid,
the edge far from the origin of a base grid); on the far side of the centre
the same circulation crosses the plane the other way and appears negative,
as it should. With ``magnet_axis=X`` the field follows the plane and the sign
is independent of the atom order already.

A profile over a plane in an arbitrary position is a base grid: ``origin`` is
one corner, ``ivec`` the edge along the field, ``jvec`` the edge along the
profile, ``lengths`` the two edge lengths, ``type=gauss`` and
``slices=[d, 0.0]``; the profile then runs from the origin along ``jvec``.

