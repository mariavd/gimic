
# GIMIC Quick-start Guide

In order to perform your first GIMIC calculation, you need to obtain the unperturbed and the magnetically perturbed density matrices in a quantum chemistry program. Refer to [url]https://gimic.readthedocs.io/en/latest/usage.html[/url] for more details. Assuming that the ``XDENS`` and ``MOL`` files containing the density matrices and the basis set data respectively are available, you can start the current density analysis. 

[image: gimic-steps.jpg]

Current density analysis should be done both visually and quantitatively. Start by doing a visual analysis of the current density field. After choosing the current pathways of interest, you can perform integrations of the current strength. 

[image: gimic-analysis.jpg]

## Visual analysis 

A cubical 3D grid is defined around the molecule such that the box extends about 8 bohr away from the atoms in the plane perpendicular to the magnetic field, and about 4 bohr above and below the molecule. 


**Problem 1.** I would like to calculate the current profile for the benzene molecule passing the bond between two carbon atoms. 

I will employ the bond type grid. I have to specify the indices of the atoms through which the plane will cross. They are specified by the keyword ``bond=[a1,a2]`` where a1 and a2 are the indices of the atoms according to the XYZ file. NOTE: that the counting starts from 1. 

Any plane is defined by three points. Two of them are the atomic coordinates. The third point will define one of the infinitely many points which pass through the chemical bond. A physically meaningful choice is to specify the plane which is perpendicular to the benzene ring. We call this point ``fixpoint``. 

Typically, the plane is placed such that it passes through the midpoint of the chemical bond. This is defined through the ``distance`` keyword. With this the integration plane is fully specified. 

However, we cannot perform calculations on an infinitely big integration plane. We need to choose a rectangle on the plane such that the current density vanishes from all sides. We call ``height`` the distance above and below the molecule. This is the side of the rectangle perpendicular to the benzene ring. The other side of the rectangle is parallel to the molecular ring. It is called ``width``.

In the minimal example, we also need to provide the ``spacing`` between the grid points on the integration plane. Three numbers are required, one for each spatial dimension. The numberical integration is performed using Gauss quadrature. This is specified with the keyword ``type=gauss``. The order of the Gauss quadrature is given as ``gauss_order=9``. 

Minimal example for the definition of a bond grid:

```
Grid(bond) {                    # define grid orthogonal to a bond 
    type=gauss                  # gauss distribution of grid points for the integration
    bond=[1,2]                      # the two atoms
    fixpoint=4                  # the third point defining the plane
    distance=1.32               # place grid 'distance' between atoms
    
    gauss_order=9               # order for gauss quadrature
    
    height=[-5.0, 5.0]
    width=[-2.2, 5.0]
    
    spacing=[0.02, 0.02, 0.02]     # spacing of points on grid (i,j,k)
}
```


