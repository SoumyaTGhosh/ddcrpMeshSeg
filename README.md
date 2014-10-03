#ddcrpMeshSegMentation: 
###Distance Dependent Chinese Restaurant Processes for Mesh Segmentation
####Author: Soumya Ghosh (sghosh <AT> cs.brown.edu)
------------------------------------
**Mesh segmentation software**. This software package includes an implementation of Gibbs sampling for the 
distance dependent Chinese restaurant process along with utility functions for motion based mesh segmentation.
The code is based on the model presented in:

**From Deformations to Parts:Motion-based Segmentation of 3D Objects**.
Soumya Ghosh, Erik B. Sudderth, Matthew Loper and Michael J. Black.
NIPS 2012. 

Please cite this paper if you use this software. Given meshes in different poses the model determines both the *number* and *extent* of segments. Additionally, it *guarantees* spatially connected segments.

**Dependencies -- You will need the following packages**:

1) Tom Minka's lightspeed toolbox:
http://research.microsoft.com/en-us/um/people/minka/software/lightspeed/ 

2) David Gleich's Matlab interface to the Boost Graph Library:
http://www.cs.purdue.edu/homes/dgleich/packages/matlab_bgl/
