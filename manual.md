# FEM_ORDER2 manual
*****

TABLE OF CONTENTS

- [FEM_ORDER2 manual](#fem_order2-manual)
- [Overview](#overview)
- [Numerical Methods](#numerical-methods)
  - [Constructing The Linear Elastic Problem](#constructing-the-linear-elastic-problem)
- [Module Brief](#module-brief)
    - [**fem_order2** in fem_order2.f90](#fem_order2-in-fem_order2f90)
      - [*Partition Data:*](#partition-data)
      - [*Cell Partition Data:*](#cell-partition-data)
      - [*Partitioned Mesh:*](#partitioned-mesh)
      - [*Solver and Solution Data:*](#solver-and-solution-data)
      - [*Scatterers:*](#scatterers)
      - [*Auxiliary Data:*](#auxiliary-data)
      - [*FEM Fundamentals:*](#fem-fundamentals)
      - [subroutine **initializeStatus**](#subroutine-initializestatus)
      - [subroutine **initializeLib**](#subroutine-initializelib)
      - [subroutine **initializeConstitution**](#subroutine-initializeconstitution)
      - [subroutine **ReducePoints**](#subroutine-reducepoints)
      - [subroutine **getVolumes**](#subroutine-getvolumes)
      - [subroutine **output_plt_mesh**](#subroutine-output_plt_mesh)
      - [subroutine **output_plt_scalar**](#subroutine-output_plt_scalar)
      - [subroutine **output_plt_scalar3x**](#subroutine-output_plt_scalar3x)
      - [subroutine **write_binary_plt_string**](#subroutine-write_binary_plt_string)
      - [subroutine **SetUpPartition**](#subroutine-setuppartition)
      - [subroutine **DoCellPartition**](#subroutine-docellpartition)
      - [subroutine **DoGeomPartition**](#subroutine-dogeompartition)
      - [subroutine **SetUpThermalBC_BLOCKED**](#subroutine-setupthermalbc_blocked)
      - [subroutine **SetUpThermal_InitializeObjects**](#subroutine-setupthermal_initializeobjects)
      - [subroutine **SetUpThermalPara**](#subroutine-setupthermalpara)
      - [subroutine **SetUpElasticBC_BLOCKED**](#subroutine-setupelasticbc_blocked)
      - [subroutine **SetUpElasticity_InitializeObjects**](#subroutine-setupelasticity_initializeobjects)
      - [subroutine **SetUpElasticityPara**](#subroutine-setupelasticitypara)
      - [subroutine **SolveThermal_Initialize**](#subroutine-solvethermal_initialize)
      - [subroutine **SolveThermal**](#subroutine-solvethermal)
      - [subroutine **SolveElasticity_Initialize**](#subroutine-solveelasticity_initialize)
      - [subroutine **SolveElasticity**](#subroutine-solveelasticity)
      - [subroutine **GetElasticityUGradient**](#subroutine-getelasticityugradient)
      - [subroutine **GetStrainStress**](#subroutine-getstrainstress)
      - [subroutine **GatherVec1**](#subroutine-gathervec1)
      - [subroutine **GatherVec3**](#subroutine-gathervec3)

# Overview

The current program is designed to solve the problem of heat transfer and heat-related elastic problem.

The method is traditional FEM performed on a 2nd-order mesh, with nodes in edges but not on in faces or volumes. The mesh is read into the program SEQUENTIALLY, where only process is responsible for that. And all the definitions of the problem, including boundary conditions an constitutional relations are also only given to process 0. 

After reading mesh and problem definitions, mesh info is partitioned and distributed to all processes, both topology and needed point coordinates. Meanwhile, the size of needed CSR is counted and distributed so that preallocation for PETSC matrix can be done. Then, to solve the FEM problems, the program conceptually performs integration, both in the volumes and faces, to obtain the stiffness matrix and right hand side vector (RHS). The integration is performed within the normalized space for each face or volume element. As we are using isoparametric method, the jacobian for transformation between physical space to the element's normalized space can be calculated via some simple linear algebraic techniques. In the end, the program uses PETSC's KSP object to solve the linear system. 

In the current implementation, while assembling the matrices and RHS, only volume integration is well parallelized while surface integration on boundary conditions is performed only by process 0. When the boundary conditions are relatively small in number, performance is unchanged. I have future plans of parallization of the boundary integrals.

# Numerical Methods

In conventional finite-element method, or a common Galerkin method, the trial functions, being identical with the bases of the discrete solution, are considered to be $C^0$. Therefore the piecewise defined polynomial bases should maintain $C^0$ continuity on the interfaces of volumes (take 3-D for example). As a result, it is better to consider the discrete DOFs to be set on the vertices, and the bases are interpolation functions that satisfy a kronecker-delta property over the vertices. For example, when each volume is a tetrahedron with 4 vertices, you can easily transform it linearly into a corner of a cartesian box, and using 3 cartesian axes you can easily define first-order bases functions for each vertex. Generally, the basis functions are defined in a normalized coordinate system as polynomials. As the transformation between the normalized space and the geometric space is generally not a linear mapping (basically as curved elements), generally the same set of basis functions are used to interpolate the mapping. Fundamentally, the actual bases are fractions rather than polynomials, but as flat-faced and near flat-faced elements are the majority, the orders of polynomial-based numerical integral are mostly decided with the situation of a linear spacial mapping.

This program uses serendipity method to infer the correct basis functions in the normal space, which is 2nd order polynomials. The exact formulae can be found in the source code.

Using the interpolation described above, a Galerkin method can be derived. For linear heat transfer and linear elasticity, the Galerkin discretion causes the ODEs: 

$$
M_t\ddot{a}_t+C_t\dot{a}_t+K_ta_t=F_t
$$

$$
M\ddot{a}+C\dot{a}+Ka=F
$$

Where the subscripts 't' denote that it's related to the heat transfer problem, or elastic problem otherwise. $a$ and $a_t$ are just DOFs or nodal values, in the elastic problems, we denote that for nodes $[i_0,i_1,...]$, $a$ is column vector $[u_0,v_0,w_0,u_1,v_1,w_1,...]^T$, where $[u,v,w]^T$ is the displacement of some point.
In common problems, only $M,M_t,K,K_t,F,F_t$ are desired, so we discuss them here first:

## Constructing The Linear Elastic Problem

The Galerkin method tells us that: 

$$
K=\int_{\Omega}{B^{T}DBdV}+\int_{\partial \Omega}{N^THNd\Gamma},\ \ \ \;
M=\int_{\Omega}{N^{T}\rho NdV},\ \ \ \;\\\ \\
F=\int_{\Omega}{f^TNdV}+\int_{\partial \Omega}{p^TNd\Gamma}+\int_{\partial \Omega}{(Hx_0)^TNd\Gamma}
$$

$\Omega$ is the domain of definition, and $dV$ denotes its differential. $\partial \Omega$ is the boundary of the domain of definition, and $d\Gamma$ denotes its differential. 

In the formulae above, $B$ is the function transforming $a$ into the interpolated strain field (using a as the right vector), where $B(x,y,z)_{ij}a_j=e(x,y,z)_i$, and $i=1,2,...,6$ is the subscript for the components of the strain vector, where the transposing and matrix products are on the $i,j$ subscripts.

$D=D_{ij}$ is the 6x6 constitutional relation transforming strain into stress, where $\sigma_i=D_{ij}e_j$.

$N$ transforms $a$ into the interpolated displacement field, it's definition is similar with B, where $N(x,y,z)_{ij}a_j=disp(x,y,z)_i\equiv[u,v,w]^T(x,y,z)_i$.

$H$ is the 3x3 matrix defining the stiffness of a linearly elastic basement, (not a scalar instead of a matrix for the stiffness could very likely be anisotropic), transforming displacement form base-point $x_0$ to facial force. $H$ is zero on the boundaries not defined as a elastic basement. $H$ should be symmetric and positive semidefinite to maintain correct physical meaning, that is, it is actually defined by its 3 eigenvalues and 3 eigenvectors.

$\rho$ is the density field, defined as a scalar field.

$f$ is the volume force, including the effect of prestress caused by temperature change.

$p$ is only not zero on known-force boundaries, it is a 3d vector field.

$Hx_0$ can be viewed as the payload force caused by elastic basement boundaries. Actually, in the program, force boundaries are performed with the elastic base boundary scheme. As the inputs are $Hx_0$ and $H$, setting $H=0$ (while $x_0\rarr\infty$, which is not the numeric input) means the boundary condition becomes a given force.

The discrete values of $K,M,F$ are integrated on the piecewise polynomials, where (including $N$, $B$) they can actually be decomposed into each element (of volume or face), (for 3d case) we define $K_e,M_e,F_e$, which are the results of integrals performing only in the volume $e$. More over, the ordering of $K_e,M_e,F_e$ are solely consisting the part of $a$ (DOFs on nodes) adjacent to the volume (the other parts are zero in the global matrices and vectors). Therefore, 


# Module Brief

### **fem_order2** in fem_order2.f90

This module is the primary container of all the data and functions directly connected to mesh management, case building, solving and output functions.

#### *Partition Data:*

Created for initial partitioning of the mesh, primarily created in **SetUpPartition**.

#### *Cell Partition Data:*

Created with the help of Partition Data, primarily created in **DoCellPartition**. Used for distributing the cell data.

#### *Partitioned Mesh:*

Distributed mesh, stored in CSR. In the columns of the CSR, denoting the index INODE, if INODE is in [0,localDOFs), then the vertex is a local vertex, its global index is INODE+indexLo. If INODE is in [localDOFs,inf), then the vertex is a ghost vertex, its global index is ghostingGlobal(INODE-localDOFs+1). To acquire nodal data on the vertex, the index INODE(0 based) is used in the local array of Vec objects, and the Vec objects should be equipped with ghost points. For example, localCoords is a 3x Vec object equipped with 3x ghost points, that is, all the elements are expanded into blocks in 3, which means INODE indicates (INODE+1)\*3-2~(INODE+1)\*3 (1 based) data in the local array of localCoords. This procedure is demonstrated in the SetUpThermalPara and SetUpElasticityPara.

#### *Solver and Solution Data:*

Solver and solution data, along with boundary-condition definitions.

Vec dofFix\<T\>Dist is distributed 1x or 3x nodal data, where if the data is not nan, the nodal dof is mathematically fixed to this value.

#### *Scatterers:*

Scatterer objects, to redistribute nodal data between: [1] original data numbering on proc0 and [2] partitioned data numbering on all procs.

Created automatically when calling GatherVec1 or GatherVec3

#### *Auxiliary Data:*

Data used temporarily but could be used in some occasion afterwards. Need a clean-up procedure. 

localAdjacencyNum and ghostAdjacencyNum are merely indicating the row sizes of the parallel matrix to minimize dynamic memory allocation and space wasting.

#### *FEM Fundamentals:*

Some tiny data predefined for FEM-related calculation.

#### subroutine **initializeStatus**

To initialize status of some auto-created objects.

Call collectively.

#### subroutine **initializeLib**

To initialize tiny data of element libs, concerning the shape functions and integration points in a regularized coordinate.

Call collectively.

#### subroutine **initializeConstitution**

To initialize tiny data concerning constitutional relations.

Call collectively.

#### subroutine **ReducePoints**

To reduce the unused points in mesh data, actually only concerning module globals.

Call on proc0.

#### subroutine **getVolumes**

To calculate cell volumes, useful in checking if input cells are bad. Only on proc0.

Call on proc0.

#### subroutine **output_plt_mesh**

To write the mesh file, with data in module globals.

Call on proc0.

#### subroutine **output_plt_scalar**

To output a scalar field, either cell-centered or nodal, in original numbering.

Call on proc0.

#### subroutine **output_plt_scalar3x**

Same as **output_plt_scalar**, only that data is blocked in 3.

Call on proc0.

#### subroutine **write_binary_plt_string**

Internally called, to write strings in tecplot scheme.

#### subroutine **SetUpPartition**

Creates partition and distributes the data to all procs.

Call collectively.

#### subroutine **DoCellPartition**

Internally called, to distribute cell data to processes.

#### subroutine **DoGeomPartition**

Internally called, to distribute node coords to processes.

#### subroutine **SetUpThermalBC_BLOCKED**

To preprocess the block-defined boundary conditions. Currently only distributes the fixed boundary condition to procs.

Call collectively.

#### subroutine **SetUpThermal_InitializeObjects**

To initialize load vector and stiffness matrix objects for thermal, and do preallocation for the matrix.

Call collectively.

#### subroutine **SetUpThermalPara**

To assemble (integrate) the stiffness matrix and the load vector. Basically doing volume and surface integration.

Call collectively.

#### subroutine **SetUpElasticBC_BLOCKED**

Same as **SetUpThermalBC_BLOCKED** but the elastic part.

Call collectively.

#### subroutine **SetUpElasticity_InitializeObjects**

Same as **SetUpThermal_InitializeObjects** but the elastic part.

Call collectively.

#### subroutine **SetUpElasticityPara**

Same as **SetUpThermalPara** but the elastic part.

Call collectively.

#### subroutine **SolveThermal_Initialize**

To initialize and configure solver objects for thermal problem.

Call collectively.

#### subroutine **SolveThermal**

To actually solve the problem.

Call collectively.

#### subroutine **SolveElasticity_Initialize**

Same as **SolveThermal_Initialize** but the elastic part.

Call collectively.

#### subroutine **SolveElasticity**

Same as **SolveThermal** but the elastic part.

Call collectively.

#### subroutine **GetElasticityUGradient**

To calculate the gradient of displacement for a elastic solution.

Call collectively.

#### subroutine **GetStrainStress**

To calculate the strain and stress using the gradient of displacement for a elastic solution.

Call collectively.

#### subroutine **GatherVec1**

To redistribute nodal data between: [1] original data numbering on proc0 and [2] partitioned data numbering on all procs, using the scatterers.

If ifreverse is true, then scatters [1] to [2], else gathers [2] to [1].

Call collectively.

#### subroutine **GatherVec3**

The 3-blocked version of **GatherVec3**.

Call collectively.

