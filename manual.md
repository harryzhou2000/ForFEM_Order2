- [FEM_ORDER2 manual](#fem_order2-manual)
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

# FEM_ORDER2 manual

## Module Brief

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

