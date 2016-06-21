cmake_minimum_required(VERSION 2.8)
project(barneshut)

# version number
set (bh_VERSION_MAJOR 1)
set (bh_VERSION_MINOR 0)

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
        "${PROJECT_SOURCE_DIR}/src/OctTree/dendro.h.in"
        "${PROJECT_BINARY_DIR}/dendro.h"
)

if(OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

if(MPI_COMPILE_FLAGS)
  set(COMPILE_FLAGS "${COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set(LINK_FLAGS "${LINK_FLAGS} ${MPI_LINK_FLAGS}")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(LINK_FLAGS "${LINK_FLAGS} -rdynamic")

option(USE_64BIT_INDICES "Use 64-Bit indices. Reverts to 32-bit if turned off" ON)
option(ALLTOALLV_FIX "Use K-way all to all v" ON)
option(HILBERT_ORDERING "Use Hilbert ordering instead of Morton" ON)
option(ENERGY_CONSERVATION "Compute total energy before and after" OFF)
option(__SILENT_MODE__ "Restrict output (used for timings)" OFF)
option (SPLITTER_SELECTION_FIX "Use the splitter fix for the treeSort" ON)
option (DIM_2 "To enable DIM2 version of Sorting. Tree sort part is tested and works wioth DIM 2 but rest of the dendro might not " OFF)
option (REMOVE_DUPLICATES "Enables the remove duplocates and multilevel Sorting. " ON)
set(KWAY 128 CACHE INT 128)
set(NUM_NPES_THRESHOLD 16 CACHE INT 16)

if(USE_64BIT_INDICES)
    add_definitions(-DUSE_64BIT_INDICES)
endif()

if(ALLTOALLV_FIX)
    add_definitions(-DALLTOALLV_FIX)
    add_definitions(-DKWAY=${KWAY})
endif()

if(HILBERT_ORDERING)
    add_definitions(-DHILBERT_ORDERING)
endif()

if(ENERGY_CONSERVATION)
  add_definitions(-DENERGY_CONSERVATION)
endif()

if(__SILENT_MODE__)
  add_definitions(-D__SILENT_MODE__)
endif()

if(SPLITTER_SELECTION_FIX)
    add_definitions(-DSPLITTER_SELECTION_FIX)
    add_definitions(-DNUM_NPES_THRESHOLD=${NUM_NPES_THRESHOLD})
endif()

if(DIM_2)
    add_definitions(-DDIM_2)
endif()

if(REMOVE_DUPLICATES)
    add_definitions(-DREMOVE_DUPLICATES)
endif()

##------
include_directories(${PROJECT_BINARY_DIR}
                    ${MPI_INCLUDE_PATH}
		    src
		    src/BarnesHut
		    src/OctTree
		    src/OctTree/tree
		    src/OctTree/test
                    )

add_executable(barneshut src/BarnesHut/barneshut.cpp
                         src/BarnesHut/planet.h
			 src/BarnesHut/msd.h
			 src/BarnesHut/OctTree.h
			 src/BarnesHut/EnergyConservation.h
			 src/BarnesHut/vector.h
			 src/BarnesHut/bh_hermite.h
			 src/BarnesHut/hermite.h
			 src/BarnesHut/barneshut.h
			 src/OctTree/binUtils.h
			 src/OctTree/Point.h
			 src/OctTree/sfcSort.h
			 src/OctTree/hcurvedata.h
			 src/OctTree/TreeNode.h
			 src/OctTree/treenode2vtk.h
			 src/OctTree/parUtils.h
			 src/OctTree/ompUtils.h
			 src/OctTree/dtypes.h
			 src/OctTree/colors.h
			 src/OctTree/test/testUtils.h
			 src/OctTree/seqUtils.h
			 src/OctTree/nodeAndValues.h
			 src/OctTree/octUtils.h
			 src/OctTree/tree/treenode2vtk.cpp
			 src/OctTree/tree/binUtils.cpp
			 src/OctTree/tree/parUtils.cpp
			 src/OctTree/tree/TreeNode.cpp
			 src/OctTree/tree/octFunctions.cpp
			 src/OctTree/tree/hcurvedata.cpp
                         )
target_link_libraries(barneshut ${MPI_LIBRARIES} m)

add_executable(generate src/BarnesHut/generate.cpp)
target_link_libraries(generate ${MPI_LIBRARIES} m)