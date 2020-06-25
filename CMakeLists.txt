###########################################################################
#                 CMake Build File for ICM-HYDRO
#
#    Written By: Zach Cobell
#
###########################################################################
#
# The CMake build system enable ICM-HYDRO to be deployed and built
# in a cross platform environment. 
#
###########################################################################
#...Set the default build type
IF(DEFINED CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of
        build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug
        Release RelWithDebInfo MinSizeRel.")
ELSE()
    SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build,
        options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release
        RelWithDebInfo MinSizeRel.")
ENDIF()
###########################################################################

###########################################################################
#  GENERAL OPTIONS
###########################################################################
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)
PROJECT(ICM_HYDRO Fortran)
ENABLE_LANGUAGE( C )

###########################################################################


###########################################################################
#...Put the static libraries in the CMakeFiles folder so they don't
#   contaminate the build directory
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/CMakeFiles)
###########################################################################


###########################################################################
#...Set the version string (used for libraries)
SET(ICMHYDRO_VERSION_MAJOR 23)
SET(ICMHYDRO_VERSION_MINOR  1)
SET(ICMHYDRO_VERSION_PATCH  3)
SET(ICMHYDRO_VERSION_STRING ${ICMHYDRO_VERSION_MAJOR}.${ICMHYDRO_VERSION_MINOR}.${ICMHYDRO_VERSION_PATCH})
###########################################################################


###########################################################################
#...Determine architecture specific parameters
INCLUDE(${CMAKE_SOURCE_DIR}/cmake/architecture.cmake)
###########################################################################


###########################################################################
# ICM_HYDRO
###########################################################################
ADD_EXECUTABLE( hydro_v${ICMHYDRO_VERSION_STRING} ChLaP.f DIN.f DIP.f DON.f 
                                                 DOP.f DOx.f 
                                                 ICM_InterpolateMonthlyToGrid.f 
                                                 ICM_InterpolateToGrid.f 
                                                 ICM_MapMonthlyToGrid.f 
                                                 ICM_MapToGrid.f 
                                                 ICM_SalinityThreshold_WM.f 
                                                 ICM_TreeConditions_Veg.f 
                                                 ICM_formatting.f 
                                                 ICM_summaries.f 
                                                 NH4.f NO3.f POP.f TOC.f 
                                                 TP.f TSSolids.f allocate_params.f 
                                                 celldChem.f celldQ.f celldSS.f 
                                                 celldSal.f celldTmp.f chemical.f 
                                                 deadA.f hydrod.f infile.f livA.f 
                                                 main.f orgN.f params.f photo.f 
                                                 salinity.f temperature.f tideBC.f 
                                                 vanRijnSediment.f waves_YV.f )
SET_TARGET_PROPERTIES(hydro_${ICMHYDRO_VERSION_STRING} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/CMakeFiles/mod)
SET_TARGET_PROPERTIES(hydro_${ICMHYDRO_VERSION_STRING} PROPERTIES COMPILE_FLAGS "${Fortran_LINELENGTH_FLAG} ${Fortran_COMPILER_SPECIFIC_FLAG}" )
