# CMakeLists.txt for event package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(FitMomentumScale)

add_compile_options(-g -Wall -Wextra -pedantic -Werror)

### where to put the libraries
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
### where to put the executables
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/include/ )

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
set(ROOT_CONFIG_DEBUG 1)
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS RIO Core Net RooFit RooFitCore MathCore Hist Spectrum)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})


include_directories(${CMAKE_CURRENT_SOURCE_DIR})
ROOT_GENERATE_DICTIONARY(G__FitMomentumScale RooTwoSidedCBShape.h LINKDEF LinkDef.h OPTIONS -p)
ROOT_GENERATE_ROOTMAP( FitMomentumScaleDict LINKDEF RooTwoSidedCBShape.h  DEPENDENCIES )
#ROOT_LINKER_LIBRARY( FitMomentumScaleDict G__FitMomentumScale.cxx RooTwoSidedCBShape.cxx  HelperClasses.cxx LIBRARIES RIO Core RooFit RooFitCore MathCore)


#---Create a shared library with geneated dictionary
add_library(FitMomentumScale SHARED G__FitMomentumScale.cxx RooTwoSidedCBShape.cxx  HelperClasses.cxx )
target_link_libraries(FitMomentumScale ${ROOT_LIBRARIES})


#---Create  a main program using the library
add_executable(performFit TestFit.C)
target_link_libraries(performFit FitMomentumScale)

