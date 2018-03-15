#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(DuneDM)

set(EXROOTANALYSIS_PATH "" CACHE FILEPATH "Path to exrootanalysis")
set(DK2NU_PATH "" CACHE FILEPATH "Path to dk2nu")

find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})

include_directories( ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS} "${EXROOTANALYSIS_PATH}/ExRootAnalysis" "${DK2NU_PATH}/tree" )
link_directories( ${ROOT_LIBRARY_DIR} ${EXROOTANALYSIS_PATH} "${DK2NU_PATH}/build/lib" )

ROOT_GENERATE_DICTIONARY("${CMAKE_CURRENT_BINARY_DIR}/DuneDMDict.pcm" "${EXROOTANALYSIS_PATH}/ExRootAnalysis/ExRootClasses.h" "${DK2NU_PATH}/tree/dk2nu.h" LINKDEF "${PROJECT_SOURCE_DIR}/LinkDef.h" )

file(GLOB DuneDM_SRCS ${PROJECT_SOURCE_DIR}/*.cpp)

#add_definitions(${ROOT_DEFINITIONS})
add_executable(DuneDM ${DuneDM_SRCS} "${CMAKE_CURRENT_BINARY_DIR}/DuneDMDict.pcm")
target_link_libraries(DuneDM ${ROOT_LIBRARIES} libExRootAnalysis.so libarmadillo.so libdk2nuTree.so)


