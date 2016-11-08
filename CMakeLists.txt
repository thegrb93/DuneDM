#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(DuneDM)
set(CMAKE_BUILD_TYPE Debug)
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")

set(EXROOTANALYSIS_PATH "" CACHE FILEPATH "Path to exrootanalysis")

find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})

include_directories( ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS} "${EXROOTANALYSIS_PATH}/ExRootAnalysis" )
link_directories( ${ROOT_LIBRARY_DIR} ${EXROOTANALYSIS_PATH} )

ROOT_GENERATE_DICTIONARY("${CMAKE_CURRENT_BINARY_DIR}/DuneDMDict.pcm" "${EXROOTANALYSIS_PATH}/ExRootAnalysis/ExRootClasses.h" LINKDEF "${PROJECT_SOURCE_DIR}/LinkDef.h" )

file(GLOB DuneDM_SRCS ${PROJECT_SOURCE_DIR}/*.cpp)

add_definitions(${ROOT_DEFINITIONS})
add_executable(DuneDM ${DuneDM_SRCS} "${CMAKE_CURRENT_BINARY_DIR}/DuneDMDict.pcm")
target_link_libraries(DuneDM ${ROOT_LIBRARIES} libExRootAnalysis.so)


