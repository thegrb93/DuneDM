#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(DuneDM)
set(CMAKE_BUILD_TYPE Debug)
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")

find_package(ROOT REQUIRED)


include_directories( ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIR} )
link_directories( ${ROOT_LIBRARY_DIR} )

#include(${ROOT_USE_FILE})
#include(CMakeParseArguments)


#----------------------------------------------------------------------------
#
# libDuneDM.so
#
#MESSAGE("--DuneDM- DuneDM section begin")
file(GLOB DuneDM_SRCS ${PROJECT_SOURCE_DIR}/*.cpp)
#file(GLOB DuneDM_HDRS ${PROJECT_SOURCE_DIR}/tree/*.h)
#list(REMOVE_ITEM DuneDM_HDRS ${PROJECT_SOURCE_DIR}/tree/LinkDef.h)

#MESSAGE("DuneDM_SRCS=${DuneDM_SRCS}")
#MESSAGE(" ")
#MESSAGE("DuneDM_HDRS=${DuneDM_HDRS}")
#MESSAGE(" ")

#set(DuneDM_LINKDEF ${PROJECT_SOURCE_DIR}/tree/LinkDef.h)
# do not include ".cxx" here!!!!
#set(DuneDM_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/DuneDMDict)

#MESSAGE("--DuneDM- about to ROOT_GENERATE_DICTIONARY ${DuneDM_DICTIONARY}")
#MESSAGE("--rwh- ROOT_GENERATE_DICTIONARY ${DuneDM_DICTIONARY} ${DuneDM_HDRS} ${DuneDM_LINKDEF} ")
#ROOT_GENERATE_DICTIONARY("${DuneDM_DICTIONARY}" "${DuneDM_HDRS}" LINKDEF "${DuneDM_LINKDEF}" )

# add dictionary to list of source files
#SET(DuneDM_SRCS ${DuneDM_SRCS} ${DuneDM_DICTIONARY})

add_executable(DuneDM ${DuneDM_SRCS})
target_link_libraries(DuneDM ${ROOT_LIBRARIES})


