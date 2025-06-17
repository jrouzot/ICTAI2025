cmake_minimum_required(VERSION 3.25)
project(omdpi VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)

# Add paths for CPLEX
set(MAIN_DIR "/home/jrouzot/")
set(CPLEX_INCLUDE_DIR "${MAIN_DIR}/cplex/include")
set(CONCERT_INCLUDE_DIR "${MAIN_DIR}/concert/include")
set(CPLEX_LIB_DIR "${MAIN_DIR}/cplex/lib/x86-64_linux/static_pic")
set(CONCERT_LIB_DIR "${MAIN_DIR}/concert/lib/x86-64_linux/static_pic")

# Include CPLEX headers
include_directories(${CPLEX_INCLUDE_DIR})
include_directories(${CONCERT_INCLUDE_DIR})

# Link CPLEX libraries
find_library(CPLEX_LIB NAMES cplex PATHS ${CPLEX_LIB_DIR} NO_DEFAULT_PATH)
find_library(ILOCPLEX_LIB NAMES ilocplex PATHS ${CPLEX_LIB_DIR} NO_DEFAULT_PATH)
find_library(CONCERT_LIB NAMES concert PATHS ${CONCERT_LIB_DIR} NO_DEFAULT_PATH)

if (CPLEX_LIB AND ILOCPLEX_LIB AND CONCERT_LIB)
    message(STATUS "CPLEX and CONCERT libraries found.")
    message(STATUS "CPLEX library path: ${CPLEX_LIB}")
    message(STATUS "ILOCPLEX library path: ${ILOCPLEX_LIB}")
    message(STATUS "CONCERT library path: ${CONCERT_LIB}")
else()
    message(FATAL_ERROR "CPLEX or CONCERT library not found.")
endif()


# TO REMOVE FOR PERFORMANCES
#[[
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wpedantic -fsanitize=address -g")
]]

add_executable(${PROJECT_NAME} main.cpp
        src/Solution/Solution.cpp
        src/Instance/Instance.cpp
        src/Simulator/Simulator.cpp
        src/Algorithm/Algorithm.cpp
        src/DownlinkCount/DownlinkCount.cpp
        src/LubySequence/LubySequence.cpp
        src/CplexFlowDataTransferModel/CplexFlowDataTransferModel.cpp)

# Link CPLEX libraries
# THE ORDER OF THE LIBRARIES MATTER !!!
target_link_libraries(${PROJECT_NAME} ${ILOCPLEX_LIB} ${CPLEX_LIB} ${CONCERT_LIB})