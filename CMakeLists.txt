cmake_minimum_required(VERSION 3.14)
project(pileup_cpp VERSION 0.1.0)

# Use C++23
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add compiler warnings
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic")
elseif(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
endif()

# Find htslib (this assumes htslib is installed in the system)
# If htslib is not found, you may need to specify its location with -DHTSLIB_ROOT=path/to/htslib
find_path(HTSLIB_INCLUDE_DIR htslib/sam.h
    HINTS
    $ENV{HTSLIB_ROOT}/include
    ${HTSLIB_ROOT}/include
    /usr/include
    /usr/local/include
)

find_library(HTSLIB_LIBRARY
    NAMES hts libhts
    HINTS
    $ENV{HTSLIB_ROOT}/lib
    ${HTSLIB_ROOT}/lib
    /usr/lib
    /usr/local/lib
)

if(NOT HTSLIB_INCLUDE_DIR OR NOT HTSLIB_LIBRARY)
    message(FATAL_ERROR "Could not find htslib. Please install htslib or set HTSLIB_ROOT.")
endif()

message(STATUS "Found htslib include dir: ${HTSLIB_INCLUDE_DIR}")
message(STATUS "Found htslib library: ${HTSLIB_LIBRARY}")

# Add include directories
include_directories(${HTSLIB_INCLUDE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})  # Add project root for argparse include

# Define source files
set(SOURCES
    main.cpp
    bam_reader.cpp
    bed_reader.cpp
    allele_counter.cpp
    result_writer.cpp
)

# Define header files
set(HEADERS
    common.h
    bam_reader.h
    bed_reader.h
    allele_counter.h
    result_writer.h
    argparse/argparse.hpp  # Add argparse header
)

# Create the executable
add_executable(pileup_cpp ${SOURCES} ${HEADERS})

# Link against htslib
target_link_libraries(pileup_cpp ${HTSLIB_LIBRARY})

# Install targets
install(TARGETS pileup_cpp DESTINATION bin)
