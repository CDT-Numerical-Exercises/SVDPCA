# Setup project
cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(BaseProject)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# set build type variable in preprocessor
if(CMAKE_BUILD_TYPE)
    string(TOUPPER "${CMAKE_BUILD_TYPE}" _upper_build_type)
    add_compile_definitions(BUILD_${_upper_build_type}=1)
endif()

add_compile_definitions(HAVE_INLINE=1)

# add include dir
include_directories(${PROJECT_SOURCE_DIR}/include)

# add libraries
link_directories(${PROJECT_SOURCE_DIR}/lib)

# add sources and headers
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)

# add executables
add_executable(svd_test svd_test.cpp ${sources} ${headers})
target_compile_definitions(svd_test PRIVATE VERBOSE_SVD=1)
add_executable(pca_test pca_test.cpp ${sources} ${headers})
add_executable(problem1 problem1.cpp ${sources} ${headers})
add_executable(problem2 problem2.cpp ${sources} ${headers})

# copy runtime dependencies
configure_file(svdpca-problem1-data.csv svdpca-problem1-data.csv COPYONLY)
configure_file(Faces.png Faces.png COPYONLY)

# add libraries
find_package(GSL REQUIRED)

# add Boost libraries
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.71.0 COMPONENTS iostreams system filesystem)
include_directories(${Boost_INCLUDE_DIRS})

# connect libraries to executables
target_link_libraries(svd_test GSL::gsl GSL::gslcblas ${Boost_LIBRARIES})
target_link_libraries(pca_test GSL::gsl GSL::gslcblas ${Boost_LIBRARIES})
target_link_libraries(problem1 GSL::gsl GSL::gslcblas ${Boost_LIBRARIES})
target_link_libraries(problem2 GSL::gsl GSL::gslcblas ${Boost_LIBRARIES})
