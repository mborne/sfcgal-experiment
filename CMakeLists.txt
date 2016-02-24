cmake_minimum_required( VERSION 2.8 )
project( SFCGALExperiment )

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_SOURCE_DIR}/cmake/Modules")

ADD_DEFINITIONS(-DBOOST_LOG_DYN_LINK=1)
SET(Boost_USE_STATIC_LIBS OFF)
find_package(Boost 1.54 COMPONENTS log REQUIRED)
find_package (Threads)

find_package( CGAL 4.3 COMPONENTS Core REQUIRED )
find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)

add_executable( morphing-distance
    morphing-distance.cpp
)

target_link_libraries(morphing-distance
    ${CGAL_LIBRARIES}
    ${Boost_LIBRARIES}
    ${MPFR_LIBRARIES}
    ${GMP_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
)
