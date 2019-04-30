cmake_minimum_required(VERSION 3.2)
project(epacts VERSION 3.4.2)
include(CMakePackageConfigHelpers)

foreach (PROG_TO_CHECK cat cut dirname gnuplot grep groff gs head make mkdir mv perl pfbtops R rm Rscript sort zcat)
    find_program(${PROG_TO_CHECK}_EXE ${PROG_TO_CHECK})
    if(NOT ${PROG_TO_CHECK}_EXE)
        message(FATAL_ERROR "${PROG_TO_CHECK} not found!")
    endif()
    message(STATUS "Found ${PROG_TO_CHECK}")
endforeach()

add_definitions(-DEPACTS_VERSION="${PROJECT_VERSION}")

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;${CMAKE_FIND_LIBRARY_SUFFIXES}") # Prefer libz.a when both are available

find_package(ZLIB REQUIRED)
find_library(HTS_LIBRARY hts)
find_package(savvy REQUIRED)
#find_package(Eigen3 REQUIRED)
find_package(Threads)

#add_custom_target(epactsR ALL COMMAND R CMD INSTALL -l ${CMAKE_BINARY_DIR} data/epactsR_3.3.2.tar.gz WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMENT "Builing epactsR ...")
#add_custom_target(mmSKAT ALL COMMAND R CMD INSTALL -l ${CMAKE_BINARY_DIR} data/mmSKAT_0.95.tar.gz WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} COMMENT "Builing mmSKAT ...")

add_custom_target(epactsR-build
                  COMMAND R CMD build  ${CMAKE_CURRENT_SOURCE_DIR}/epactsR
                  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                  COMMENT "Builing epactsR ...")
add_custom_target(epactsR-install
                  ALL
                  DEPENDS epactsR-build
                  #${CMAKE_COMMAND} -E env CPLUS_INCLUDE_PATH=${CMAKE_SOURCE_DIR}/cget/include LIBRARY_PATH=${CMAKE_SOURCE_DIR}/cget/lib R CMD INSTALL -l ${CMAKE_BINARY_DIR} epactsR_3.3.3.tar.gz
                  COMMAND R CMD INSTALL -l ${CMAKE_BINARY_DIR} --configure-vars='INCLUDE_DIR=${CGET_PREFIX}/include LINK_DIR=${CGET_PREFIX}/lib' epactsR_${PROJECT_VERSION}.tar.gz
                  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                  COMMENT "Installing epactsR ...")

add_custom_target(mmSKAT-build
                  COMMAND R CMD build ${CMAKE_CURRENT_SOURCE_DIR}/mmSKAT
                  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                  COMMENT "Builing mmSKAT ...")
add_custom_target(mmSKAT-install
                  ALL
                  DEPENDS mmSKAT-build
                  COMMAND R CMD INSTALL -l ${CMAKE_BINARY_DIR} mmSKAT_0.95.tar.gz
                  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                  COMMENT "Installing mmSKAT ...")

add_custom_target(epacts-version
                  ALL
                  COMMAND printf '\#!/usr/bin/perl -wl\\nprint STDOUT \"%s\"\;' ${PROJECT_VERSION} > epacts-version)

add_library(libtabix src/bedidx.c src/knetfile.c src/index.c)
set_target_properties(libtabix PROPERTIES PREFIX "")
target_link_libraries(libtabix ${HTS_LIBRARY} ${ZLIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

add_executable(pEmmax src/pEmmax.cpp src/Parameters.cpp src/Error.cpp src/PhredHelper.cpp src/codonHelper.cpp src/fVcf.h)
target_link_libraries(pEmmax libtabix savvy)
target_include_directories(pEmmax PRIVATE ${CMAKE_SOURCE_DIR}/src)
#target_include_directories(pEmmax PRIVATE ${CMAKE_SOURCE_DIR}/src ${EIGEN3_INCLUDE_DIR})

add_executable(vcfast src/vcfast.cpp src/Error.cpp src/PhredHelper.cpp src/Parameters.cpp src/codonHelper.cpp src/genomePosition.cpp)
target_link_libraries(vcfast libtabix savvy)

add_executable(chaps src/chaps.cpp src/Error.cpp src/PhredHelper.cpp src/Parameters.cpp)
target_link_libraries(chaps libtabix savvy)

add_executable(anno src/Anno.cpp src/IO.cpp)
target_link_libraries(anno libtabix)

add_executable(bgzip src/bgzip.c src/bedidx.c src/knetfile.c src/kstring.c)
target_link_libraries(bgzip ${HTS_LIBRARY} ${ZLIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

add_executable(tabix src/tabix.c src/bedidx.c src/knetfile.c src/kstring.c src/index.c)
target_link_libraries(tabix ${HTS_LIBRARY} ${ZLIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

add_executable(parse-sample-ids src/parse_sample_ids.cpp)
target_link_libraries(parse-sample-ids savvy)

file(GLOB_RECURSE EPACTS_PERL_SCRIPTS "scripts/ep*" "scripts/wGetOptions.pm")
file(GLOB_RECURSE EPACTS_R_SCRIPTS  "data/*.R")
file(GLOB_RECURSE EPACTS_DATA_FILES "data/*")

install(FILES ${EPACTS_PERL_SCRIPTS} ${CMAKE_BINARY_DIR}/epacts-version
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
        DESTINATION bin)
install(FILES ${EPACTS_DATA_FILES} DESTINATION share/EPACTS)
install(TARGETS pEmmax vcfast chaps anno bgzip tabix parse-sample-ids RUNTIME DESTINATION bin)
install(DIRECTORY ${CMAKE_BINARY_DIR}/epactsR DESTINATION lib) #install(DIRECTORY ${CMAKE_BINARY_DIR}/epactsR DESTINATION lib)
install(DIRECTORY ${CMAKE_BINARY_DIR}/mmSKAT DESTINATION lib) #install(DIRECTORY ${CMAKE_BINARY_DIR}/mmSKAT DESTINATION lib)

set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set(CPACK_DEBIAN_PACKAGE_DEPENDS "ghostscript, r-base, groff, gnuplot") 
set(CPACK_RPM_PACKAGE_REQUIRES "ghostscript, R, groff, gnuplot")
include(CPack)
