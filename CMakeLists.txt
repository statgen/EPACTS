cmake_minimum_required(VERSION 3.2)
project(epacts)
add_definitions(-DEPACTS_VERSION=1.0.0)

set(CMAKE_CXX_STANDARD 14)

find_package(ZLIB REQUIRED)
find_library(HTS_LIBRARY hts)
find_package(savvy REQUIRED)

add_library(libtabix src/bedidx.c src/knetfile.c src/index.c)
set_target_properties(libtabix PROPERTIES PREFIX "")
target_link_libraries(libtabix ${ZLIB_LIBRARIES} ${HTS_LIBRARY})

add_executable(pEmmax src/pEmmax.cpp src/Parameters.cpp src/Error.cpp src/PhredHelper.cpp src/codonHelper.cpp src/fVcf.h)
target_link_libraries(pEmmax libtabix savvy)
target_include_directories(pEmmax PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_executable(vcfast src/vcfast.cpp src/Error.cpp src/PhredHelper.cpp src/Parameters.cpp src/codonHelper.cpp src/genomePosition.cpp)
target_link_libraries(vcfast libtabix savvy)

add_executable(chaps src/chaps.cpp src/Error.cpp src/PhredHelper.cpp src/Parameters.cpp)
target_link_libraries(chaps libtabix savvy)

add_executable(anno src/Anno.cpp src/IO.cpp)
target_link_libraries(anno libtabix)

add_executable(bgzip src/bgzip.c src/bedidx.c src/knetfile.c src/kstring.c)
target_link_libraries(bgzip ${ZLIB_LIBRARIES} ${HTS_LIBRARY})

add_executable(tabix src/tabix.c src/bedidx.c src/knetfile.c src/kstring.c src/index.c)
target_link_libraries(tabix ${ZLIB_LIBRARIES} ${HTS_LIBRARY})

install(TARGETS pEmmax vcfast chaps anno bgzip tabix RUNTIME DESTINATION bin)