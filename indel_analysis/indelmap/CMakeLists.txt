##########################################################################
# CMakeLists.txt
#########################################################################*/
cmake_minimum_required(VERSION 2.6)

execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

project (indelmap)
set (CMAKE_CXX_STANDARD 11)

set (INDELMAP_OUTPUT_DIR ${CMAKE_SOURCE_DIR}/bin CACHE STRING "Where do you want the output binaries to install to?" )

set ( SRC_FILES indel.cpp indel.h oligo.cpp oligo.h readmap.cpp readmap.h gen.cpp gen.h)
set ( TEST_SRC_FILES test.cpp test.h )

add_executable( indelmap ${SRC_FILES} indelmap.cpp)
add_executable( indelmaptest ${SRC_FILES} ${TEST_SRC_FILES} )
add_executable( indelgen ${SRC_FILES} indelgen.cpp )
add_executable( indelgen_i1 ${SRC_FILES} indelgen_i1.cpp )
add_executable( indelmh ${SRC_FILES} indelmh.cpp )
add_definitions("-DGIT_COMMIT_HASH=${GIT_COMMIT_HASH}")

configure_file(
  ${CMAKE_SOURCE_DIR}/version.h.in
  ${CMAKE_BINARY_DIR}/generated/version.h
)
include_directories(${CMAKE_BINARY_DIR}/generated)

install( TARGETS indelmap 
		 DESTINATION ${INDELMAP_OUTPUT_DIR} )
install( TARGETS indelgen 
		 DESTINATION ${INDELMAP_OUTPUT_DIR} )
install( TARGETS indelgen_i1
		 DESTINATION ${INDELMAP_OUTPUT_DIR} )
install( TARGETS indelmh 
		 DESTINATION ${INDELMAP_OUTPUT_DIR} )
