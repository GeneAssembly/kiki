cmake_minimum_required(VERSION 2.6)
project(ki)

set(IS_CRAY $ENV{CRAY_MPICH2_DIR})

if (IS_CRAY)   
  set(CMAKE_C_COMPILER "cc")

elseif (${MPE})			  # long-term:  cmake .. -DMPE=1   
  set(CMAKE_C_COMPILER "mpecc")
elseif ($ENV{MPE})		  # short-term: MPE=1 cmake .. 
  set(CMAKE_C_COMPILER "mpecc")

else ()
  set(CMAKE_C_COMPILER "mpicc")

endif ()


if (${MPE})
  set(CMAKE_C_FLAGS "-mpilog -O3")
elseif ($ENV{MPE})
  set(CMAKE_C_FLAGS "-mpilog -O3")

elseif (${DEBUG}) 
  set(CMAKE_C_FLAGS "-Wall -g -O0")
elseif ($ENV{DEBUG}) 
  set(CMAKE_C_FLAGS "-Wall -g -O0")

elseif (${NDEBUG}) 
  set(CMAKE_C_FLAGS "-O3 -DNDEBUG")
elseif ($ENV{NDEBUG}) 
  set(CMAKE_C_FLAGS "-O3 -DNDEBUG")

else ()
  set(CMAKE_C_FLAGS "-O3")

endif ()

if (IS_CRAY)
  set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")      # avoid -rdynamic error on Cray
elseif (${DEBUG}) 
elseif ($ENV{DEBUG}) 
else ()
  # set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -finline-functions -funroll-loops")
endif ()

if (${SIMPLE_HASH})
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DSIMPLE_HASH")
elseif ($ENV{SIMPLE_HASH})
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DSIMPLE_HASH")
endif ()

# if (APPLE)
#   set(CMAKE_C_LINK_FLAGS "-lm")
# elseif (${CMAKE_SYSTEM} MATCHES "ppc")
#   set(CMAKE_C_LINK_FLAGS "-lm")
# else ()
#   set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -lm")
# endif ()


message ("Compiling flags: ${CMAKE_C_COMPILER} ${CMAKE_C_FLAGS}")
message ("Linking flags:   ${CMAKE_C_LINK_FLAGS}")

set(CORE_SRC
  ./extern.h
  ./debug.h
  ./debug.c
  ./mem.h
  ./mem.c
  ./io.h
  ./io.c
  ./comm.h
  ./comm.c
  ./seq.h
  ./seq.c
  ./hash.h
  ./hash_chain.h
  ./hash_chain.c
  ./hash_oa.h
  ./hash_oa.c
  ./util.h
  ./util.c
  ./ki.c
  ./ki.h
  ./raiphy.c
  ./raiphy.h
  )

# message ("${CORE_SRC}")

add_executable(ki     ${CORE_SRC} ./asm.c)
add_executable(demo   ${CORE_SRC} ./demo.c)
add_executable(search ${CORE_SRC} ./search.c)
add_executable(rai    ${CORE_SRC} ./rai.c)
add_executable(rait   ${CORE_SRC} ./rait.c)


#if (${CMAKE_SYSTEM} MATCHES "ppc")
#  find_library(M_LIB mass)
#else ()
  find_library(M_LIB m)
#endif ()

target_link_libraries(ki     ${M_LIB})
target_link_libraries(demo   ${M_LIB})
target_link_libraries(search ${M_LIB})
target_link_libraries(rai    ${M_LIB})
target_link_libraries(rait   ${M_LIB})
