################################################################################
# Custom settings for compiler flags and similar
################################################################################

if (UNIX)

  set ( ADDITIONAL_CXX_FLAGS )
  set ( ADDITIONAL_CXX_DEBUG_FLAGS )
  set ( ADDITIONAL_CXX_RELEASE_FLAGS )
  set ( ADDITIONAL_CXX_RELWITHDEBINFO_FLAGS )
  
  set ( ADDITIONAL_C_FLAGS )
  set ( ADDITIONAL_C_DEBUG_FLAGS )
  set ( ADDITIONAL_C_RELEASE_FLAGS )
  set ( ADDITIONAL_C_RELWITHDEBINFO_FLAGS )

  ################################################################################
  # Defaults
  ################################################################################

  # add our standard flags for Template inclusion
  list(APPEND ADDITIONAL_CXX_FLAGS          "-DINCLUDE_TEMPLATES" )
  list(APPEND ADDITIONAL_C_FLAGS            "-DINCLUDE_TEMPLATES" )
  
  # Increase the template depth as this might be exceeded from time to time
  IF( NOT CMAKE_SYSTEM MATCHES "SunOS*")
    list(APPEND ADDITIONAL_CXX_FLAGS          "-ftemplate-depth-100" )
  ENDIF()


  
  ################################################################################
  # Build/Release Defines
  ################################################################################
  IF( NOT CMAKE_SYSTEM MATCHES "SunOS*")
    list(APPEND ADDITIONAL_CXX_DEBUG_FLAGS          "-DDEBUG" )
    list(APPEND ADDITIONAL_CXX_RELEASE_FLAGS        "-DNDEBUG" )
    list(APPEND ADDITIONAL_CXX_RELWITHDEBINFO_FLAGS "-DDEBUG" )    
    
    list(APPEND ADDITIONAL_C_DEBUG_FLAGS            "-DDEBUG" )
    list(APPEND ADDITIONAL_C_RELEASE_FLAGS          "-DNDEBUG" )
    list(APPEND ADDITIONAL_C_RELWITHDEBINFO_FLAGS   "-DDEBUG" )
  ENDIF()  
  
  ################################################################################
  # Warnings
  ################################################################################
  
  if ("${CMAKE_CXX_COMPILER_ID}" MATCHES ".*Clang.*")
      list(APPEND ADDITIONAL_CXX_FLAGS "-Weverything")
      list(APPEND ADDITIONAL_CXX_FLAGS "-Wno-c++98-compat")
      list(APPEND ADDITIONAL_CXX_FLAGS "-Wno-padded")
      list(APPEND ADDITIONAL_CXX_FLAGS "-Wno-old-style-cast")
      list(APPEND ADDITIONAL_CXX_FLAGS "-Wno-documentation-unknown-command")
      list(APPEND ADDITIONAL_CXX_FLAGS "-Wno-unreachable-code-return")
      # enable later:
      list(APPEND ADDITIONAL_CXX_FLAGS "-Wno-sign-conversion")
      list(APPEND ADDITIONAL_CXX_FLAGS "-Wno-deprecated")
      list(APPEND ADDITIONAL_CXX_FLAGS "-Wno-weak-vtables")
  endif()
  if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
      list(APPEND ADDITIONAL_CXX_FLAGS "-Wall")
      list(APPEND ADDITIONAL_CXX_FLAGS "-Wextra")
      list(APPEND ADDITIONAL_CXX_FLAGS "-Wpedantic")
      list(APPEND ADDITIONAL_CXX_FLAGS "-Wshadow")
      list(APPEND ADDITIONAL_CXX_FLAGS "-Wpointer-arith")
      list(APPEND ADDITIONAL_CXX_FLAGS "-Wcast-qual")
      list(APPEND ADDITIONAL_CXX_FLAGS "-Wconversion")
      list(APPEND ADDITIONAL_CXX_FLAGS "-Wtype-limits")
      list(APPEND ADDITIONAL_CXX_FLAGS "-Wsign-compare")
  endif()
  
  ################################################################################
  # STL Vector checks
  ################################################################################
  
  # Pre initialize stl vector check variable
  if ( NOT STL_VECTOR_CHECKS )
    set ( STL_VECTOR_CHECKS false CACHE BOOL "Include full stl vector checks in debug mode (This option is only used in debug Mode!)" )
  endif ( NOT STL_VECTOR_CHECKS )
  
  # Add a flag to check stl vectors in debugging mode
  if ( STL_VECTOR_CHECKS AND NOT CMAKE_SYSTEM MATCHES "SunOS*"  )
    list(APPEND ADDITIONAL_CXX_DEBUG_FLAGS          "-D_GLIBCXX_DEBUG" )
    list(APPEND ADDITIONAL_CXX_DEBUG_FLAGS          "-D_GLIBCXX_DEBUG_PEDANTIC")
    list(APPEND ADDITIONAL_CXX_RELWITHDEBINFO_FLAGS "-D_GLIBCXX_DEBUG" )
    list(APPEND ADDITIONAL_CXX_RELWITHDEBINFO_FLAGS "-D_GLIBCXX_DEBUG_PEDANTIC")
    
    list(APPEND ADDITIONAL_C_DEBUG_FLAGS            "-D_GLIBCXX_DEBUG" )
    list(APPEND ADDITIONAL_C_DEBUG_FLAGS            "-D_GLIBCXX_DEBUG_PEDANTIC")
    list(APPEND ADDITIONAL_C_RELWITHDEBINFO_FLAGS   "-D_GLIBCXX_DEBUG" )
    list(APPEND ADDITIONAL_C_RELWITHDEBINFO_FLAGS   "-D_GLIBCXX_DEBUG_PEDANTIC")
  endif()

  ################################################################################
  # Process the additional flags:
  ################################################################################

  # Add the debug flags
  foreach( flag ${ADDITIONAL_CXX_FLAGS} ${ADDITIONAL_CXX_DEBUG_FLAGS} )
    list (FIND ${CMAKE_CXX_FLAGS_DEBUG} ${flag} _index)
    if (${_index} EQUAL -1)
      set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${flag} ")
    endif()
  endforeach()

  # Add the release flags
  foreach( flag ${ADDITIONAL_CXX_FLAGS} ${ADDITIONAL_CXX_RELEASE_FLAGS} )
    list (FIND ${CMAKE_CXX_FLAGS_RELEASE} ${flag} _index)
    if (${_index} EQUAL -1)
      set( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${flag} ")
    endif()
  endforeach()

  # Add the release with debug info flags
  foreach( flag ${ADDITIONAL_CXX_FLAGS} ${ADDITIONAL_CXX_RELWITHDEBINFO_FLAGS} )
    list (FIND ${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${flag} _index)
    if (${_index} EQUAL -1)
      set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${flag} ")
    endif()
  endforeach()

  # Add the debug flags
  foreach( flag ${ADDITIONAL_C_FLAGS} ${ADDITIONAL_C_DEBUG_FLAGS} )
    list (FIND ${CMAKE_C_FLAGS_DEBUG} ${flag} _index)
    if (${_index} EQUAL -1)
      set( CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${flag} ")
    endif()
  endforeach()

  # Add the release flags
  foreach( flag ${ADDITIONAL_C_FLAGS} ${ADDITIONAL_C_RELEASE_FLAGS} )
      list (FIND ${CMAKE_C_FLAGS_RELEASE} ${flag} _index)
    if (${_index} EQUAL -1)
      set( CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${flag} ")
    endif()
  endforeach()

  # Add the release with debug info flags
  foreach( flag ${ADDITIONAL_C_FLAGS} ${ADDITIONAL_C_RELWITHDEBINFO_FLAGS} )
    list (FIND ${CMAKE_C_FLAGS_RELWITHDEBINFO} ${flag} _index)
    if (${_index} EQUAL -1)
      set( CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} ${flag} ")
    endif()
  endforeach()

  #TODO : Test and remove it?!
  IF( CMAKE_SYSTEM MATCHES "SunOS*")
    set (CMAKE_CFLAGS_RELEASE "-xO3")
    set (CMAKE_CXX_FLAGS_RELEASE "-xO3")        
  endif ( CMAKE_SYSTEM MATCHES "SunOS*" ) 

  ################################################################################
  # C++ 11 support
  ################################################################################

  # On apple, if we have c++ 11 support, we enable it automatically here
#  if (APPLE)
#
#    include(CheckCXXCompilerFlag)
#    CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
#    CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
#
#    if(COMPILER_SUPPORTS_CXX11)
#       set( CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11")
#       set( CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11")
#       set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11")
#    elseif(COMPILER_SUPPORTS_CXX0X)
#       set( CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG} -std=c++0x")
#       set( CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_DEBUG} -std=c++0x")
#       set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_DEBUG} -std=c++0x")
#    else()
#        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Falling back to non C++11 mode. If you encounter errors, please use a different C++ compiler.")
#    endif()
#
#  endif()

  
endif ()
