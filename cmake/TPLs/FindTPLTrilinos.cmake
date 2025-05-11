INCLUDE(TribitsTplDeclareLibraries)

# TRIBITS_TPL_DECLARE_LIBRARIES( Trilinos
#   REQUIRED_HEADERS Epetra_Comm.h
#   REQUIRED_LIBS_NAMES "epetra"
#   )

if(Trilinos_LIBRARY_DIRS)
    set (Trilinos_DIR ${Trilinos_DIR} "${Trilinos_LIBRARY_DIRS}/cmake/Trilinos")
endif()
if (Trilinos_INCLUDE_DIRS)
    set (Trilinos_DIR ${Trilinos_DIR} ${Trilinos_INCLUDE_DIRS})
endif()

MESSAGE("Trilinos_DIR: ${Trilinos_DIR}")

# Here I am looking for TrilinosConfig.cmake and I will import it.
find_package (Trilinos NO_MODULE HINTS ${Trilinos_DIR})

# Stop cmake if Trilinos is not found.
if (NOT Trilinos_FOUND)
  message (FATAL_ERROR "Could not find Trilinos!")
endif ()

# Echo trilinos build info just for fun
MESSAGE("\nFound Trilinos!  Here are the details: ")
MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
MESSAGE("   Trilinos_LIBRARY_DIRS = ${Trilinos_LIBRARY_DIRS}")
MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
MESSAGE("   Trilinos_TPL_LIBRARY_DIRS = ${Trilinos_TPL_LIBRARY_DIRS}")
MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
MESSAGE("End of Trilinos details\n")

# Attempt to determine Trilinos library directories if not set by TrilinosConfig.cmake
set(FEDD_Trilinos_LIBRARY_DIRS "")
if(DEFINED Trilinos_LIBRARY_DIRS AND NOT "${Trilinos_LIBRARY_DIRS}" STREQUAL "")
  set(FEDD_Trilinos_LIBRARY_DIRS ${Trilinos_LIBRARY_DIRS})
  message(STATUS "FindTPLTrilinos.cmake: Using Trilinos_LIBRARY_DIRS from TrilinosConfig.cmake: ${FEDD_Trilinos_LIBRARY_DIRS}")
else()
  if(DEFINED Trilinos_INSTALL_DIR AND NOT "${Trilinos_INSTALL_DIR}" STREQUAL "")
    set(potential_lib_dir "${Trilinos_INSTALL_DIR}/lib")
    if(EXISTS "${potential_lib_dir}")
      set(FEDD_Trilinos_LIBRARY_DIRS "${potential_lib_dir}")
      message(STATUS "FindTPLTrilinos.cmake: Trilinos_LIBRARY_DIRS not set by TrilinosConfig.cmake. Inferred as: ${FEDD_Trilinos_LIBRARY_DIRS}")
    else()
      message(WARNING "FindTPLTrilinos.cmake: Trilinos_LIBRARY_DIRS not set and inferred path ${potential_lib_dir} does not exist. Library paths might be missing.")
    endif()
  else()
    message(WARNING "FindTPLTrilinos.cmake: Trilinos_LIBRARY_DIRS not set and Trilinos_INSTALL_DIR not available to infer. Library paths might be missing.")
  endif()
endif()


# if("${Trilinos_VERSION_MAJOR}" GREATER 10)
#   set (HAVE_TRILINOS_GT_10_6 TRUE)
#   message (STATUS "Using Trilinos > 10.6 : " ${Trilinos_VERSION_MAJOR} "." ${Trilinos_VERSION_MINOR})
# else()
#  if("${Trilinos_VERSION_MAJOR}" SMALLER 10)
#    message (STATUS "Using Trilinos <= 10.6 : " ${Trilinos_VERSION_MAJOR} "." ${Trilinos_VERSION_MINOR})
#  else()
#    if("${Trilinos_VERSION_MINOR}" GREATER 6)
#     set (HAVE_TRILINOS_GT_10_6 TRUE)
#     message (STATUS "Using Trilinos > 10.6 : " ${Trilinos_VERSION_MAJOR} "." ${Trilinos_VERSION_MINOR})
#    else()
#     message (STATUS "Using Trilinos <= 10.6 : " ${Trilinos_VERSION_MAJOR} "." ${Trilinos_VERSION_MINOR})
#    endif()
#  endif()
# endif ()

# Here it will be better just to raise a warning or have if(USE_TRILINOS_COMPILERS)
# Make sure to use same compilers and flags as Trilinos
IF(NOT ${CMAKE_CXX_COMPILER} STREQUAL ${Trilinos_CXX_COMPILER})
  MESSAGE(STATUS "the selected compiler differs from Trilinos CXX compiler")
  MESSAGE(STATUS "CMAKE_CXX_COMPILER:    " ${CMAKE_CXX_COMPILER})
  MESSAGE(STATUS "Trilinos_CXX_COMPILER: " ${Trilinos_CXX_COMPILER})
ENDIF()
IF(NOT ${CMAKE_C_COMPILER} STREQUAL ${Trilinos_C_COMPILER})
  MESSAGE(STATUS "the selected compiler differs from Trilinos C compiler")
  MESSAGE(STATUS "CMAKE_C_COMPILER:    " ${CMAKE_C_COMPILER})
  MESSAGE(STATUS "Trilinos_C_COMPILER: " ${Trilinos_C_COMPILER})
ENDIF()
IF(NOT ${CMAKE_Fortran_COMPILER} STREQUAL ${Trilinos_Fortran_COMPILER})
  MESSAGE(STATUS "the selected compiler differs from Trilinos Fortran compiler")
  MESSAGE(STATUS "CMAKE_Fortran_COMPILER:    " ${CMAKE_C_COMPILER})
  MESSAGE(STATUS "Trilinos_Fortran_COMPILER: " ${Trilinos_C_COMPILER})
ENDIF()

# Optional Packages (to be moved outside with COMPONENTS ...)
list (APPEND XLib_OPTIONAL_Trilinos_PKGS
  "NOX" "Thyra" "Rythmos" "Teko" "Stratimikos" "Isorropia" "ShyLU" "Zoltan2" "MueLu")

# Required packages (to be moved outside, like REQUIRED COMPONENTS ...)
list (APPEND XLib_REQUIRED_Trilinos_PKGS
  "Belos" "Epetra" "EpetraExt" "ShyLU_DDFROSch" "Stratimikos" "Teko" "Teuchos" "Thyra" "Tpetra" "Xpetra")

# Start scanning Trilinos configuration
foreach (TYPE IN ITEMS "OPTIONAL" "REQUIRED")
  foreach (PKG IN LISTS XLib_${TYPE}_Trilinos_PKGS)
    # Look for PKG
    list (FIND Trilinos_PACKAGE_LIST "${PKG}" PKG_FOUND)
    if (PKG_FOUND GREATER -1)
      # Found! Let's announce it!
      message (STATUS "Trilinos :: ${PKG} Found!")
      list (APPEND XLib_Trilinos_LIBRARIES "${${PKG}_LIBRARIES}")
      list (APPEND XLib_Trilinos_TPL_INCLUDE_DIRS "${${PKG}_TPL_INCLUDE_DIRS}")
      list (APPEND XLib_Trilinos_TPL_LIST "${${PKG}_TPL_LIST}")
      string (TOUPPER ${PKG} UPKG)
      set (${UPKG}_FOUND True)
      set (HAVE_TRILINOS_${UPKG} True)
    else ()
      if (TYPE STREQUAL "REQUIRED")
        message (FATAL_ERROR "Trilinos :: ${PKG} NOT Found!")
      else ()
        message (WARNING "Trilinos :: ${PKG} NOT Found! Some test might not compile properly ...")
      endif ()
    endif ()
  endforeach (PKG)
endforeach (TYPE)

# Cleaning duplicates
list (REVERSE XLib_Trilinos_TPL_LIST)
list (REMOVE_DUPLICATES XLib_Trilinos_TPL_LIST)
list (REVERSE XLib_Trilinos_TPL_LIST)
list (REVERSE XLib_Trilinos_LIBRARIES)
list (REMOVE_DUPLICATES XLib_Trilinos_LIBRARIES)
list (REVERSE XLib_Trilinos_LIBRARIES)
list (REVERSE Trilinos_TPL_LIBRARIES)
list (REMOVE_DUPLICATES Trilinos_TPL_LIBRARIES)
list (REVERSE Trilinos_TPL_LIBRARIES)
set (XLib_Trilinos_TPL_LIBRARIES ${Trilinos_TPL_LIBRARIES})

list (REMOVE_DUPLICATES XLib_Trilinos_TPL_INCLUDE_DIRS)

list (APPEND XLib_Trilinos_INCLUDE_DIRS
  ${Trilinos_INCLUDE_DIRS}
  ${XLib_Trilinos_TPL_INCLUDE_DIRS})

# Construct XLib_Trilinos_LIBS carefully
set(XLib_Trilinos_LIBS "")
if(NOT "${FEDD_Trilinos_LIBRARY_DIRS}" STREQUAL "")
  foreach(lib_dir IN LISTS FEDD_Trilinos_LIBRARY_DIRS) # Handle if it's a list
    list(APPEND XLib_Trilinos_LIBS "-L${lib_dir}")
  endforeach()
endif()

foreach (LIB IN LISTS XLib_Trilinos_LIBRARIES) # XLib_Trilinos_LIBRARIES are package library names
  list(APPEND XLib_Trilinos_LIBS "-l${LIB}")
endforeach (LIB)
# XLib_Trilinos_TPL_LIBRARIES is set from Trilinos_TPL_LIBRARIES earlier.
# These are assumed to be linker arguments already (e.g. from TPLs *of* Trilinos)
if(XLib_Trilinos_TPL_LIBRARIES)
  list(APPEND XLib_Trilinos_LIBS ${XLib_Trilinos_TPL_LIBRARIES})
endif()

# TPLs
foreach (TPL IN ITEMS "ParMETIS" "Boost" "LAPACK" "BLAS" "UMFPACK" "SuperLU" "SuperLUDist" "HDF5")
    list (FIND XLib_Trilinos_TPL_LIST ${TPL} TPL_FOUND)
  if (TPL_FOUND GREATER -1)
    string (TOUPPER ${TPL} UTPL)
    set (${UTPL}_IS_IN_TRILINOS True)
  endif()
endforeach (TPL)

# Filling variables needed by the TriBITS system
set (TPL_Trilinos_INCLUDE_DIRS ${XLib_Trilinos_INCLUDE_DIRS})
set (TPL_Trilinos_LIBRARY_DIRS ${FEDD_Trilinos_LIBRARY_DIRS}) # Use the determined path; will be "" if not found, satisfying AssertDefined
set (TPL_Trilinos_LIBRARIES Trilinos::all_selected_libs)



