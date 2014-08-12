cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Create main project.
project(gotm Fortran)

# Set version information.
set(gotm_VERSION_MAJOR 4)
set(gotm_VERSION_MINOR 1)
set(gotm_VERSION_PATCH 0)
set(gotm_VERSION "${gotm_VERSION_MAJOR}.${gotm_VERSION_MINOR}.${gotm_VERSION_PATCH}")

# Generate include file with version information.
configure_file(../include/version.h.in version.h)

# Specify default build type for single-type systems (not VS)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set (CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Customize compiler flags
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set (CMAKE_Fortran_FLAGS "-ffree-line-length-none")
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray")
  set (CMAKE_Fortran_FLAGS "-N 1023")
endif()

# Use use position-independent code (-fPIC) everywhere if building shared libraries
if(BUILD_SHARED_LIBS)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

# Set default installation prefix.
if(WIN32)
  if(DEFINED ENV{LOCALAPPDATA})
    set(DEFAULT_PREFIX "$ENV{LOCALAPPDATA}/gotm")
  else()
    set(DEFAULT_PREFIX "$ENV{APPDATA}/gotm")
  endif()
else()
  set(DEFAULT_PREFIX "$ENV{HOME}/local/gotm")
endif()
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  # Still on CMake default - override with our own default.
  set(CMAKE_INSTALL_PREFIX ${DEFAULT_PREFIX} CACHE PATH "Directory to install GOTM in" FORCE)
else()
  # Just set the doc string for the variable.
  set(CMAKE_INSTALL_PREFIX ${DEFAULT_PREFIX} CACHE PATH "Directory to install GOTM in")
endif()

# Global preprocessor definitions
add_definitions(-DFORTRAN95)

option(EXTRA_OUTPUT "Include additional turbulence diagnostics in output" OFF)
mark_as_advanced(EXTRA_OUTPUT)
if(EXTRA_OUTPUT)
  add_definitions(-DEXTRA_OUTPUT)
endif()

# Global include directories
include_directories("${PROJECT_SOURCE_DIR}/../include"
                    "${PROJECT_SOURCE_DIR}/util"
                    "${CMAKE_BINARY_DIR}/modules"
                   )
# Libraries added in alphabetic order - dependencies are dealt with below 
add_library(airsea OBJECT
            airsea/airsea.F90
            airsea/airsea_fluxes.F90
            airsea/airsea_variables.F90
            airsea/back_radiation.F90
            airsea/fairall.F90
            airsea/humidity.F90
            airsea/kondo.F90
            airsea/short_wave_radiation.F90
           )
set_property(TARGET airsea APPEND PROPERTY COMPILE_DEFINITIONS GUSTINESS)

#option(USE_ICE "Include ICE support" OFF)
#if(USE_ICE)
#  add_library(ice OBJECT
#             )
#endif(USE_ICE)

add_library(input OBJECT
            input/input.F90
           )

add_library(meanflow OBJECT
            meanflow/buoyancy.F90
            meanflow/convectiveadjustment.F90
            meanflow/coriolis.F90
            meanflow/extpressure.F90
            meanflow/friction.F90
            meanflow/intpressure.F90
            meanflow/meanflow.F90
            meanflow/salinity.F90
            meanflow/shear.F90
            meanflow/stratification.F90
            meanflow/temperature.F90
            meanflow/uequation.F90
            meanflow/updategrid.F90
            meanflow/vequation.F90
           )

add_library(observations OBJECT
            observations/analytical_profile.F90
            observations/const_NNS.F90
            observations/const_NNT.F90
            observations/observations.F90
           )

add_library(output OBJECT
            output/asciiout.F90
            output/ncdfout.F90
            output/output.F90
           )

add_library(turbulence
            turbulence/algebraiclength.F90
            turbulence/alpha_mnb.F90
            turbulence/cmue_a.F90
            turbulence/cmue_b.F90
            turbulence/cmue_c.F90
            turbulence/cmue_d.F90
            turbulence/cmue_ma.F90
            turbulence/cmue_rf.F90
            turbulence/cmue_sg.F90
            turbulence/compute_cpsi3.F90
            turbulence/compute_rist.F90
            turbulence/dissipationeq.F90
            turbulence/epsbalgebraic.F90
            turbulence/fk_craig.F90
            turbulence/genericeq.F90
            turbulence/gotm_lib_version.F90
            turbulence/internal_wave.F90
            turbulence/ispralength.F90
            turbulence/kbalgebraic.F90
            turbulence/kbeq.F90
            turbulence/kpp.F90
            turbulence/lengthscaleeq.F90
            turbulence/potentialml.F90
            turbulence/production.F90
            turbulence/q2over2eq.F90
            turbulence/r_ratio.F90
            turbulence/tkealgebraic.F90
            turbulence/tkeeq.F90
            turbulence/turbulence.F90
            turbulence/algebraiclength.F90
            turbulence/alpha_mnb.F90
            turbulence/cmue_a.F90
            turbulence/cmue_b.F90
            turbulence/cmue_c.F90
            turbulence/cmue_d.F90
            turbulence/cmue_ma.F90
            turbulence/cmue_rf.F90
            turbulence/cmue_sg.F90
            turbulence/compute_cpsi3.F90
            turbulence/compute_rist.F90
            turbulence/dissipationeq.F90
            turbulence/epsbalgebraic.F90
            turbulence/fk_craig.F90
            turbulence/genericeq.F90
            turbulence/gotm_lib_version.F90
            turbulence/internal_wave.F90
            turbulence/ispralength.F90
            turbulence/kbalgebraic.F90
            turbulence/kbeq.F90
            turbulence/kpp.F90
            turbulence/lengthscaleeq.F90
            turbulence/potentialml.F90
            turbulence/production.F90
            turbulence/q2over2eq.F90
            turbulence/r_ratio.F90
            turbulence/tkealgebraic.F90
            turbulence/tkeeq.F90
            turbulence/turbulence.F90
            turbulence/variances.F90
           )

add_library(util
            util/adv_center.F90
            util/convert_fluxes.F90
            util/diff_center.F90
            util/diff_face.F90
            util/eqstate.F90
            util/gridinterpol.F90
            util/lagrange.F90
            util/ode_solvers.F90
            util/time.F90
            util/tridiagonal.F90
            util/util.F90
           )

add_library(gotm OBJECT
            gotm/gotm.F90
           )

option(SEAGRASS "Enable seagrass module" OFF)
if(SEAGRASS)
  add_definitions(-DSEAGRASS)
  add_library(seagrass OBJECT
              extras/seagrass/seagrass.F90
             )
  set(EXTRA_SRCS ${EXTRA_SRCS} $<TARGET_OBJECTS:seagrass>)
  add_dependencies(seagrass output meanflow)
  add_dependencies(gotm seagrass)
endif()

option(BIO "Enable legacy biogeochemical modules" OFF)
if(BIO)
  add_definitions(-DBIO)
  add_library(bio OBJECT
              extras/bio/bio.F90
              extras/bio/bio_fluxes.F90
              extras/bio/bio_save.F90
              extras/bio/bio_alloc_info.F90
              extras/bio/bio_var.F90
              extras/bio/template/bio_template.F90
              extras/bio/npzd/bio_npzd.F90
              extras/bio/cl/bio_cl.F90
              extras/bio/iow/bio_iow.F90
              extras/bio/fasham/bio_fasham.F90
              extras/bio/sed/bio_sed.F90
              #extras/bio/mab/bio_mab.F90
              #extras/bio/rolm/bio_rolm.F90
              extras/bio/npzd_fe/bio_npzd_fe.F90
              extras/bio/photo/bio_photo.F90
              extras/bio/mangan/bio_mangan.F90
             )
  set_property(TARGET bio APPEND PROPERTY COMPILE_DEFINITIONS NO_0D_BIO BIO_TEMPLATE BIO_NPZD BIO_CL BIO_IOW BIO_FASHAM BIO_SED BIO_NPZD_FE BIO_PHOTO BIO_MANGAN)
  set(EXTRA_SRCS ${EXTRA_SRCS} $<TARGET_OBJECTS:bio>)
  add_dependencies(bio util output)
  add_dependencies(gotm bio)
endif()

target_link_libraries(turbulence util)
add_dependencies(input util)
add_dependencies(observations util input)
add_dependencies(airsea util input)
add_dependencies(meanflow util observations)
add_dependencies(output util meanflow turbulence observations airsea)
add_dependencies(gotm input output util meanflow turbulence observations airsea)

set_property(TARGET util turbulence PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

option(USE_NetCDF "Enable output in NetCDF format" ON)
if(USE_NetCDF)
  find_package(NetCDF REQUIRED)
  add_definitions(-DNETCDF_FMT -DREAL_4B=real\(4\))
  include_directories("${NetCDF_INCLUDE_DIRS}")
endif(USE_NetCDF)

# For now set to off as default - should maybe be changed
# Use FABM_PREFIX to locate include and lib directories
option(USE_FABM "Include FABM support" ON)
if(USE_FABM)
  find_package(FABM REQUIRED)
  add_library(gotm_fabm
              fabm/gotm_fabm.F90
              fabm/gotm_fabm_input.F90
              fabm/gotm_fabm_output.F90
             )
  add_dependencies(gotm_fabm output)
  target_link_libraries(gotm_fabm util)
  set_property(TARGET gotm_fabm APPEND PROPERTY INCLUDE_DIRECTORIES "${FABM_INCLUDE_DIRS}")
  set_property(TARGET gotm_fabm PROPERTY Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
  add_definitions(-D_FABM_)

  add_dependencies(gotm gotm_fabm)
  set_property(TARGET gotm APPEND PROPERTY INCLUDE_DIRECTORIES "${FABM_INCLUDE_DIRS}")
endif(USE_FABM)

# Build GOTM executable.
add_executable(gotm_exe
               gotm/main.F90
               $<TARGET_OBJECTS:gotm>
               $<TARGET_OBJECTS:airsea>
               $<TARGET_OBJECTS:meanflow>
               $<TARGET_OBJECTS:observations>
               $<TARGET_OBJECTS:output>
               $<TARGET_OBJECTS:input>
               ${EXTRA_SRCS}
              )
set_property(TARGET gotm_exe PROPERTY OUTPUT_NAME "gotm")
target_link_libraries(gotm_exe turbulence util)

if(USE_NetCDF)
  target_link_libraries(gotm_exe "${NetCDF_LIBRARIES}")
  if (NetCDF_STATIC_MSVC_BUILD)
    set_property(TARGET gotm_exe PROPERTY LINK_FLAGS_DEBUG "/NODEFAULTLIB:\"libcmtd\"")
    target_compile_options(gotm_exe PRIVATE "/libs:static")
  endif()
endif(USE_NetCDF)

if(USE_FABM)
  target_link_libraries(gotm_exe gotm_fabm ${FABM_LIBRARIES})
endif(USE_FABM)

# Install GOTM library and include files, unless the host overruled this.
install(TARGETS gotm_exe DESTINATION bin)
install(TARGETS turbulence util DESTINATION lib)
if(USE_FABM)
  install(TARGETS gotm_fabm DESTINATION lib)
endif()
install(DIRECTORY ${CMAKE_BINARY_DIR}/modules/ DESTINATION include)
