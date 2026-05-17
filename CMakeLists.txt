cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

# CMP0074: find_package() uses <PackageName>_ROOT variables
cmake_policy(SET CMP0074 NEW)

PROJECT(epic
  LANGUAGES CXX
  DESCRIPTION "DD4hep Geometry Description of the EPIC Experiment"
  )

# Set version based on git
include(cmake/git_version.cmake)
set_git_version(CMAKE_PROJECT_VERSION)

# C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD 17 CACHE STRING "Set the C++ standard to be used")
if(NOT CMAKE_CXX_STANDARD MATCHES "17|20|23")
  message(FATAL_ERROR "Unsupported C++ standard: ${CMAKE_CXX_STANDARD}")
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Build type
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type configuration" FORCE)
  message(STATUS "Setting default build type: ${CMAKE_BUILD_TYPE}")
endif()

# Error on all warnings
if(NOT CMAKE_BUILD_TYPE STREQUAL "Release")
  add_compile_options(-Wall -Wextra -Werror -pedantic)
endif()

# Export compile commands as json for run-clang-tidy
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Install to the top directory by default
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/prefix CACHE PATH "Install in $src/prefix/ directory by default" FORCE)
endif()

# Install into GNU standard directories
include(GNUInstallDirs)

# Dependencies
find_package(DD4hep 1.27 REQUIRED COMPONENTS DDCore DDRec)
find_package(fmt REQUIRED)
find_package(IRT2 2.1.1)

#-----------------------------------------------------------------------------------

# Address sanitizer
option(USE_ASAN "Compile with address sanitizer" OFF)
if(${USE_ASAN})
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    execute_process(
      COMMAND ${CMAKE_CXX_COMPILER} -print-file-name=libclang_rt.asan-${CMAKE_SYSTEM_PROCESSOR}.so
      OUTPUT_VARIABLE PRELOAD_SANITIZER_LIB_ASAN
      OUTPUT_STRIP_TRAILING_WHITESPACE)
  else()
    execute_process(
      COMMAND ${CMAKE_CXX_COMPILER} -print-file-name=libasan.so
      OUTPUT_VARIABLE PRELOAD_SANITIZER_LIB_ASAN
      OUTPUT_STRIP_TRAILING_WHITESPACE)
  endif()
  get_filename_component(ASAN_RUNTIME_DIR "${PRELOAD_SANITIZER_LIB_ASAN}" DIRECTORY)
  list(APPEND CMAKE_BUILD_RPATH "${ASAN_RUNTIME_DIR}")
  list(PREPEND PRELOAD_SANITIZER_LIB ${PRELOAD_SANITIZER_LIB_ASAN})
endif()

# Undefined behavior sanitizer
option(USE_UBSAN "Compile with undefined behavior sanitizer" OFF)
if(${USE_UBSAN})
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    execute_process(
      COMMAND ${CMAKE_CXX_COMPILER} -print-file-name=libclang_rt.ubsan_standalone-${CMAKE_SYSTEM_PROCESSOR}.so
      OUTPUT_VARIABLE PRELOAD_SANITIZER_LIB_UBSAN
      OUTPUT_STRIP_TRAILING_WHITESPACE)
  else()
    execute_process(
      COMMAND ${CMAKE_CXX_COMPILER} -print-file-name=libubsan.so
      OUTPUT_VARIABLE PRELOAD_SANITIZER_LIB_UBSAN
      OUTPUT_STRIP_TRAILING_WHITESPACE)
  endif()
  get_filename_component(UBSAN_RUNTIME_DIR "${PRELOAD_SANITIZER_LIB_UBSAN}" DIRECTORY)
  list(APPEND CMAKE_BUILD_RPATH "${UBSAN_RUNTIME_DIR}")
  # With Clang, the ASAN runtime already integrates UBSAN; preloading
  # libclang_rt.ubsan_standalone alongside libclang_rt.asan causes runtime
  # conflicts (overlapping interceptors). Skip the UBSAN preload when ASAN is
  # also active — ASAN handles both.
  if(NOT (USE_ASAN AND CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
    list(APPEND PRELOAD_SANITIZER_LIB ${PRELOAD_SANITIZER_LIB_UBSAN})
  endif()
endif()

#-----------------------------------------------------------------------------------

function(dd4hep_generate_rootmap library)

  if(APPLE)
    set(ENV_VAR DYLD_LIBRARY_PATH)
  else()
    set(ENV_VAR LD_LIBRARY_PATH)
  endif()
  set(${ENV_VAR}_VALUE $<TARGET_FILE_DIR:${library}>:$<TARGET_FILE_DIR:DD4hep::DD4hepGaudiPluginMgr>:$ENV{${ENV_VAR}})

  set(rootmapfile ${CMAKE_SHARED_MODULE_PREFIX}${library}.components)

  # Build the listcomponents invocation. When sanitizers are active we must
  # LD_PRELOAD the runtime so ASAN initialises before any allocation occurs;
  # loading it lazily via dlopen (RPATH) after DD4hep has already allocated
  # memory causes a crash. To prevent LeakSanitizer from scanning the entire
  # DD4hep/ROOT heap at process exit (which hangs the build), we force
  # detect_leaks=0 for this helper invocation only.
  string(JOIN ":" LD_PRELOAD_VALUE ${PRELOAD_SANITIZER_LIB})
  if(LD_PRELOAD_VALUE)
    set(EXTRA_ENV
      "LD_PRELOAD=${LD_PRELOAD_VALUE}"
      "ASAN_OPTIONS=detect_leaks=0:verify_asan_link_order=0"
      "LSAN_OPTIONS=detect_leaks=0"
    )
  endif()

  add_custom_command(OUTPUT ${rootmapfile}
                     DEPENDS ${library}
                     COMMAND ${CMAKE_COMMAND} -E env
                       "${ENV_VAR}=${${ENV_VAR}_VALUE}"
                       ${EXTRA_ENV}
                       $<TARGET_FILE:DD4hep::listcomponents> -o ${rootmapfile} $<TARGET_FILE_NAME:${library}>
                     WORKING_DIRECTORY ${LIBRARY_OUTPUT_PATH}
                     COMMAND_EXPAND_LISTS
                     )

  add_custom_target(Components_${library} ALL DEPENDS ${rootmapfile})
  SET( install_destination "lib" )
  if( CMAKE_INSTALL_LIBDIR )
    SET( install_destination ${CMAKE_INSTALL_LIBDIR} )
  endif()
  install(FILES $<TARGET_FILE_DIR:${library}>/${rootmapfile}
    DESTINATION ${install_destination}
  )
endfunction()

#-----------------------------------------------------------------------------------

set(a_lib_name ${PROJECT_NAME})

set(DD4HEP_SET_RPATH TRUE)
dd4hep_configure_output(INSTALL ${CMAKE_INSTALL_PREFIX})
dd4hep_set_compiler_flags()

file(GLOB sources CONFIGURE_DEPENDS src/*.cpp)

dd4hep_add_plugin(${a_lib_name}
  SOURCES ${sources}
  USES ROOT::Core ROOT::Gdml
  )
target_link_libraries(${a_lib_name}
  PUBLIC DD4hep::DDCore DD4hep::DDRec fmt::fmt
  )
set(WITH_IRT1_SUPPORT TRUE CACHE STRING "Enable constant definitions used by EICrecon integration in IRT1")
if(WITH_IRT1_SUPPORT)
    target_compile_definitions(${a_lib_name} PRIVATE WITH_IRT1_SUPPORT)
endif()
if(IRT2_FOUND)
    target_compile_definitions(${a_lib_name} PRIVATE WITH_IRT2_SUPPORT)
    target_link_libraries(${a_lib_name} PRIVATE IRT2)
endif()

#-----------------------------------------------------------------------------------

# Address sanitizer
if(${USE_ASAN})
  target_compile_options(${a_lib_name} PRIVATE -fsanitize=address -fno-omit-frame-pointer -g -O1)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_link_options(${a_lib_name} PRIVATE -fsanitize=address -shared-libasan)
  else()
    target_link_options(${a_lib_name} PRIVATE -fsanitize=address)
    target_link_libraries(${a_lib_name} PRIVATE asan)
  endif()
  install(FILES .github/asan.supp .github/lsan.supp
          DESTINATION ${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME})
endif()

# Undefined behavior sanitizer
if(${USE_UBSAN})
  foreach(sanitizer undefined float-divide-by-zero unsigned-integer-overflow
                    implicit-conversion local-bounds nullability)
    check_cxx_compiler_flag("-fsanitize=${sanitizer}"
                            CXX_COMPILER_HAS_sanitize_${sanitizer})
    if(CXX_COMPILER_HAS_sanitize_${sanitizer})
      target_compile_options(${a_lib_name} PRIVATE -fsanitize=${sanitizer} -g -O1)
      if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        target_link_options(${a_lib_name} PRIVATE -fsanitize=${sanitizer} -shared-libsan)
      else()
        target_link_options(${a_lib_name} PRIVATE -fsanitize=${sanitizer})
        target_link_libraries(${a_lib_name} PRIVATE ubsan)
      endif()
    endif()
  endforeach()
  install(FILES .github/ubsan.supp
          DESTINATION ${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME})
endif()

#-----------------------------------------------------------------------------------
# Parse jinja templates: once by default, and once for all yml files
set(TEMPLATE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/templates)
set(TEMPLATE_XML ${PROJECT_NAME}.xml.jinja2)

file(GLOB CONFIG_YMLS ${CMAKE_CURRENT_SOURCE_DIR}/configurations/*.yml)
foreach(config_yml ${CONFIG_YMLS})
    get_filename_component(config ${config_yml} NAME_WE)
    add_custom_target(${config}.xml ALL)
    add_custom_command(TARGET ${config}.xml POST_BUILD
        COMMAND
            ${CMAKE_CURRENT_SOURCE_DIR}/bin/make_detector_configuration
                --dir ${TEMPLATE_DIR}
                --template ${TEMPLATE_XML}
                --config ${config_yml}
                --output ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_${config}.xml
        COMMENT "Creating configuration ${config} for ${PROJECT_NAME}"
        )
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_${config}.xml
        DESTINATION ${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}/
        )
endforeach()

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}_full.xml
    DESTINATION ${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}/
    RENAME ${PROJECT_NAME}.xml
    )

#-----------------------------------------------------------------------------------
# Install the detector description files.
install(DIRECTORY compact/
    DESTINATION ${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}/compact
    FILES_MATCHING PATTERN "*.xml" PATTERN "*.gdml"
    )

#-----------------------------------------------------------------------------------
# Install the detector calibration files.
install(DIRECTORY calibrations/
    DESTINATION ${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}/calibrations
    )

#-----------------------------------------------------------------------------------
# Configure and install compact/version.xml
#
configure_file(templates/version.xml.in ${CMAKE_CURRENT_BINARY_DIR}/version.xml @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/version.xml
    DESTINATION ${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}/compact
)

#-----------------------------------------------------------------------------------
# Configure and install detector setup script
#
# When built with -DUSE_ASAN=ON / -DUSE_UBSAN=ON, embed sanitizer environment
# variable exports into thisepic.sh so that reusable workflows (and local users)
# automatically get the right suppression paths when they source this script.
# Suppression files are installed alongside the plugin under ${DETECTOR_PATH}.
# Note: bracket strings preserve ${...} literally (no CMake expansion); configure_file
# with @ONLY passes them through so the shell expands $DETECTOR_PATH at runtime.
set(ASAN_SETUP_LINES "")
set(UBSAN_SETUP_LINES "")
if(USE_ASAN)
  set(ASAN_SETUP_LINES [[
export ASAN_OPTIONS=suppressions=${DETECTOR_PATH}/asan.supp:malloc_context_size=20:detect_leaks=1:verify_asan_link_order=0:detect_stack_use_after_return=1:detect_odr_violation=1:new_delete_type_mismatch=0:intercept_tls_get_addr=0
export LSAN_OPTIONS=suppressions=${DETECTOR_PATH}/lsan.supp
]])
endif()
if(USE_UBSAN)
  set(UBSAN_SETUP_LINES [[
export UBSAN_OPTIONS=suppressions=${DETECTOR_PATH}/ubsan.supp:print_stacktrace=1:silence_unsigned_overflow=1:report_error_type=1
]])
endif()
configure_file(templates/thisepic.sh.in ${CMAKE_CURRENT_BINARY_DIR}/thisepic.sh @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/thisepic.sh
    DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(FILES cmake/setup.sh
    DESTINATION ${CMAKE_INSTALL_PREFIX}
)

# install programs
install(PROGRAMS bin/g4MaterialScan_to_csv
  DESTINATION ${CMAKE_INSTALL_BINDIR}
)
