# Determine version from git describe cmake-lint: disable=C0103
macro(set_git_version VERSION)
  if(NOT Git_Found)
    find_package(Git)
  endif()

  if(GIT_EXECUTABLE)
    # Generate a git-describe version string from Git repository tags
    execute_process(
      COMMAND ${GIT_EXECUTABLE} describe --tags --dirty --match "v*"
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
      OUTPUT_VARIABLE GIT_DESCRIBE_VERSION
      RESULT_VARIABLE GIT_DESCRIBE_ERROR_CODE
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(NOT GIT_DESCRIBE_ERROR_CODE)
      set(${VERSION} ${GIT_DESCRIBE_VERSION})
    endif()
  endif()

  # Final fallback: Just use a bogus version string that is semantically older
  # than anything else and spit out a warning to the developer.
  if(NOT DEFINED ${VERSION})
    set(${VERSION} v0.0.0-unknown)
    message(
      WARNING
        "Failed to determine VERSION from Git tags. Using default version \"${${VERSION}}\"."
    )
  endif()
endmacro()
