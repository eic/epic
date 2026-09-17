# Determine version from git describe cmake-lint: disable=C0103
#
# Sets two variables in the calling scope:
#   ${VERSION}       - CMake-compatible version (X.Y.Z or X.Y.Z.N), suitable
#                      for passing to the project() command.
#   ${VERSION}_FULL  - Full git-describe string (e.g. 26.06.0-10-ge600dc77d),
#                      suitable for display / install scripts.
macro(set_git_version VERSION)
  if(NOT Git_Found)
    find_package(Git)
  endif()

  if(GIT_EXECUTABLE)
    # Generate a git-describe version string from Git repository tags
    execute_process(
      COMMAND ${GIT_EXECUTABLE} describe --tags --dirty --match "*.*.*"
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
      OUTPUT_VARIABLE GIT_DESCRIBE_VERSION
      RESULT_VARIABLE GIT_DESCRIBE_ERROR_CODE
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(NOT GIT_DESCRIBE_ERROR_CODE)
      set(${VERSION}_FULL ${GIT_DESCRIBE_VERSION})
    endif()
  endif()

  # Final fallback: Just use a bogus version string that is semantically older
  # than anything else and spit out a warning to the developer.
  if(NOT DEFINED ${VERSION}_FULL)
    set(${VERSION}_FULL "0.0.0-unknown")
    message(
      WARNING
        "Failed to determine VERSION from Git tags. Using default version \"${${VERSION}_FULL}\"."
    )
  endif()

  # Strip a leading "v" (e.g. v1.2.3 -> 1.2.3)
  string(REGEX REPLACE "^v" "" _git_version_stripped "${${VERSION}_FULL}")

  # Extract only the leading X.Y.Z or X.Y.Z.N numeric part so the result is a
  # valid CMake version string that can be passed to project().
  string(REGEX MATCH "^([0-9]+\\.[0-9]+(\\.[0-9]+(\\.[0-9]+)?)?)" _git_version_clean "${_git_version_stripped}")

  if(_git_version_clean)
    set(${VERSION} "${_git_version_clean}")
  else()
    set(${VERSION} "0.0.0")
    message(
      WARNING
        "Could not parse a numeric version from \"${${VERSION}_FULL}\". Falling back to \"${${VERSION}}\"."
    )
  endif()
  unset(_git_version_stripped)
  unset(_git_version_clean)
endmacro()
