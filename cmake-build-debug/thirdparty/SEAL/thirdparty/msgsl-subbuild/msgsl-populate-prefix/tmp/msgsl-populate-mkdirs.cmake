# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

# If CMAKE_DISABLE_SOURCE_CHANGES is set to true and the source directory is an
# existing directory in our source tree, calling file(MAKE_DIRECTORY) on it
# would cause a fatal error, even though it would be a no-op.
if(NOT EXISTS "/home/tl/HE3DB/cmake-build-debug/thirdparty/SEAL/thirdparty/msgsl-src")
  file(MAKE_DIRECTORY "/home/tl/HE3DB/cmake-build-debug/thirdparty/SEAL/thirdparty/msgsl-src")
endif()
file(MAKE_DIRECTORY
  "/home/tl/HE3DB/cmake-build-debug/thirdparty/SEAL/thirdparty/msgsl-build"
  "/home/tl/HE3DB/cmake-build-debug/thirdparty/SEAL/thirdparty/msgsl-subbuild/msgsl-populate-prefix"
  "/home/tl/HE3DB/cmake-build-debug/thirdparty/SEAL/thirdparty/msgsl-subbuild/msgsl-populate-prefix/tmp"
  "/home/tl/HE3DB/cmake-build-debug/thirdparty/SEAL/thirdparty/msgsl-subbuild/msgsl-populate-prefix/src/msgsl-populate-stamp"
  "/home/tl/HE3DB/cmake-build-debug/thirdparty/SEAL/thirdparty/msgsl-subbuild/msgsl-populate-prefix/src"
  "/home/tl/HE3DB/cmake-build-debug/thirdparty/SEAL/thirdparty/msgsl-subbuild/msgsl-populate-prefix/src/msgsl-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/tl/HE3DB/cmake-build-debug/thirdparty/SEAL/thirdparty/msgsl-subbuild/msgsl-populate-prefix/src/msgsl-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/tl/HE3DB/cmake-build-debug/thirdparty/SEAL/thirdparty/msgsl-subbuild/msgsl-populate-prefix/src/msgsl-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
