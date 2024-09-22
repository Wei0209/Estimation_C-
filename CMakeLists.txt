CMAKE_MINIMUM_REQUIRED(VERSION 3.0.0)

## project name
project(hi)
# add_executable(hi2 hello.cpp fVector.h fVector.cpp)
# add_executable(PE1 testPE.cpp  fVector.cpp fMatrix.cpp ParamEstimator.cpp fVector.h fMatrix.h ParamEstimator.h)
add_executable(go VectorMatrixTestAP.cpp fMatrix.h fMatrix.cpp fVector.h fVector.cpp )


# ## Compiler settings
# set(CMAKE_BUILD_TYPE "Debug")
# #set(CMAKE_BUILD_TYPE "Release")

# IF (CMAKE_BUILD_TYPE STREQUAL "Debug")
# #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O0")
# set(CMAKE_CONFIGURATION_TYPES "Debug")
# ENDIF()


# IF (CMAKE_BUILD_TYPE STREQUAL "Release")
# #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -03")
# set(CMAKE_CONFIGURATION_TYPES "Release")
# ENDIF()
# ## End of compiler settings
# target_compile_definitions(foo PUBLIC DEBUG)

# ##Compiler output settings
# aux_source_directory(./ SRC_LIST)
# MESSAGE(STATUS "${SRC_LIST}")
# add_executable(ex1 hello1.cpp)

# #this part is so the Modern CMak book can verify this example builds.
# #For you code. you'll probably wants test too.
# enable_testing()
# add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
# add_executable(ex1 VectorMatrixTestAP.cpp fMatrix.h fVector.h)