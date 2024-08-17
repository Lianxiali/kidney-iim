## ---------------------------------------------------------------------
##
## Copyright (c) 2020 - 2020 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------

PROJECT(main3d)

CMAKE_MINIMUM_REQUIRED(VERSION 3.15.0)

SET(IBAMR_ROOT /home/michael/software/IIM/IBAMR/build)


# change the path of the executable 
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR})

# List your source files here - this example has just one.
SET(SOURCE_FILES example.cpp FeedbackForcer.cpp VelocityBcCoefs.cpp BcData.cpp)
ADD_EXECUTABLE(${PROJECT_NAME} ${SOURCE_FILES})

FIND_PACKAGE(IBAMR 0.15.0 REQUIRED)
TARGET_LINK_LIBRARIES(${PROJECT_NAME} IBAMR::IBAMR3d)
# IBAMR saves the flags it used to compile - you can reuse them if you want to
SET(CMAKE_CXX_FLAGS ${IBAMR_CXX_FLAGS})


# SET(IBAMR_DIR $HOME/software/IBAMR/debug)

SET(CMAKE_BUILD_TYPE "Debug")
# SET(CMAKE_LIBRARY_PATH /home/michael/software/IBAMR/lib)

