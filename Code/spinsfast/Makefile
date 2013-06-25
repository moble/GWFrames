#################################################################
#
#  Before running make, set up a congfiguration file in build/
#  then set environment variable "build" to point at that file.
#
#################################################################


# define the software modules to be built
MODULES := code obj example python

# location for programs / object files
EXECDIR := bin
OBJDIR := obj
LIBDIR:=lib

ifndef build
  build:=$(error Environment variable "build" undefined, should point to configuration file in build/, as in "export build=build/config.mk")UNDEFINED
endif


# Primary target is all, machine specific info from env. var. $build
all : alltargets $(build)


# Use pattern sub-string replacement to include module directories
INC := $(patsubst %,-I%, $(MODULES))
INC += -Iinclude

# Include machine specifics
include $(build)

# Variables to hold targets
EXEC :=  
OBJ := 
LIB := 
# Get the module.mk files and include them, collecting targets and subprogs
MODMK := $(patsubst %,%/module.mk, $(MODULES))
include $(MODMK)

TARGETS := $(EXEC) $(LIB) $(OBJ) 


alltargets : $(TARGETS)

clean : 
	rm $(TARGETS)
