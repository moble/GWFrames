# Tell 'make' not to look for files with the following names
.PHONY : all cpp doc clean

# By default, just call the python build process
all :
	make -C Code

# If needed, we can also make object files to use in other C++ programs
cpp :
	make -C Code cpp

# This rebuilds the documentation, assuming doxygen is working
doc :
	make -C Docs

# This just cleans out the builds, etc., in Code
clean :
	make -C Code clean
