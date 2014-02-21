# Tell 'make' not to look for files with the following names
.PHONY : all cpp doc clean

# By default, just call the python build process
all :
	$(MAKE) -C Code

# If needed, we can also $(MAKE) object files to use in other C++ programs
cpp :
	$(MAKE) -C Code cpp

# This rebuilds the documentation, assuming doxygen is working
doc :
	$(MAKE) -C Docs

# This just cleans out the builds, etc., in Code
clean :
	$(MAKE) -C Code clean

MikeHappy :
	$(MAKE) -C Code MikeHappy
