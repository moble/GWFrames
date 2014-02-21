# Tell 'make' not to look for files with the following names
.PHONY : all cpp doc

# By default, just call the python build process
all :
	$(MAKE) -C Code

# If needed, we can also make object files to use in other C++ programs
cpp :
	$(MAKE) -C Code cpp

# This rebuilds the documentation, assuming doxygen is working
doc :
	$(MAKE) -C Docs


MikeHappy :
	$(MAKE) -C Code MikeHappy
