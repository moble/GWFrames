# Make targets in ./Code
# build_message: prints helpful build documentation
# install_user:  just call the python user build process
# cpp:           make object files to use in other C++ programs
# clean:         cleans out the builds, etc., in Code
# MikeHappy:     self-explanatory
CODE_TARGETS := build_message install_user cpp clean MikeHappy

# Tell 'make' not to look for files with the following names
.PHONY : $(CODE_TARGETS) doc

.DEFAULT_GOAL := build_message

$(CODE_TARGETS):
	$(MAKE) -C Code $@

# This rebuilds the documentation, assuming doxygen is working
doc :
	$(MAKE) -C Docs
