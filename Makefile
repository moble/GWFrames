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
	@echo ""
	@echo "========================================================"
	@echo "NOTE:"
	@echo "Unfortunately, this documentation system is not perfect."
	@echo "If SWIG reports an error when building the module, just"
	@echo "edit Docs/GWFrames_Doc.i, go to the line(s) in the error"
	@echo "message, and remove the offending clause(s)."
	@echo "========================================================"
