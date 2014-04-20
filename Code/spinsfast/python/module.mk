LIB+=$(LIBDIR)/spinsfast.so

$(LIBDIR)/spinsfast.so : python/spinsfast_module.c python/setup.py $(LIBDIR)/libspinsfast.a
	python python/setup.py build --build-base $(OBJDIR) install --install-lib $(LIBDIR)