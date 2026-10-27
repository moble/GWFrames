LIB+=$(LIBDIR)/libspinsfast.a

DEP_libspinsfast:=$(OBJ)

#DEP_libspinsfast:=$(OBJDIR)/alm.o  $(OBJDIR)/spinsfast_backward_Gmm.o  $(OBJDIR)/spinsfast_backward_transform.o  $(OBJDIR)/spinsfast_forward_Imm.o $(OBJDIR)/spinsfast_forward_Jmm.o  $(OBJDIR)/spinsfast_forward_transform.o  $(OBJDIR)/wigner_d_halfpi.o


$(LIBDIR)/libspinsfast.a : $(DEP_libspinsfast)
	ar -r $@ $(DEP_libspinsfast)
