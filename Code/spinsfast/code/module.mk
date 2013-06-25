OBJ+=$(OBJDIR)/alm.o $(OBJDIR)/wigner_d_halfpi_Risbo.o $(OBJDIR)/wigner_d_halfpi_TN.o $(OBJDIR)/wigner_d_halfpi_methods.o $(OBJDIR)/spinsfast_backward_transform.o $(OBJDIR)/spinsfast_backward_Gmm.o  $(OBJDIR)/spinsfast_backward_Gmm_alm2iqu.o $(OBJDIR)/spinsfast_forward_transform.o $(OBJDIR)/spinsfast_forward_Imm.o  $(OBJDIR)/spinsfast_forward_Jmm.o $(OBJDIR)/spinsfast_forward_transform_from_Imm.o $(OBJDIR)/spinsfast_forward_transform_eo.o $(OBJDIR)/spinsfast_forward_transform_iqu2alm.o $(OBJDIR)/dump_cimage.o $(OBJDIR)/healpix_convert.o


$(OBJDIR)/alm.o : code/alm.c include/alm.h
	$(CC) $(INC) -o $@ -c $<

$(OBJDIR)/dump_cimage.o : code/dump_cimage.c
	$(CC) $(INC) -o $@ -c $<

$(OBJDIR)/wigner_d_halfpi_Risbo.o : code/wigner_d_halfpi_Risbo.c include/wigner_d_halfpi.h include/wigner_d_halfpi_Risbo.h
	$(CC) $(INC) -o $@ -c $<

$(OBJDIR)/wigner_d_halfpi_TN.o : code/wigner_d_halfpi_TN.c include/wigner_d_halfpi.h include/wigner_d_halfpi_TN.h
	$(CC) $(INC) -o $@ -c $<

$(OBJDIR)/wigner_d_halfpi_methods.o : code/wigner_d_halfpi_methods.c include/wigner_d_halfpi.h include/wigner_d_halfpi_methods.h
	$(CC) $(INC) -o $@ -c $<

$(OBJDIR)/spinsfast_backward_transform.o : code/spinsfast_backward_transform.c include/spinsfast_backward.h
	$(CC) $(INC) -o $@ -c $<

$(OBJDIR)/spinsfast_backward_Gmm.o : code/spinsfast_backward_Gmm.c  include/spinsfast_backward.h
	$(CC) $(INC) -o $@ -c $<

$(OBJDIR)/spinsfast_backward_Gmm_alm2iqu.o : code/spinsfast_backward_Gmm_alm2iqu.c  include/spinsfast_backward.h
	$(CC) $(INC) -o $@ -c $<

$(OBJDIR)/spinsfast_forward_transform.o : code/spinsfast_forward_transform.c include/spinsfast_forward.h include/wigner_d_halfpi*.h 
	$(CC) $(INC) -o $@ -c $<

$(OBJDIR)/spinsfast_forward_transform_eo.o : code/spinsfast_forward_transform_eo.c include/spinsfast_forward.h include/wigner_d_halfpi*.h 
	$(CC) $(INC) -o $@ -c $<

$(OBJDIR)/spinsfast_forward_transform_iqu2alm.o : code/spinsfast_forward_transform_iqu2alm.c include/spinsfast_forward.h include/wigner_d_halfpi*.h 
	$(CC) $(INC) -o $@ -c $<

$(OBJDIR)/spinsfast_forward_transform_from_Imm.o : code/spinsfast_forward_transform_from_Imm.c include/spinsfast_forward.h include/wigner_d_halfpi*.h 
	$(CC) $(INC) -o $@ -c $<


$(OBJDIR)/spinsfast_forward_Imm.o : code/spinsfast_forward_Imm.c  include/spinsfast_forward.h
	$(CC) $(INC) -o $@ -c $<

$(OBJDIR)/spinsfast_forward_Jmm.o : code/spinsfast_forward_Jmm.c  include/spinsfast_forward.h
	$(CC) $(INC) -o $@ -c $<

 $(OBJDIR)/healpix_convert.o :  code/healpix_convert.c 
	$(CC) $(INC) -o $@ -c $<
