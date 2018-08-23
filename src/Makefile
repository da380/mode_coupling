BINDIR=../bin
OBJDIR=../obj
MODDIR=../mod

CC=gcc
F95=gfortran
COMP_F95_FLAGS=-ffixed-line-length-none -fbounds-check -mcmodel=medium -finit-local-zero
MODFLAG=-J$(MODDIR)
OPTFLAGS=-O3
F77FLAGS=
#LAPACK_LINK=-L/opt/acml4.4.0/gfortran64_mp/lib -lacml_mp -lacml_mv
LAPACK_LINK= -L/home/chaos/moulik/libtage -llapack -lblas

CFLAGS=$(OPTFLAGS)
F95FLAGS=$(COMP_F95_FLAGS) $(OPTFLAGS)

CCOMPILE=$(CC) -c $(CFLAGS)
F95COMPILE=$(F95) -c $(MODFLAG) $(F95FLAGS)
F77COMPILE=$(F95) -c $(F77FLAGS) $(F95FLAGS)
F95LINK=$(F95) $(MODFLAG) $(F95FLAGS)

$(OBJDIR)/%.o:%.f
	$(F77COMPILE) -o $@ $+
$(OBJDIR)/%.o:%.F
	$(F77COMPILE) -o $@ $+
$(OBJDIR)/%.o:%.f90
	$(F95COMPILE) -o $@ $+
$(OBJDIR)/%.o:%.c
	$(CCOMPILE) -o $@ $+


.PHONY: all mode_list mdcpl_sp mdcpl dspec_pre dspec_cal dspec_pro QnmS

all: mode_list mdcpl_sp mdcpl dspec_pre dspec_cal dspec_pro QnmS

mode_list: $(BINDIR)/mode_list
mdcpl_sp: $(BINDIR)/mdcpl_sp
mdcpl: $(BINDIR)/mdcpl
dspec_pre: $(BINDIR)/dspec_pre
dspec_cal: $(BINDIR)/dspec_cal
dspec_pro: $(BINDIR)/dspec_pro
QnmS: $(BINDIR)/QnmS





$(BINDIR)/mode_list: $(OBJDIR)/csubs.o $(OBJDIR)/fgetenv.o $(OBJDIR)/getgnl.o $(OBJDIR)/getunx.o $(OBJDIR)/chekgl.o $(OBJDIR)/chekcl.o $(OBJDIR)/istlen.o $(OBJDIR)/cat2s.o $(OBJDIR)/check.o $(OBJDIR)/rcreao.o $(OBJDIR)/ffstat.o $(OBJDIR)/rewtp.o $(OBJDIR)/opnflc_cut.o $(OBJDIR)/opnfil.o $(OBJDIR)/openfl.o $(OBJDIR)/byswap.o $(OBJDIR)/bffis.o $(OBJDIR)/bffi.o $(OBJDIR)/seteig.o $(OBJDIR)/rspln.o $(OBJDIR)/krunge.o $(OBJDIR)/modl.o $(OBJDIR)/fetch_dum.o $(OBJDIR)/indexx.o $(OBJDIR)/mode_list.o
	$(F95LINK) -o $@ $+


$(BINDIR)/mdcpl_sp:  $(OBJDIR)/csubs.o \
	$(OBJDIR)/istlen.o \
	$(OBJDIR)/cat2s.o \
	$(OBJDIR)/check.o \
	$(OBJDIR)/rcreao.o \
	$(OBJDIR)/ffstat.o \
	$(OBJDIR)/rewtp.o \
	$(OBJDIR)/nrtype.o \
	$(OBJDIR)/opnflc_cut.o \
	$(OBJDIR)/opnfil.o \
	$(OBJDIR)/openfl.o \
	$(OBJDIR)/byswap.o \
	$(OBJDIR)/bffis.o \
	$(OBJDIR)/bffi.o \
	$(OBJDIR)/seteig.o \
	$(OBJDIR)/rspln.o \
	$(OBJDIR)/krunge.o \
	$(OBJDIR)/modl.o \
	$(OBJDIR)/wig0.o \
	$(OBJDIR)/bcoff1.o \
	$(OBJDIR)/wig2.o \
	$(OBJDIR)/grinta.o \
	$(OBJDIR)/grintb.o \
	$(OBJDIR)/corfac.o \
	$(OBJDIR)/gravk.o \
	$(OBJDIR)/memdcpl_new.o \
	$(OBJDIR)/fetch.o \
	$(OBJDIR)/rdmdl.o \
	$(OBJDIR)/heread.o \
	$(OBJDIR)/mskmdl.o \
	$(OBJDIR)/rsple.o \
	$(OBJDIR)/splh.o \
	$(OBJDIR)/splhsetup.o \
	$(OBJDIR)/intmdcpl_new.o \
	$(OBJDIR)/intpltnnew.o \
	$(OBJDIR)/intpltnnew_sep.o \
	$(OBJDIR)/fgetenv.o \
	$(OBJDIR)/getgnl.o \
	$(OBJDIR)/chekgl.o \
	$(OBJDIR)/getunx.o \
	$(OBJDIR)/chekcl.o \
	$(OBJDIR)/mdcpl_sp.o
	$(F95LINK) -o $@ $+


$(BINDIR)/mdcpl:  $(OBJDIR)/csubs.o \
	$(OBJDIR)/istlen.o \
	$(OBJDIR)/cat2s.o \
	$(OBJDIR)/check.o \
	$(OBJDIR)/rcreao.o \
	$(OBJDIR)/ffstat.o \
	$(OBJDIR)/rewtp.o \
	$(OBJDIR)/nrtype.o \
	$(OBJDIR)/opnflc_cut.o \
	$(OBJDIR)/opnfil.o \
	$(OBJDIR)/openfl.o \
	$(OBJDIR)/byswap.o \
	$(OBJDIR)/bffis.o \
	$(OBJDIR)/bffi.o \
	$(OBJDIR)/seteig.o \
	$(OBJDIR)/rspln.o \
	$(OBJDIR)/krunge.o \
	$(OBJDIR)/modl.o \
	$(OBJDIR)/wig0.o \
	$(OBJDIR)/bcoff1.o \
	$(OBJDIR)/wig2.o \
	$(OBJDIR)/grinta.o \
	$(OBJDIR)/grintb.o \
	$(OBJDIR)/corfac.o \
	$(OBJDIR)/gravk.o \
	$(OBJDIR)/memdcpl_new.o \
	$(OBJDIR)/fetch.o \
	$(OBJDIR)/rdmdl.o \
	$(OBJDIR)/heread.o \
	$(OBJDIR)/mskmdl.o \
	$(OBJDIR)/rsple.o \
	$(OBJDIR)/splh.o \
	$(OBJDIR)/splhsetup.o \
	$(OBJDIR)/intmdcpl_new.o \
	$(OBJDIR)/intpltnnew.o \
	$(OBJDIR)/intpltnnew_sep.o \
	$(OBJDIR)/fgetenv.o \
	$(OBJDIR)/getgnl.o \
	$(OBJDIR)/chekgl.o \
	$(OBJDIR)/getunx.o \
	$(OBJDIR)/chekcl.o \
	$(OBJDIR)/mdcpl.o
	$(F95LINK) -o $@ $+

$(BINDIR)/dspec_pre: $(OBJDIR)/nrtype.o \
	$(OBJDIR)/nrutil.o \
	$(OBJDIR)/module_util.o \
	$(OBJDIR)/module_fourier.o \
	$(OBJDIR)/dspec_pre.o
	$(F95LINK) -o $@ $+

$(BINDIR)/QnmS: $(OBJDIR)/nrtype.o \
	$(OBJDIR)/nrutil.o \
	$(OBJDIR)/module_util.o \
	$(OBJDIR)/module_fourier.o \
	$(OBJDIR)/fgetenv.o \
	$(OBJDIR)/getgnl.o \
	$(OBJDIR)/chekgl.o \
	$(OBJDIR)/getunx.o \
	$(OBJDIR)/chekcl.o \
	$(OBJDIR)/csubs.o \
	$(OBJDIR)/rdsource.o \
	$(OBJDIR)/julday.o \
	$(OBJDIR)/lpyr.o \
	$(OBJDIR)/QnmS.o
	$(F95LINK) -o $@ $+


$(BINDIR)/dspec_cal: $(OBJDIR)/csubs.o \
	$(OBJDIR)/istlen.o \
	$(OBJDIR)/cat2s.o \
	$(OBJDIR)/check.o \
	$(OBJDIR)/rcreao.o \
	$(OBJDIR)/ffstat.o \
	$(OBJDIR)/opnflc_cut.o \
	$(OBJDIR)/opnfil.o \
	$(OBJDIR)/openfl.o \
	$(OBJDIR)/byswap.o \
	$(OBJDIR)/bffis.o \
	$(OBJDIR)/bffi.o \
	$(OBJDIR)/bffibs4.o \
	$(OBJDIR)/fgetenv.o \
	$(OBJDIR)/getgnl.o \
	$(OBJDIR)/chekgl.o \
	$(OBJDIR)/getunx.o \
	$(OBJDIR)/chekcl.o \
	$(OBJDIR)/gpremlib.o \
	$(OBJDIR)/csveca.o \
	$(OBJDIR)/crvecd.o \
	$(OBJDIR)/soprema.o \
	$(OBJDIR)/reprema.o \
	$(OBJDIR)/module_util.o \
	$(OBJDIR)/dspec_cal.o 
	$(F95LINK) -o $@ $+ $(LAPACK_LINK)


$(BINDIR)/dspec_pro: \
	$(OBJDIR)/csubs.o \
	$(OBJDIR)/istlen.o \
	$(OBJDIR)/cat2s.o \
	$(OBJDIR)/check.o \
	$(OBJDIR)/rcreao.o \
	$(OBJDIR)/ffstat.o \
	$(OBJDIR)/nrutil.o \
	$(OBJDIR)/module_util.o \
	$(OBJDIR)/module_fourier.o \
	$(OBJDIR)/gpremlib.o \
	$(OBJDIR)/opnflc_cut.o \
	$(OBJDIR)/opnfil.o \
	$(OBJDIR)/openfl.o \
	$(OBJDIR)/byswap.o \
	$(OBJDIR)/bffis.o \
	$(OBJDIR)/bffi.o \
	$(OBJDIR)/bffibs4.o \
	$(OBJDIR)/fgetenv.o \
	$(OBJDIR)/getgnl.o \
	$(OBJDIR)/chekgl.o \
	$(OBJDIR)/getunx.o \
	$(OBJDIR)/chekcl.o \
	$(OBJDIR)/dspec_pro.o
	$(F95LINK) -o $@ $+  $(LAPACK_LINK)

