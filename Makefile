# Multiple-machine Makefile
#   puts machine-specific *.o files in Obj sub-directories
#   1st-style targets requires a "make" command with pattern-matching 
#     capability, e.g. "make" command on Suns and Gnu "gmake"
#   2nd-style targets should work with "make" on any Unix machine

SHELL = /bin/sh
.IGNORE:

#
# Files
#

ROOT =		lj

SRC_a =		lja.f
SRC_f =		ljf.f
SRC_s =		ljs.f

SRC_a_unix =	$(SRC_a) parliba_unix.f
SRC_a_ncube =	$(SRC_a) parliba_ncube.f
SRC_a_gamma =	$(SRC_a) parliba_gamma.f
SRC_a_mos =	$(SRC_a) parliba_mos.f
SRC_a_osf =	$(SRC_a) parliba_osf.f
SRC_a_t3d =	parliba_t3d.f
SRC_a_mpi =	$(SRC_a) parliba_mpi.f

SRC_f_unix =	$(SRC_f) parlibf_unix.f
SRC_f_ncube =	$(SRC_f) parlibf_ncube.f
SRC_f_gamma =	$(SRC_f) parlibf_gamma.f
SRC_f_mos =	$(SRC_f) parlibf_mos.f
SRC_f_osf =	$(SRC_f) parlibf_osf.f
SRC_f_t3d =	parlibf_t3d.f
SRC_f_mpi =	$(SRC_f) parlibf_mpi.f

SRC_s_unix =	$(SRC_s) parlibs_unix.f
SRC_s_ncube =	$(SRC_s) parlibs_ncube.f
SRC_s_gamma =	$(SRC_s) parlibs_gamma.f
SRC_s_mos =	$(SRC_s) parlibs_mos.f
SRC_s_osf =	$(SRC_s) parlibs_osf.f
SRC_s_t3d =	parlibs_t3d.f
SRC_s_mpi =	$(SRC_s) parlibs_mpi.f

#
# Machine-specific definitions
#

unix_DEFS =	"F77 = f77" \
		"F77FLAGS = -O" \
		"CC = cc" \
		"CCFLAGS = -O" \
		"LINK = f77" \
		"LINKFLAGS = -O" \
		"USRLIB = " \
		"SYSLIB = " \
		"SIZE = size"

ncube_DEFS =	"F77 = $(NCUBE_PATH)/ncc" \
		"F77FLAGS = -O" \
		"CC = $(NCUBE_PATH)/ncc" \
		"CCFLAGS = -O" \
		"LINK = $(NCUBE_PATH)/ncc" \
		"LINKFLAGS = -g -Ncomm 200000 -Nheap 65536" \
		"USRLIB = " \
		"SYSLIB = -lf -lm -ln_lib" \
		"SIZE = $(NCUBE_PATH)/nsize"

gamma_DEFS =	"F77 = $(GAMMA_PATH)/if77" \
		"F77FLAGS = -O3 -Knoieee" \
		"CC = $(GAMMA_PATH)/icc" \
		"CCFLAGS = -O3 -Knoieee" \
		"LINK = $(GAMMA_PATH)/if77" \
		"LINKFLAGS = -Knoieee -node" \
		"USRLIB = " \
		"SYSLIB = " \
		"SIZE = $(GAMMA_PATH)/size860"

mos_DEFS =	"F77 = $(MOS_PATH)/sif77" \
		"F77FLAGS = -O3 -Knoieee" \
		"CC = $(MOS_PATH)/sicc" \
		"CCFLAGS = -O3 -Knoieee" \
		"LINK = $(MOS_PATH)/sif77" \
		"LINKFLAGS = -Knoieee" \
		"USRLIB = " \
		"SYSLIB = " \
		"SIZE = $(PARAGON_PATH)/size860"

osf_DEFS =	"F77 = $(PARAGON_PATH)/if77" \
		"F77FLAGS = -O3 -Knoieee" \
		"CC = $(PARAGON_PATH)/icc" \
		"CCFLAGS = -O3 -Knoieee" \
		"LINK = $(PARAGON_PATH)/if77" \
		"LINKFLAGS = -Knoieee -nx" \
		"USRLIB = " \
		"SYSLIB = " \
		"SIZE = $(PARAGON_PATH)/size860"

t3d_DEFS =	"F77 = /mpp/bin/cft77" \
		"F77FLAGS = -Ccray-t3d -dp -oaggress" \
		"CC = /mpp/bin/cc" \
		"CCFLAGS = -Tcray-t3d -c" \
		"LINK = /mpp/bin/mppldr" \
		"LINKFLAGS = " \
		"USRLIB = " \
		"SYSLIB = " \
		"SIZE = /mpp/bin/mppsize"

mpi_DEFS =	"F77 = $(MOS_PATH)/sif77" \
		"F77FLAGS = -O3 -Knoieee" \
		"CC = $(MOS_PATH)/sicc" \
		"CCFLAGS = -O3 -Knoieee" \
		"LINK = $(MOS_PATH)/sif77" \
		"LINKFLAGS = -Knoieee" \
		"USRLIB = " \
		"SYSLIB = -lmpi" \
		"SIZE = $(PARAGON_PATH)/size860"

#
# Default
#

help:
	@echo 'Type "(g)make prefix.suffix"
	@echo '   where prefix is one of:'
	@echo '      a		(atom decomposition)'
	@echo '      f		(force decomposition)'
	@echo '      s		(spatial decomposition)'
	@echo '   and suffix is one of:'
	@echo '      unix	(for single processor - workstation or vector)'
	@echo '      ncube	(for nCUBE 2)'
	@echo '      gamma	(for Intel iPSC/860)'
	@echo '      mos	(for Intel Paragon w/ SUNMOS)'
	@echo '      osf	(for Intel Paragon w/ OSF)'
	@echo '      t3d	(for Cray T3D)'
	@echo '      mpi	(for universal MPI)'

#
# Parallel Targets
#

a.unix:
	@$(MAKE) "PREFIX = a" "SUFFIX = unix" recurse

a.ncube:
	@$(MAKE) "PREFIX = a" "SUFFIX = ncube" recurse

a.gamma:
	@$(MAKE) "PREFIX = a" "SUFFIX = gamma" recurse

a.mos:
	@$(MAKE) "PREFIX = a" "SUFFIX = mos" recurse

a.osf:
	@$(MAKE) "PREFIX = a" "SUFFIX = osf" recurse

a.mpi:
	@$(MAKE) "PREFIX = a" "SUFFIX = mpi" recurse

f.unix:
	@$(MAKE) "PREFIX = f" "SUFFIX = unix" recurse

f.ncube:
	@$(MAKE) "PREFIX = f" "SUFFIX = ncube" recurse

f.gamma:
	@$(MAKE) "PREFIX = f" "SUFFIX = gamma" recurse

f.mos:
	@$(MAKE) "PREFIX = f" "SUFFIX = mos" recurse

f.osf:
	@$(MAKE) "PREFIX = f" "SUFFIX = osf" recurse

f.mpi:
	@$(MAKE) "PREFIX = f" "SUFFIX = mpi" recurse

s.unix:
	@$(MAKE) "PREFIX = s" "SUFFIX = unix" recurse

s.ncube:
	@$(MAKE) "PREFIX = s" "SUFFIX = ncube" recurse

s.gamma:
	@$(MAKE) "PREFIX = s" "SUFFIX = gamma" recurse

s.mos:
	@$(MAKE) "PREFIX = s" "SUFFIX = mos" recurse

s.osf:
	@$(MAKE) "PREFIX = s" "SUFFIX = osf" recurse

s.mpi:
	@$(MAKE) "PREFIX = s" "SUFFIX = mpi" recurse

a.t3d:
	@if [ ! -d Obj_t3d ]; then mkdir Obj_t3d; fi
	@cp -p *.[fh] Makefile_lower Obj_t3d
	@cd Obj_t3d; \
	$(MAKE) -f Makefile_lower $(t3d_DEFS) \
		"OBJDIR = Obj_t3d" \
		"OBJ = $(SRC_a:.f=.o) \
		       $(SRC_a_t3d:.f=.o)" \
		"EXE = ../$(ROOT)a_t3d" \
		../$(ROOT)a_t3d
	@rm Obj_t3d/*.[fh] Obj_t3d/Makefile_lower

f.t3d:
	@if [ ! -d Obj_t3d ]; then mkdir Obj_t3d; fi
	@cp -p *.[fh] Makefile_lower Obj_t3d
	@cd Obj_t3d; \
	$(MAKE) -f Makefile_lower $(t3d_DEFS) \
		"OBJDIR = Obj_t3d" \
		"OBJ = $(SRC_f:.f=.o) \
		       $(SRC_f_t3d:.f=.o)" \
		"EXE = ../$(ROOT)f_t3d" \
		../$(ROOT)f_t3d
	@rm Obj_t3d/*.[fh] Obj_t3d/Makefile_lower

s.t3d:
	@if [ ! -d Obj_t3d ]; then mkdir Obj_t3d; fi
	@cp -p *.[fh] Makefile_lower Obj_t3d
	@cd Obj_t3d; \
	$(MAKE) -f Makefile_lower $(t3d_DEFS) \
		"OBJDIR = Obj_t3d" \
		"OBJ = $(SRC_s:.f=.o) \
		       $(SRC_s_t3d:.f=.o)" \
		"EXE = ../$(ROOT)s_t3d" \
		../$(ROOT)s_t3d
	@rm Obj_t3d/*.[fh] Obj_t3d/Makefile_lower

#
# Cleans
#

clean_all:
	rm -r Obj_*

clean_unix:
	rm -r Obj_unix

clean_ncube:
	rm -r Obj_ncube

clean_gamma:
	rm -r Obj_gamma

clean_mos:
	rm -r Obj_mos

clean_osf:
	rm -r Obj_osf

clean_t3d:
	rm -r Obj_t3d

clean_mpi:
	rm -r Obj_mpi

#
# Recursive make
#

recurse:
	@if [ ! -d Obj_$(SUFFIX) ]; then mkdir Obj_$(SUFFIX); fi
	@$(MAKE) $($(SUFFIX)_DEFS) \
		"OBJDIR = Obj_$(SUFFIX)" \
		"OBJ = $(SRC_$(PREFIX)_$(SUFFIX):%.f=Obj_$(SUFFIX)/%.o)" \
		"EXE = $(ROOT)$(PREFIX)_$(SUFFIX)" \
		$(ROOT)$(PREFIX)_$(SUFFIX)

EXE =	NULL
$(EXE): $(OBJ)
	$(LINK) $(LINKFLAGS) $(OBJ) $(USRLIB) $(SYSLIB) -o $@
	$(SIZE) $@

#
# Compilation rules
#

$(OBJDIR)/%.o:%.f
	$(F77) $(F77FLAGS) -c -o $@ $<

$(OBJDIR)/%.o:%.c
	$(CC) $(CCFLAGS) -c -o $@ $<

#
# Individual dependencies
#

$(OBJDIR)/lja.o:	lja.h
$(OBJDIR)/ljf.o:	ljf.h
$(OBJDIR)/ljs.o:	ljs.h
