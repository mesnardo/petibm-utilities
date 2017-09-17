include ../Makefile.inc

all: vorticity2d vorticity3d

vorticity2d: vorticity2d.o  chkopts
	-$(CLINKER) -o vorticity2d $< $(PETSC_LIB)
	-$(RM) $<

vorticity3d: vorticity3d.o  chkopts
	-$(CLINKER) -o vorticity3d $< $(PETSC_LIB)
	$(RM) $<

clean-vorticity:
	-$(RM) vorticity2d
	-$(RM) vorticity3d

.PHONY: all clean-vorticity
