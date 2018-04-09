bin_PROGRAMS = \
	petibm-crop

petibm_crop_SOURCES = \
	crop.cpp
petibm_crop_CPPFLAGS = \
	-I$(top_srcdir)/include \
	$(PETSC_CPPFLAGS)
petibm_crop_LDADD = \
	$(top_builddir)/src/libpetibm-utilities.la \
	$(PETSC_LDFLAGS) $(PETSC_LIBS)
