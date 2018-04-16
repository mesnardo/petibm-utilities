bin_PROGRAMS = \
	petibm-interpolate

petibm_interpolate_SOURCES = \
	interpolate.cpp
petibm_interpolate_CPPFLAGS = \
	-I$(top_srcdir)/include \
	$(PETSC_CPPFLAGS)
petibm_interpolate_LDADD = \
	$(top_builddir)/src/libpetibm-utilities.la \
	$(PETSC_LDFLAGS) $(PETSC_LIBS)
