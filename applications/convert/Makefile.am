bin_PROGRAMS = \
	petibm-convert

petibm_convert_SOURCES = \
	convert.cpp
petibm_convert_CPPFLAGS = \
	-I$(top_srcdir)/include \
	$(PETSC_CPPFLAGS)
petibm_convert_LDADD = \
	$(top_builddir)/src/libpetibm-utilities.la \
	$(PETSC_LDFLAGS) $(PETSC_LIBS)
