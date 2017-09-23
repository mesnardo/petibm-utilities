bin_PROGRAMS = petibm-interpolation2d \
               petibm-interpolation3d \
               petibm-interpolation2d-test \
               petibm-interpolation3d-test

petibm_interpolation2d_SOURCES = interpolation2d.cpp
petibm_interpolation2d_CPPFLAGS = -I$(top_srcdir)/include
petibm_interpolation2d_LDADD = $(top_builddir)/src/libpetibm-utilities.la

petibm_interpolation3d_SOURCES = interpolation3d.cpp
petibm_interpolation3d_CPPFLAGS = -I$(top_srcdir)/include
petibm_interpolation3d_LDADD = $(top_builddir)/src/libpetibm-utilities.la

petibm_interpolation2d_test_SOURCES = test2d.cpp
petibm_interpolation2d_test_CPPFLAGS = -I$(top_srcdir)/include
petibm_interpolation2d_test_LDADD = $(top_builddir)/src/libpetibm-utilities.la

petibm_interpolation3d_test_SOURCES = test3d.cpp
petibm_interpolation3d_test_CPPFLAGS = -I$(top_srcdir)/include
petibm_interpolation3d_test_LDADD = $(top_builddir)/src/libpetibm-utilities.la
