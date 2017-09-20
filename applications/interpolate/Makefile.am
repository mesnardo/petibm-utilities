bin_PROGRAMS = petibm-interpolate2d \
               petibm-interpolate3d \
               petibm-interpolate2d-test \
               petibm-interpolate3d-test

petibm_interpolate2d_SOURCES = interpolate2d.cpp
petibm_interpolate2d_CPPFLAGS = -I$(top_srcdir)/include
petibm_interpolate2d_LDADD = $(top_builddir)/src/libpetibm-utilities.la

petibm_interpolate3d_SOURCES = interpolate3d.cpp
petibm_interpolate3d_CPPFLAGS = -I$(top_srcdir)/include
petibm_interpolate3d_LDADD = $(top_builddir)/src/libpetibm-utilities.la

petibm_interpolate2d_test_SOURCES = test2d.cpp
petibm_interpolate2d_test_CPPFLAGS = -I$(top_srcdir)/include
petibm_interpolate2d_test_LDADD = $(top_builddir)/src/libpetibm-utilities.la

petibm_interpolate3d_test_SOURCES = test3d.cpp
petibm_interpolate3d_test_CPPFLAGS = -I$(top_srcdir)/include
petibm_interpolate3d_test_LDADD = $(top_builddir)/src/libpetibm-utilities.la
