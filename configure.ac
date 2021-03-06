# file: configure.ac
# author: Olivier Mesnard (mesnardo@gwu.edu)
# brief: Input file for Autoconf.

# set-up
AC_PREREQ([2.69])
AC_INIT([petibm-utilities], [0.1dev], [mesnardo@gwu.edu])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign subdir-objects])
AC_CONFIG_SRCDIR([src/field.cpp])
AC_CONFIG_HEADERS([config/config.h])
AM_SILENT_RULES([no])
AM_MAINTAINER_MODE([disable])
AC_LANG([C++])
LT_PREREQ([2.2])
LT_INIT([dlopen])

# check for programs
AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_RANLIB
AC_PROG_INSTALL
AC_PROG_MKDIR_P
AC_PROG_MAKE_SET
AC_PROG_LN_S

# check for typedefs, structures, and compiler characteristics

# check for header files

# check for library functions


AC_SUBST([BUILDDIR], [`pwd`])
CPPFLAGS_="-I$prefix/include"
LDFLAGS_="-L$prefix/lib -Wl,-rpath,$prefix/lib"
CPPFLAGS_PREPEND($CPPFLAGS_)
LDFLAGS_PREPEND($LDFLAGS_)

CONFIGURE_PETSC

# check version of PetIBM
AC_ARG_WITH([petibm],
            AS_HELP_STRING([--with-petibm=VERSION],
                           [set PetIBM version (0.2 or 0.3)]),
            [PETIBM_VERSION=$withval],
            [PETIBM_VERSION=0.3])
if test "$PETIBM_VERSION" = "0.3"; then
  CPPFLAGS_APPEND("-DPETIBM_0_3")
elif test "$PETIBM_VERSION" = "0.2"; then
  CPPFLAGS_APPEND("-DPETIBM_0_2")
fi

# list of Makefiles to generate
AC_CONFIG_FILES([Makefile
                 include/Makefile
                 src/Makefile
                 applications/Makefile
                 applications/convert/Makefile
                 applications/crop/Makefile
                 applications/interpolate/Makefile
                 applications/vorticity/Makefile])

# output message
AC_OUTPUT
echo
echo "====================="
echo "Configuration summary"
echo "====================="
AC_MSG_NOTICE([PetIBM-utilities appears to have been configured successfully

Prefix: '${prefix}'

CXX: ${CXX}
CXXFLAGS: ${CXXFLAGS}

CPPFLAGS: ${CPPFLAGS}

LDFLAGS: ${LDFLAGS}

LIBS: ${LIBS}

PETSC_DIR: ${PETSC_DIR}
PETSC_ARCH: ${PETSC_ARCH}

Now type 'make @<:@<target>@:>@'
  where the optional @<:@<target>@:>@ is:
    all            - build all binaries
    install        - install everything
])
