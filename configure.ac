#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(IMPERIAL, 1.0, mei@usc.edu)
AC_CONFIG_SRCDIR([src])
AM_INIT_AUTOMAKE()

# Checks for programs.
AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST

dnl 
dnl Setup environment so dependencies are used in build
dnl 
CFLAGS="$CFLAGS"
CPPFLAGS="-I$prefix/include $CPPFLAGS"
LDFLAGS="-L$prefix/lib -L$prefix/lib64 $LDFLAGS"

# Checks for library functions.
dnl PROJ
AC_ARG_ENABLE([proj],
    [AC_HELP_STRING([--enable-proj],
        [install Proj @<:@default=yes@:>@])],
        [if test "$enableval" = "yes" ; then install_proj=yes; else install_proj=no; fi],
        [install_proj=yes])
AM_CONDITIONAL([INSTALL_PROJ], [test "$install_proj" = yes])
AC_ARG_WITH([proj-incdir],
    [AC_HELP_STRING([--with-proj-incdir], [directory containing Proj header files @<:@default=no@:>@])],
        [with_proj_incdir="$withval"],
        [with_proj_incdir=no])
AC_SUBST(with_proj_incdir)
AC_ARG_WITH([proj-libdir],
    [AC_HELP_STRING([--with-proj-libdir], [directory containing Proj library @<:@default=no@:>@])],
        [with_proj_libdir="$withval"],
        [with_proj_libdir=no])
AC_SUBST(with_proj_libdir)

dnl proj (required)
if test "$install_proj" = yes ; then
  dnl sqlite (required by proj)
  if test "$install_sqlite" = yes ; then
    SQLITE3_INCLUDES="-I$prefix/include"
    SQLITE3_LDFLAGS="-L$prefix/lib"
  else
    if test "$with_sqlite_incdir" != no; then
      SQLITE3_INCLUDES="-I$with_sqlite_incdir"
    fi
    if test "$with_sqlite_libdir" != no; then
      SQLITE3_LDFLAGS="-L$with_sqlite_libdir"
    fi
    SCEC_SQLITE3_HEADER
    SCEC_SQLITE3_LIB
  fi
  AC_SUBST(SQLITE3_INCLUDES)
  AC_SUBST(SQLITE3_LDFLAGS)

  dnl tiff (required by proj)
  if test "$install_tiff" = yes ; then
    TIFF_INCLUDES="-I$prefix/include"
    TIFF_LDFLAGS="-L$prefix/lib"
  else
    if test "$with_tiff_incdir" != no; then
      TIFF_INCLUDES="-I$with_tiff_incdir"
    fi
    if test "$with_tiff_libdir" != no; then
      TIFF_LDFLAGS="-L$with_tiff_libdir"
    fi
    SCEC_TIFF_HEADER
    SCEC_TIFF_LIB
  fi
  AC_SUBST(TIFF_INCLUDES)
  AC_SUBST(TIFF_LDFLAGS)

  PROJ_INCLUDES="-I$prefix/include"
  PROJ_LDFLAGS="-L$prefix/lib"
else
  SCEC_PROJ_HEADER
  SCEC_PROJ_LIB
fi
if test "$with_proj_incdir" != no ; then
  PROJ_INCLUDES="-I$with_proj_incdir"
fi
if test "$with_proj_libdir" != no ; then
  PROJ_LDFLAGS="-L$with_proj_libdir -lproj"
fi
AC_SUBST(PROJ_INCLUDES)
AC_SUBST(PROJ_LDFLAGS)

AC_SUBST(CFLAGS)
AC_SUBST(CPPFLAGS)
AC_SUBST(LDFLAGS)

AC_CONFIG_FILES([Makefile
		data/Makefile
		src/Makefile
		tests/Makefile])
AC_OUTPUT
