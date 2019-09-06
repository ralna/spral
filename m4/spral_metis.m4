# Note: loosely based on immxdb_lib_metis.m4 and acx_blas.m4
AC_DEFUN([SPRAL_METIS], [
AC_MSG_CHECKING(for METIS library)
AC_REQUIRE([AC_PROG_CC])
#
# User hints...
#
AC_ARG_WITH([metis],
   [AC_HELP_STRING([--with-metis=<lib>],
   [user METIS library <lib>])])
case $with_metis in
   -* | */* | *.a | *.so | *.so.* | *.o) METIS_LIBS="$with_metis" ;;
    *) METIS_LIBS="-L$with_metis -lmetis" ;;
esac

# Get fortran linker names for function of interest
AC_F77_FUNC(metis_nodend)

spral_metis_ok=no

# Check supplied location
if test $spral_metis_ok = no; then
if test "x$METIS_LIBS" != x; then
   save_LIBS="$LIBS"; LIBS="$METIS_LIBS $LIBS -lm"
   AC_MSG_CHECKING([for $metis_nodend in $METIS_LIBS])
   AC_TRY_LINK_FUNC($metis_nodend, [spral_metis_ok=yes], [METIS_LIBS=""])
   AC_MSG_RESULT($spral_metis_ok)
   LIBS="$save_LIBS"
fi
fi

AC_SUBST(METIS_LIBS)

# Try just -lmetis
if test $spral_metis_ok = no; then
   AC_CHECK_LIB(metis, $metis_nodend, [spral_metis_ok=yes; METIS_LIBS="-lmetis"], [], [-lm])
fi

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$spral_metis_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_METIS,1,[Define if you have a MeTiS library.]),[$1])
        :
else
        spral_metis_ok
        $2
fi

# Determine which metis interface to compile
save_LIBS="$LIBS"; LIBS="$METIS_LIBS $LIBS -lm"
AC_MSG_CHECKING([version of METIS])
AC_TRY_LINK_FUNC([METIS_Free],
                 [
                  METIS_VERSION="5"
                  AC_MSG_RESULT("version 5")
                  ],
                 [
                  METIS_VERSION="4"
                  AC_MSG_RESULT("version 4")
                  ])
AC_SUBST(METIS_VERSION)
LIBS="$save_LIBS"

AC_DEFINE(SPRAL_HAVE_METIS_H, [0])

# Check for metis include
AC_ARG_WITH([metis_inc_dir],
        [AC_HELP_STRING([--with-metis-inc-dir=<inc-dir>],
        [user METIS include directory <inc-dir>])],
        [
                metis_inc_dir="$withval"
        ],
        [metis_inc_dir=no])

if test "$metis_inc_dir" != "no" ; then
   CFLAGS="-I$metis_inc_dir $CFLAGS"
   FFLAGS="-I$metis_inc_dir $FFLAGS"
   AC_CHECK_HEADERS([metis.h], [AC_DEFINE([SPRAL_HAVE_METIS_H], [1], [Define to 1 if you have metis.h.])])

   AC_CHECK_SIZEOF([idx_t], [], [[#include <metis.h>]])
fi

])dnl SPRAL_METIS

