# Note: loosely based on immxdb_lib_metis.m4 and acx_blas.m4
AC_DEFUN([SPRAL_GTG], [
AC_MSG_CHECKING(for GTG library)
AC_REQUIRE([AC_PROG_CC])
#
# User hints...
#
AC_ARG_WITH([gtg],
   [AC_HELP_STRING([--with-gtg=<lib>],
   [use GTG library <lib>])])
case $with_gtg in
   -* | */* | *.a | *.so | *.so.* | *.o) GTG_LIBS="$with_gtg" ;;
    *) GTG_LIBS="-L$with_gtg -lgtg" ;;
esac

spral_gtg_lib_ok=no
spral_gtg_inc_ok=no

# Check supplied location
if test $spral_gtg_lib_ok = no; then
if test "x$GTG_LIBS" != x; then
   save_LIBS="$LIBS"; LIBS="$GTG_LIBS $LIBS -lm"
   GTG_INCLUDE="$with_gtg"
   AC_MSG_CHECKING([for initTrace in $GTG_LIBS])
   AC_TRY_LINK_FUNC([initTrace], [spral_gtg_lib_ok=yes], [GTG_LIBS="";GTG_INCLUDE=""])
   AC_MSG_RESULT($spral_gtg_lib_ok)
   LIBS="$save_LIBS"
fi
fi

# Check supplied location if its gtg src directory
if test $spral_gtg_lib_ok = no; then
if test "x$with_gtg" != x; then
   GTG_LIBS="-L$with_gtg/src/.libs -lgtg"
   GTG_INCLUDE="-I$with_gtg/inc"
   save_LIBS="$LIBS"; LIBS="$GTG_LIBS $LIBS -lm"
   AC_MSG_CHECKING([for initTrace in $GTG_LIBS])
   AC_TRY_LINK_FUNC([initTrace], [spral_gtg_lib_ok=yes], [GTG_LIBS="";GTG_INCLUDE=""])
   AC_MSG_RESULT($spral_gtg_lib_ok)
   LIBS="$save_LIBS"
fi
fi

# Try just -lgtg
if test $spral_gtg_lib_ok = no; then
   AC_CHECK_LIB(gtg, initTrace, [spral_gtg_lib_ok=yes; GTG_LIBS="-lgtg";GTG_INCLUDE=""], [], [-lm])
fi

save_CPPFLAGS="$CPPFLAGS"; CPPFLAGS="$CPPFLAGS $GTG_INCLUDE"
AC_CHECK_HEADER([GTG.h], [spral_gtg_inc_ok=yes], [GTG_INCLUDE=""])
CPPFLAGS="$save_CPPFLAGS"

AC_SUBST(GTG_INCLUDE)
AC_SUBST(GTG_LIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$spral_gtg_lib_ok$spral_gtg_inc_ok" = xyesyes; then
        ifelse([$1],,AC_DEFINE(HAVE_GTG,1,[Define if you have a GTG library.]),[$1])
        :
else
        $2
fi

])dnl SPRAL_GTG
