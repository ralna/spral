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

spral_gtg_ok=no

# Check supplied location
if test $spral_gtg_ok = no; then
if test "x$GTG_LIBS" != x; then
   save_LIBS="$LIBS"; LIBS="$GTG_LIBS $LIBS -lm"
   AC_MSG_CHECKING([for initTrace in $GTG_LIBS])
   AC_TRY_LINK_FUNC([initTrace], [spral_gtg_ok=yes], [GTG_LIBS=""])
   AC_MSG_RESULT($spral_gtg_ok)
   LIBS="$save_LIBS"
fi
fi

AC_SUBST(GTG_LIBS)

# Try just -lgtg
if test $spral_gtg_ok = no; then
   AC_CHECK_LIB(gtg, initTrace, [spral_gtg_ok=yes; GTG_LIBS="-lgtg"], [], [-lm])
fi

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$spral_gtg_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_GTG,1,[Define if you have a GTG library.]),[$1])
        :
else
        spral_gtg_ok
        $2
fi

])dnl SPRAL_GTG
