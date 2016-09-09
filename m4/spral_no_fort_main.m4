# ifort fails to link if we have a C main without special flags.
# This macro will detect required flags and set NO_FORT_MAIN to contain them.
AC_DEFUN([SPRAL_NO_FORT_MAIN], [
AC_MSG_CHECKING([flags to link C main with $FC])

# Compile to .o using C
AC_LANG_PUSH(C)
AC_LANG_CONFTEST([AC_LANG_SOURCE([ int main(void) { return 0; } ])])
(eval "$ac_compile") 2>conftest.err
ac_status=$?
AS_IF([test $ac_status != 0],[
   echo "$ac_compile" >> config.log
   cat conftest.err >> config.log
   echo "Failed program was" >> config.log
   cat conftest.c >> config.log
   ])
AC_LANG_POP(C)

# Try and link with Fortran
AC_LANG_PUSH(Fortran)
spral_no_fort_main_done=""
spral_link='$FC -o conftest$ac_exeext $FCFLAGS $LDFLAGS $ac_fcflags_srcext conftest.$ac_cv_objext $LIBS >&5'
# No flags
NO_FORT_MAIN="" >> config.log
echo "Trying '$NO_FORT_MAIN'" >> config.log
(eval "$spral_link") 2>conftest.err
ac_status=$?
AS_IF([test $ac_status != 0],[
   echo `echo "$spral_link"` >> config.log
   cat conftest.err >> config.log
   echo "Failed program was" >> config.log
   cat conftest.c >> config.log
   ])
AS_IF([test $ac_status = 0], [
   AC_MSG_RESULT(none);spral_no_fort_main_done="yes"
   ])
# -nofor_main
AS_IF([test "x$spral_no_fort_main_done" == "x"],[
   NO_FORT_MAIN="-nofor_main"
   echo "Trying '$NO_FORT_MAIN'" >> config.log
   save_FCFLAGS=$FCFLAGS
   FCFLAGS="$FCFLAGS $NO_FORT_MAIN"
   (eval "$spral_link") 2>conftest.err
   ac_status=$?
   AS_IF([test $ac_status != 0],[
      echo `echo "$spral_link"` >> config.log
      cat conftest.err >> config.log
      echo "Failed program was" >> config.log
      cat conftest.c >> config.log
      ])
   AS_IF([test $ac_status = 0], [
      AC_MSG_RESULT($NO_FORT_MAIN)
      spral_no_fort_main_done="yes"
      ])
   FCFLAGS=$save_FCFLAGS
   ])
AC_LANG_POP(Fortran)

# Clean up
rm -rf conftest*

# Fail if we didn't succeed
AS_IF([test "x$spral_no_fort_main_done" == "x"],[
      NO_FORT_MAIN=""
      AC_MSG_ERROR("failed to determine how to link")
      ])
AC_SUBST(NO_FORT_MAIN)

])dnl SPRAL_NO_FORT_MAIN
