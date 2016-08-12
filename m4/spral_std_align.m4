# Some versions of libstdc++ fail to defined std::align.
# This macro sets HAVE_STD_ALIGN if it is defined.
AC_DEFUN([SPRAL_STD_ALIGN], [
AC_MSG_CHECKING(for std::align)
AC_REQUIRE([AC_PROG_CXX])
AC_LANG_PUSH(C++)

AC_TRY_LINK([#include <memory>],
[std::size_t alignment, size, space; void* ptr;
std::align(alignment, size, ptr, space);],
AC_MSG_RESULT(yes);AC_DEFINE(HAVE_STD_ALIGN, 1, [Define to 1 if you have std::align().]),
AC_MSG_RESULT(no)
)

AC_LANG_POP(C++)
])dnl SPRAL_STD_ALIGN
