AC_INIT([gmmem], [0.1], [pikalaxalt@gmail.com])
AM_INIT_AUTOMAKE
AC_PROG_CC
AC_CYGWIN
AC_MINGW32
AC_CHECK_TYPES([int, unsigned int, double])
AC_C_INLINE
AC_DEFINE([HAVE_INLINE], [1])
AC_CHECK_LIB([blas], [dger], [], [AC_MSG_ERROR([Please install BLAS])])
AC_CHECK_LIB([gsl], [gsl_blas_dger], [], [AC_MSG_ERROR([Please install GSL with BLAS support])])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
