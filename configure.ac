AC_INIT([gmmem], [0.1], [https://github.com/pikalaxalt/gmmem/issues])
AM_INIT_AUTOMAKE
AC_PROG_CC
AC_CYGWIN
AC_MINGW32
AC_CHECK_TYPES([int, unsigned int, double])
AC_C_INLINE
if test "$ac_cv_c_inline" != no ; then
  AC_DEFINE(HAVE_INLINE,1)
  AC_SUBST(HAVE_INLINE)
fi
AX_BLAS([AC_DEFINE([HAVE_BLAS])], [AC_MSG_ERROR([Please install BLAS or OpenBLAS])])
AC_CHECK_LIB([m],[cos], [], [AC_MSG_ERROR([Please install GSL with BLAS support])])
AC_CHECK_LIB([gslcblas],[cblas_dgemm], [], [AC_MSG_ERROR([Please install GSL with BLAS support])])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm], [], [AC_MSG_ERROR([Please install GSL with BLAS support])])
AX_PATH_GSL(2.4, [], [AC_MSG_ERROR([GSL found but no multivariate Gaussian!])])
AC_ARG_ENABLE([debug], [AS_HELP_STRING([--enable-debug], [Enable debug prints])], [AC_DEFINE([DEBUG])])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
