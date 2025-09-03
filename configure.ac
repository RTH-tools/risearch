AC_INIT([RIsearchUnified], [1.0], [you@example.com])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_PROG_CC
AC_ARG_WITH([matpath],
  [AS_HELP_STRING([--with-matpath=PATH], [Path to RIsearch tables (default: /usr/share/risearch/tables)])],
  [MATPATH="$with_matpath"],
  [MATPATH="/usr/share/risearch/tables"]
)

AC_DEFINE_UNQUOTED([MATPATH], ["$MATPATH"], [Path to RIsearch tables])
AC_PROG_CC
AC_OPENMP
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile RIsearch1/src/Makefile RIsearch2/src/Makefile])
AC_OUTPUT
