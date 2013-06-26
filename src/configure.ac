AC_PREREQ([2.67])
AC_INIT([Msh Package], [1.0])

AC_PROG_CXX
AC_LANG(C++)

AC_CHECK_PROG([HAVE_MKOCTFILE], [mkoctfile], [yes], [no])
if [test $HAVE_MKOCTFILE = "no"]; then
  AC_MSG_ERROR([mkoctfile required to install $PACKAGE_NAME])
fi

AC_CHECK_HEADERS([dolfin.h], [ ],
  [AC_MSG_WARN([dolfin required to install $PACKAGE_NAME])]
)
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
