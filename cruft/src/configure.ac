AC_PREREQ([2.67])
AC_INIT([Msh Package], [1.0])

AC_PROG_CXX
AC_LANG(C++)

AC_CHECK_PROG([HAVE_MKOCTFILE], [mkoctfile], [yes], [no])
if [test $HAVE_MKOCTFILE = "no"]; then
  AC_MSG_ERROR([mkoctfile required to install $PACKAGE_NAME])
fi

AC_CHECK_HEADER([dolfin.h], 
  [AC_SUBST(ac_dolfin_cpp_flags,-DHAVE_DOLFIN_H) AC_SUBST(ac_dolfin_ld_flags,-ldolfin)],
  [AC_MSG_WARN([dolfin headers could not be found, som functionalities will be disabled, don't worry your package will still be working, though.])]
 )

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
