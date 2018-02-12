#!/bin/sh

set -x
rm -f aclocal.m4 compile config.guess config.log config.status config.sub configure \
	depcomp install-sh Makefile.in missing INSTALL COPYING
aclocal && autoconf && automake --add-missing
