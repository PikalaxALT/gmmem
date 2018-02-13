#!/bin/sh

set -x
./envclean.sh
aclocal && autoconf && automake --add-missing
