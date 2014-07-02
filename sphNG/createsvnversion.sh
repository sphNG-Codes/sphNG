#!/bin/sh

echo '      CHARACTER*20 svnversion' > svn_version.h

# Check whether subversion is available
which svnversion 2>&1 > /dev/null

# Set the version number using svnversion if available or set to unknown otherwise
if [ $? -eq 0 ] ; then
    echo '      PARAMETER (svnversion = "-'`svnversion`'")' >> svn_version.h
else
    echo '      PARAMETER (svnversion = "-'unknown'")' >> svn_version.h
fi

