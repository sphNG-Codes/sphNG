#!/bin/csh

# Check whether subversion is available
which svnversion > /dev/null

# Set the version number using svnversion if available or set to unknown otherwise
echo '      CHARACTER*20 svnversion' >! svn_version.h
if ($? == 0) then 
    echo '      PARAMETER (svnversion = "-'`svnversion`'")' >> svn_version.h
else
    echo '      PARAMETER (svnversion = "-'unknown'")' >> svn_version.h
endif
