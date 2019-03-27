
# A simple metaheuristics framework to couple optimization and simulation models

master: [![master build status](https://ci.appveyor.com/api/projects/status/9h5snds4ymuv4ynq/branch/master?svg=true)](https://ci.appveyor.com/project/jmp75/wila/branch/master) testing: [![testing build status](https://ci.appveyor.com/api/projects/status/9h5snds4ymuv4ynq/branch/testing?svg=true)](https://ci.appveyor.com/project/jmp75/wila/branch/devel)

## Purpose

This framework is designed to incorporate the state of the art in metaheuristics software frameworks, yet limiting the software complexity to users who are interested in applying it without advanced knowledge of software or optimisation research. It has been used mainly to calibrate environmental models, mostly hydrology models.

The purpose of this framework is to define a set of programming interface, rather than replicate optimisation algorithms found in other optimisation frameworks.

## License

This software is released under the LGPL v3. See LICENSE.

## Requirements

This is written in C++ using some features from the C++11 standard. It is known to compile with:

* GCC version 4.8 or above
* Visual C++ as from Visual Studio 2015, or above.

## Getting started

This framework is (almost) header-only, so you do not build a shared binary library.

### Setting dependencies

You will need the following dependencies available:

#### Windows

* [Boost](http://www.boost.org/). You may choose to set up boost on your machine as you wish and sort out the compilation/path settings as you wish. The following is but one approach, is prescriptive, albeit one based on experience to limit some annoyances.
 * Download the Boost installers from e.g. [1.61](https://sourceforge.net/projects/boost/files/boost-binaries/1.61.0) at the time of writing. You will find files such as boost_1_61_0-msvc-12.0-64.exe (or more recent than "12.0" depending on the version of visual studio you use. Download for 64 and/or 32 bits.
 * Install at least one of the boost binary distro, for instance (and most likely) the 64 bit one (boost_1_61_0-msvc-12.0-64.exe). Let's assume for the rest of this document that you install to F:\local\boost_1_61_0\. 
 * We will create a subset of Boost in a folder F:\local\boost.  The following batch script creates it with enough such that the 'wila' unit tests can be compiled.

```
set BOOST_VERSION=1_61
set BOOST_DIR=F:\local\boost_%BOOST_VERSION%_0\
set LIB_BOOST=%BOOST_DIR%lib64-msvc-12.0\
set MYBOOST=F:\local\boost\
set MYBOOST_HDR=%MYBOOST%boost\
set MYBOOST_LIB=%MYBOOST%lib\64\ 
if not exist %MYBOOST_HDR% mkdir %MYBOOST_HDR%
if not exist %MYBOOST_LIB% mkdir %MYBOOST_LIB%

set MYCP=xcopy

set COPYOPTIONS=/Y /R /D

:: Copy the header files:
robocopy %BOOST_DIR%boost\    %MYBOOST_HDR%   /MIR

:: Copy necessary binaries:
%MYCP% %LIB_BOOST%boost_thread-vc120-mt-%BOOST_VERSION%.dll          %MYBOOST_LIB%   %COPYOPTIONS%
%MYCP% %LIB_BOOST%boost_thread-vc120-mt-gd-%BOOST_VERSION%.dll       %MYBOOST_LIB%   %COPYOPTIONS%
%MYCP% %LIB_BOOST%libboost_thread-vc120-mt-%BOOST_VERSION%.lib       %MYBOOST_LIB%   %COPYOPTIONS%
%MYCP% %LIB_BOOST%libboost_thread-vc120-mt-gd-%BOOST_VERSION%.lib    %MYBOOST_LIB%   %COPYOPTIONS%
                                                                    
%MYCP% %LIB_BOOST%boost_system-vc120-mt-%BOOST_VERSION%.lib          %MYBOOST_LIB%   %COPYOPTIONS%
%MYCP% %LIB_BOOST%boost_system-vc120-mt-gd-%BOOST_VERSION%.dll       %MYBOOST_LIB%   %COPYOPTIONS%
%MYCP% %LIB_BOOST%libboost_system-vc120-mt-%BOOST_VERSION%.lib       %MYBOOST_LIB%   %COPYOPTIONS%
%MYCP% %LIB_BOOST%libboost_system-vc120-mt-gd-%BOOST_VERSION%.lib    %MYBOOST_LIB%   %COPYOPTIONS%
                                                                    
%MYCP% %LIB_BOOST%boost_date_time-vc120-mt-%BOOST_VERSION%.dll       %MYBOOST_LIB%   %COPYOPTIONS%
%MYCP% %LIB_BOOST%boost_date_time-vc120-mt-gd-%BOOST_VERSION%.dll    %MYBOOST_LIB%   %COPYOPTIONS%
%MYCP% %LIB_BOOST%libboost_date_time-vc120-mt-%BOOST_VERSION%.lib    %MYBOOST_LIB%   %COPYOPTIONS%
%MYCP% %LIB_BOOST%libboost_date_time-vc120-mt-gd-%BOOST_VERSION%.lib %MYBOOST_LIB%   %COPYOPTIONS%
                                                                    
%MYCP% %LIB_BOOST%boost_chrono-vc120-mt-%BOOST_VERSION%.dll          %MYBOOST_LIB%   %COPYOPTIONS%
%MYCP% %LIB_BOOST%boost_chrono-vc120-mt-gd-%BOOST_VERSION%.dll       %MYBOOST_LIB%   %COPYOPTIONS%
%MYCP% %LIB_BOOST%libboost_chrono-vc120-mt-%BOOST_VERSION%.lib       %MYBOOST_LIB%   %COPYOPTIONS%
%MYCP% %LIB_BOOST%libboost_chrono-vc120-mt-gd-%BOOST_VERSION%.lib    %MYBOOST_LIB%   %COPYOPTIONS%
```

* You need to install Threadpool (http://threadpool.sourceforge.net/), You can fetch the source code from [this fork on github](https://github.com/jmp75/threadpool/tree/master) using git.

You will find in the top directory of wila the file [./wila.props.in](./wila.props.in). Copy it to a file named wila.props under "My Documents", and adjust the values BoostPath and BoostThreadPool to the relevant locations you chose.

```xml
  <PropertyGroup Label="UserMacros">
    <BoostPath>F:/local/boost/</BoostPath>
    <BoostThreadPool>F:/src/github_jm/threadpool/</BoostThreadPool>
  </PropertyGroup>
```

Open "wila\tests\mhcpptest.vcxproj". Make sure you are in the appropriate configuration (Debug/x64) and it should compile.

#### Linux

* [Boost](http://www.boost.org/)
* Threading Building Blocks (https://www.threadingbuildingblocks.org/). 
* Threadpool (http://threadpool.sourceforge.net/, You can fetch the source code from [this fork on github](https://github.com/jmp75/threadpool/tree/master) and place in the top of your source tree for CMake to find it)

```sh
cmake .
make package
```

Alternately build as a shared library:

```sh
cmake -DAS_SHARED_LIB=True .
make
make install
```

### Documentation

You can find some introductory sample [in the doc folder, Getting Started](./doc/GettingStarted.md)

### Unit tests

The step above will also try and build the test application. To run the tests:

```sh
./wila_tests
```

### Further Documentation

TODO will probably set up a github page

## What's the name of the repo about

This is from Slavic mythology, as the description was a metaphor somewhat appropriate for the field of optimisation.

"In Slavic mythology, there is a form of nymph which lies somewhere between a ghost and a fairy. \[...\] They can either blend into the wind as incorporeal shapes — translucent and intangible — or they can become solid, touching, and being touched, by the natural world around them. \[...\] If such an easily enticed man were to go searching for a Wila, he would most likely find her in places similar to those which the fairies and nymphs prefer—on hill tops or mounds, or in the center of a ring of trees."

* Ryan Stone, [Beware the Wandering Wilas](http://www.ancient-origins.net/myths-legends-europe/beware-wandering-wilas-002273)
* [Supernatural beings in Slavic folklore - Vila](https://en.wikipedia.org/wiki/Supernatural_beings_in_Slavic_folklore#Vila) on wikipedia.org
