A simple metaheuristics framework to couple optimization and simulation models
==============================================

# Purpose

This framework is designed to incorporate the state of the art in metaheuristics software frameworks, yet limiting the software complexity to users who are interested in applying it without advanced knowledge of software or optimisation research. It has been used mainly to calibrate environmental models, mostly hydrology models.

The purpose of this framework is to define a set of programming interface, rather than replicate optimisation algorithms found in other optimisation frameworks.

# License

This software is released under the LGPL v3. See LICENSE.

# Requirements

This is written in C++ using some features from the C++11 standard. It is known to compile with:

* GCC version 4.8 or above
* Visual C++ as from Visual Studio 2013 update 5, or above.

# Getting started

## Build library and packages

You will need the following dependencies available:

* Boost
* Threading Building Blocks (https://www.threadingbuildingblocks.org/)
* threadpool (http://threadpool.sourceforge.net/, download and place in the top
  of your source tree for CMake to find it)

```sh
cmake .
make package
```

## Unit tests

The step above will also try and build the test application. To run the tests:

```sh
./wila_tests
```

## Further Documentation

TODO will probably set up a github page

# What's the name of the repo about?

"In Slavic mythology, there is a form of nymph which lies somewhere between a ghost and a fairy. \[...\] They can either blend into the wind as incorporeal shapes — translucent and intangible — or they can become solid, touching, and being touched, by the natural world around them. \[...\] If such an easily enticed man were to go searching for a Wila, he would most likely find her in places similar to those which the fairies and nymphs prefer—on hill tops or mounds, or in the center of a ring of trees."

* Ryan Stone, [Beware the Wandering Wilas](http://www.ancient-origins.net/myths-legends-europe/beware-wandering-wilas-002273)
* [Supernatural beings in Slavic folklore - Vila](https://en.wikipedia.org/wiki/Supernatural_beings_in_Slavic_folklore#Vila) on wikipedia.org

