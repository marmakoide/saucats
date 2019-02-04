[![Build Status](https://travis-ci.com/marmakoide/saucats.svg?branch=master)](https://travis-ci.com/marmakoide/saucats) [![MIT License badge](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/marmakoide/saucats/blob/master/LICENSE)

# saucats
A C++ library to define and manipulate geometrical shapes using distance fields

This a header-only library (C++14 dialect) and depends on the [Eigen library](http://eigen.tuxfamily.org).

## Getting Started

### Prerequisites

* *saucats* depends on [Eigen library](http://eigen.tuxfamily.org).
* The build system assumes you have Python 2.7 or later installed. 
* *saucats* compiles with GCC 6 on Linux.

### Installing

*saucats* is an header-only library, therefore, installing it can be simply
copying the header files wherever you want them. However, you might want to
build the test cases and run them, which requires a compilation step.


The build system, waf, is part of the source tree. You first run the configuration
step, that will check that your environment fullfill all the requirements to run
*saucats* and build the test cases.
```
./waf configure
```


Then, you can build the debug release, to get the test cases
```
./waf build_debug
```

If you do not want to compile the test cases
```
./waf build_release
```

You can then install the headers and the test cases
```
./waf install_debug
```

Or you can install just the headers
```
./waf install_release
```

## Authors

* **Alexandre Devert** - *Initial work* - [marmakoide](https://github.com/marmakoide)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details


