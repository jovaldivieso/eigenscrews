# EigenScrews 🔩

**EigenScrews** is an open source C++ library for performing screw calculus for rigid body kinematics, implemented using the [Eigen](https://gitlab.com/libeigen/eigen) linear algebra library.

It is a direct C++ port of the `Screws.m` Mathematica package developed by Richard Murray and Sudipto Sur at Caltech. The library follows the mathematical framework presented in [1].

## Dependencies

- **C++17** (or later)
- **CMake** 3.25+ (Tested on 4.2.1)
- **Eigen** 3.4+ or 5.x (Tested on 5.0.1)

## Build and Usage

### Compilation

```bash
mkdir build && cd build
cmake ..
make
```

### Running Examples

To verify the library and see basic usage:

```bash
./examples/screws_demo
```

### Integration

To use this library in your own C++ project:

```cpp
#include <eigenscrews/screws.hpp>
```

## Credits and Attribution

This software is a C++ derivation of the original **Screws.m** Mathematica package.

**Original Software**: `Screws.m` (Mathematica)  
**Authors**: Richard Murray and Sudipto Sur  
**Institution**: Division of Engineering and Applied Science, California Institute of Technology  
**Source**: [Caltech CDS - Mathematica Software for Screw Calculus](https://www.cds.caltech.edu/~murray/books/MLS/software.html)

## References
[1] R. M. Murray, Z. Li, and S. S. Sastry, *A Mathematical Introduction to Robotic Manipulation*, CRC Press, 1994.
