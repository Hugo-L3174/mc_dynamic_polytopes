# mc_dynamic_polytopes

This library intends to wrap external polytope computations libraries to compute a DCM dynamic balance region for [mc_rtc](https://jrl.cnrs.fr/mc_rtc/index.html), intending to be used in a controller and providing GUI display of the regions.
For now the computations are wrapped only around the [politopix library](https://www.i2m.u-bordeaux.fr/politopix)

## Dependencies

- [mc_rtc](https://github.com/jrl-umi3218/mc_rtc)
- the github port of [politopix](https://github.com/Hugo-L3174/politopix)
- [qhull](http://www.qhull.org/)

## Installation

After installing the dependencies,

```sh
git clone git@github.com:Hugo-L3174/mc_dynamic_polytopes.git
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=RelWithDebInfo
make -j8
sudo make install
```

## Usage and exemple

In your `mc_rtc` controller's `CMakeLists.txt`, link the library `mc_dynamic_polytopes::DynamicPolytope`.

[Here](https://github.com/Hugo-L3174/polytopeController) is an example of a simple controller using this library.
