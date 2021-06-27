## Tiny PGA
![CMake Test Status](https://github.com/sergehog/tiny_pga/actions/workflows/master.yml/badge.svg?branch=master)

Tiny, header-only, Differentiable and Template-optimized C++ library for Projective Geometric Algebra (PGA)

(at least attempt to make it)

This is my own research and study project, where I'm trying to develop practical tools and algorithms using PGA.

Initialy, focus was in creating my own highly-templated header-only library, however later it looked like it's getting too hard to templetize every element of the multivector.
Thus, I temporally started to use `pga3d.h` implementation from https://bivector.net site.
And now it looks like it works pretty well with my Automatic Differentiation (AD) library, and probably I'll continue using it, instead writing my own PGA.


