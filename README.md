# Numerical Analysis

Numerical Analysis contains some basic numerical methods for solving equations, integrating functions and parameterizing curves with Bezier Polynomials.

## Installation

* Copy ``` include ``` and ``` libxx ``` folders located in ``` Release ``` folder in a directory ``` dependencies\numerical-analysis ``` inside your solution directory.
* Add directory ``` $(SolutionDir)dependencies\numerical-analysis\include ``` in the field ``` Project Properties -> C/C++ -> General -> Additional Include Directories ```.
* Add libraries ``` cpp-extended.lib;math-extended.lib;polynomial.lib;numerical-analysis.lib ``` in the field ``` Project Properties -> Linker -> Input -> Additional Dependencies ```.
* You are ready to use Numerical Analysis functions.
