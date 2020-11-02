# Numerical Analysis

Numerical Analysis contains some basic numerical methods for solving equations, integrating functions and parameterizing curves with Bezier Polynomials.

## Installation

### Visual Studio Configuration

* Copy ``` include ``` and ``` libxx ``` folders located in ``` Release ``` folder of the repository and paste in a directory ``` dependencies\numerical-analysis ``` inside your solution directory.
* Open your solution in Visual Studio. In ``` Solution Explorer ``` right click your Project and select ``` Properties ```. Select ``` Configuration ``` and ``` Platform ``` of your Project.
* Add directory ``` $(SolutionDir)dependencies\numerical-analysis\include ``` in the field ``` Project Properties -> C/C++ -> General -> Additional Include Directories ```.
* Add directory ``` $(SolutionDir)dependencies\numerical-analysis\libxx ``` in the field ``` Project Properties -> Linker -> General -> Additional Library Directories ```.
* Add libraries ``` cpp-extended.lib;polynomial.lib;numerical-analysis.lib; ``` in the field ``` Project Properties -> Linker -> Input -> Additional Dependencies ```.
* You are ready to use Numerical Analysis functions!
