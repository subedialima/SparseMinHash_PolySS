# Approximate Nearest Neighbor with Real World Geometry

This repository contains the source code for my work this summer as part of the Data Science Across Disciplines REU at Marquette University. 
More information on that program can be found [here](https://reu.cs.mu.edu/index.php/Main_Page) with project specific information available on both the [project description page](https://reu.cs.mu.edu/index.php/Approximate_Nearest_Neighbor_Searches_of_Real_World_Geometry) and my [personal wiki page](https://reu.cs.mu.edu/index.php/User:Schwennesen).

This research hopes to use locality sensitive hashing to find the nearest neighbors of polygons using Jaccard Distance as our metric. 
The project is primarily written in C++ using the reentrant C API for the [GEOS library](https://libgeos.org/) for the hashing algorithms.
Data visualization is done with Python using the [Shapely](https://shapely.readthedocs.io/en/stable/manual.html) along with matplotlib.

The project currently has multiple directories.

1. `bin` Where the output of the make files are placed. 
2. `data` Where the input data is stored and the output data is written. Due to the large nature of the datasets, complete sets are not available in this repository but can be found at the [UCR Star](https://star.cs.ucr.edu/) website. A few small subsets of their data may be present in the repository.
3. `example` When I started to configure the relevant libraries like GEOS I used these small, simple programs to build my understanding of the API and ensure that the library was configured correctly on the system. 
4. `include` The header files for the gnuplot connector and version 3.10.3 of GEOS.
5. `lib` The shared libraries for GEOS.
6. `src` The actual source code for the nearest neighbor search code, written in C++11. 
7. `util` A collection of small scripts which I found helpful over the course of the project. 

The `data`, `example`, `include`, `src` and `util` directories all have their own README for more information. 

Within the code itself, each function, class, type, etc in a header file has a description of its respective parameters, return values, members, etc, which I hope will make understanding the code easier. 
Also, please remember that I learned C first and as such am more comfortable working with pointers than standard library types, although I have tried to use standard library data structures and C++ practices whenever possible. I apologize if I missed something.

To get started, build the project with `make all` in this directory.
