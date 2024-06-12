# Include

This is were all of the header files for the project are located. 

While not in this repository, if you installed a local copy of GEOS, it is assumed by the 
makefile that those headers can be found in this directory as well. As far as I know, a 
system wide copy of GEOS will be found automatically if present. 

## Files:

### `geoutil.h` 

This is were all of my geometric related utility functions are. These are small, abstractable helper functions which interact directly with GEOS shapes. The functions defined in this header are listed below:

- `centerGeometry` takes an input `GEOSGeometry*` and returns a pointer to a new geometry which has been translated to that it's origin is at (0, 0). This is necessary so that distinct polygons say on a map will overlap and thus have a Jaccard Distance of less than one. 
- `createGrid` takes an input `GEOSGeometry*`, checks that it is a rectangle and a than creates a grid of smaller rectangles over that one with an adjustable (via parameter) number of cells per size.
- `fillRatio` takes two polygons, a base and an overlay shape and returns the percentage of the base shape which is "filled" by the overlay shape.
- `getCenter` Returns the centroid of the input `GEOSGeometry*` in a easier to work with format, namely two doubles, rather than as a coordinate sequence which you have to parse through.
- `jaccardDistance` computes the Jaccard Distance between two polygons and returns it as a double.
- `LSHHash` Computes a locality sensitive hash for a shape of a set length using the algorithm outlined in “Simple and Efficient Weighted Minwise Hashing” by A. Shrivastava.
- `minimumBoundingRectangle` takes a list of polygons and returns the minimum bounding rectangle for the entire set of polygons.

### `parse_geodata.h`

Here are all of my tools to read in data to the GEOS library. 

- `read_csv` reads the data format described in the README of the `data/` directory under `csv`. Not used in the final project implementation. 
- `read_wkt` reads the data format described in the README of the `data/` directory under `wkt`. This is the primary format used for input to the project. It handles all of the file reading and returns a vector of read geometries already filtered to be only polygons.
- `read_wkt_parallel` functions very similarly to `read_wkt` but uses multi-threading to accelerate the parsing of the WKT. Also not used in the final project in favor of more complex input threads, but is used in `brute_force.cpp` as a proof of concept. 

### `thread.h`

There is the body functions for all three types of threads that this project supports. Since I 
don't think that I can get the direct output from thread functions using the C++ standard library, 
each thread has all of its inputs wrapped into a `struct` and all outputs are likewise placed into 
a `struct` which is created and otherwise managed by the main thread. 

1. Input Threads. These read in data from files which have been split by the main thread (see `util.h`) and perform all the needed preprocessing on the polygons, such as centering them and tracking a local bounding box which is later used to create a global bounding box on the area of interest. 
2. Construction Threads. These take the processed output of the input threads and first construct sketches of each of the polygons. Once all of the sketches are generated, the threads wait at a barrier for the `gridFilterer` class to process all of the cell usage data and remove cells from the grid which fail to meet some usage requirements. This is done to reduce the dimensionality of the data at the cost of reduced accuracy. After that, the construction threads resume and revise the sketches to drop the same cells. Then locality sensitive hashes are computed and the polygons are inserted into the hash maps. 
3. Query Threads. These threads read from the query file, then follow a process similar to the construction threads. The polygons are sketched and locality sensitive hashed. However, instead of being inserted into the hash maps, we fetch the bucket from each hash table matching the hash of the polygon and perform a linear search for the nearest neighbors of that polygon, which forms the output of this thread. 

### `util.h`

A collection of utility functions which do not deal with GEOS shapes directly. These are mostly 
related to input and output operations like specialized print functions. This is also were my 
`hashMap` class is defined. I choose to define it here, since it uses a template type to be a more
general hash map which I happen to use for GEOS shapes. 

#### I/O Related

- `aprintf` stands for atomic `printf` and takes a `printf` style format string and other parameters which it uses to write into a temporary buffer before using the `write` system call to get atomic output. It does require knowing the number of characters to be printed. 
- `geosErrorHandler` writes GEOS errors out to standard error in an atomic manner. It assumes that these error messages are no longer than 256 characters long. If they are longer, the output is truncated.
- `geosMessageHandler` writes GEOS messages to standard out in an atomic manner. Like `geosErrorHandler`, it assumes these messages are not more than 256 characters and truncates if they are.
- `integerLength` returns the number of characters that an integer value would required to be printed in full. Very helpful when used with `aprintf`.
- `splitFile` uses the unix command `split` to split a long file into a number of smaller files while respective newlines. It does this in a round-robin style which is very useful for breaking up data files for this project. 

#### Hash Map

The `HashMap` class is a one-to-many hash map which uses a vector of integers as a key and a template type as the values. 
Since it is one-to-many, one key can have any number of values mapped to it. 
This map is also special since it offers stronger concurrency than the C++ standard library `unordered_map` by placing a mutex on each bucket. 
Thus, multiple threads can be inserting into the map so long as these insertions are in different buckets. 

Reading from the hash map can happen whenever without restriction. 
While there is some dangerous things that *could* happen here (basically if some is inserted into the bucket while we're reading from it), they aren't a factor in this case all of the insertions happen before any of the queries for this project. 

Be warned that there is no way to remove an element from the map after it has been inserted without deleting the whole map, since that functionality is not needed for this project. 
