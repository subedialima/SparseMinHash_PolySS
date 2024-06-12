// Implementation of the two query methods: LSH and brute force

#include "query.h"
#include <mpi.h>



void queryThreadFunct(
                     unordered_map<string, GEOSGeometry*> centeredQueryGeosWithID,
                     unordered_map<string, pair<GEOSGeometry*, vector<double>>> geometrySketchMap,
                     vector<HashMap<std::string>*> hashMaps,
                     int numNeighbors, 
                     int hashLength, 
                     vector<unsigned int>& seeds, 
                     int rank,
                    // std::vector<std::pair<std::pair<std::string, GEOSGeom_t*>, std::vector<std::pair<double, std::pair<std::string, GEOSGeom_t*>>>>>& LSHoutput,
                     ResultType& LSHoutput,
                     vector<vector<GEOSGeometry*>>& global_MBR_grid,
                      const vector<bool>& cellsToKeep) {

    // Create a new GEOS context handler
    GEOSContextHandle_t ctx = GEOS_init_r();
    GEOSContext_setNoticeHandler_r(ctx, geosMessageHandler);
    GEOSContext_setErrorHandler_r(ctx, geosErrorHandler);
    
   
    // Generate sketches for Query Geometries
   

unordered_map<string, pair<GEOSGeometry*, vector<double>>> QuerygeometrySketchMap;

// // Start Hashing time for polygon 
  double startTimeQuerysketch = MPI_Wtime();
// // 
// // // // Create sketch for each geometry in the unordered_map centeredQueryGeosWithID
// for (const auto& pair : centeredQueryGeosWithID) {
//     const string& id = pair.first; // Unique ID of the geometry
//     GEOSGeometry* geometry = pair.second; // Geometry pointer
// 
//     vector<double> sketchResult = sketch(ctx, &global_MBR_grid, geometry);
//     if (sketchResult.empty()) {
//         cerr << "Warning: Sketch generated for geometry with ID " << id << " is empty." << endl;
//     } else {
//         QuerygeometrySketchMap[id] = make_pair(geometry, sketchResult);
//     }
//   //  std::cout<< " - query Sketch Size: " << sketchResult.size() << std::endl;
//     
// }
// 
// // OVERALL time for querying 
//     double endTimeQuerySketch = MPI_Wtime();
//     double totalTimeQuerySketch = endTimeQuerySketch - startTimeQuerysketch;
//     
//  double maxTotalTimeQuerySketch;
// 
// // Use MPI_Reduce to find the maximum hashing time across all processes
// MPI_Reduce(&totalTimeQuerySketch, &maxTotalTimeQuerySketch, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
// 
// // The master process (rank 0) prints the maximum hashing time
// if (rank == 0) 
// {
//     std::cout << "Maximum time for query sketch creation across all processes: " << maxTotalTimeQuerySketch << " seconds" << std::endl;
// }
 


// Create and trim sketch for each geometry in the unordered_map centeredQueryGeosWithID
for (const auto& pair : centeredQueryGeosWithID) {
    const string& id = pair.first; // Unique ID of the geometry
    GEOSGeometry* geometry = pair.second; // Geometry pointer

    vector<double> fullSketch = sketch(ctx, &global_MBR_grid, geometry); // Generate the full sketch first
    vector<double> trimmedSketch;
    for (size_t i = 0; i < fullSketch.size(); ++i) {
        if (cellsToKeep[i]) {
            trimmedSketch.push_back(fullSketch[i]); 
            // Keep only the cells marked by cellsToKeep
        }
         
    }
   // std::cout << "ID: " << id << " - Trimmed Sketch Size: " << trimmedSketch.size() << std::endl;
    if (trimmedSketch.empty())
     {
        cerr << "Warning: Trimmed sketch generated for query geometry with ID " << id << " is empty." << endl;
    } else 
    {
        QuerygeometrySketchMap[id] = make_pair(geometry, trimmedSketch);
    }
}


// OVERALL time for querying 
    double endTimeQuerySketch = MPI_Wtime();
    double totalTimeQuerySketch = endTimeQuerySketch - startTimeQuerysketch;
    
 double maxTotalTimeQuerySketch;

// Use MPI_Reduce to find the maximum hashing time across all processes
MPI_Reduce(&totalTimeQuerySketch, &maxTotalTimeQuerySketch, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

// The master process (rank 0) prints the maximum hashing time
if (rank == 0) 
{
    std::cout << "Maximum time for query sketch creation across all processes: " << maxTotalTimeQuerySketch << " seconds" << std::endl;
}
 
// //printing sketch along with its id
//  for (const auto& item : QuerygeometrySketchMap) {
//      const string& id = item.first; // Unique ID
//      const auto& geometryAndSketch = item.second; // Pair of geometry pointer and sketch vector
//       const vector<double>& sketchResult = geometryAndSketch.second; // Extracting the sketch vector
//  // 
// //       cout << "ID: " << id << " - Sketch: ";
// //     for (double value : sketchResult) {
// //         cout << value << " ";
//     //}
//     cout << endl; // End of a single sketch
//  }
//  
//Similar process as before for generating sketches for query geometries

// Start Hashing time for polygon 
  double startTimeQuery = MPI_Wtime();
  
double hashingStart, hashingEnd, hashingTotal = 0.0;
double queryingStart, queryingEnd, queryingTotal = 0.0;
double queryingSortEnd, queryingSortTotal =0.0;
double hashmap_key_search_start, hashmap_key_search_end, hashmap_key_search_total = 0.0;


for (const auto& item : QuerygeometrySketchMap)
 {
 const string& queryID = item.first;
  GEOSGeometry* queryGeo = item.second.first;
 const vector<double>& sketch = item.second.second;
// const vector<double>& querySketch = item.second.second; 

 //vector<pair<double, pair<string, GEOSGeometry*>>> nn; // Nearest neighbors
  vector<pair<double, string>> nn;  // Nearest neighbors storing only unique id's for preventing memory leak
 
 set<string> visited; // Set to keep track of visited IDs


for (unsigned int h = 0; h < hashMaps.size(); h++) {

      hashingStart = MPI_Wtime();
    //vector<int> hash = LSHHashs(const_cast<vector<double>*>(&sketch), hashLength, seeds[h]);
     std::vector<int> hash = LSHHash(const_cast<std::vector<double>*>(&sketch), hashLength, &seeds);

//Stop timing for hashing and accumulate the duration
        hashingEnd = MPI_Wtime();
        hashingTotal += hashingEnd - hashingStart;

        // Start timing for querying using MPI_Wtime
        queryingStart = MPI_Wtime();
  // cout << "Hash vector for query hash map " << h << ": ";
//         for (int value : hash) {
//             cout << value << " ";
//         }
//         cout << endl;
 hashmap_key_search_start = MPI_Wtime();
    const list<string>* bucket = hashMaps[h]->get(hash);
    hashmap_key_search_end = MPI_Wtime();
     hashmap_key_search_total += hashmap_key_search_end - hashmap_key_search_start;
    if (bucket != nullptr) {
   // cout << "bucket is not empty" <<endl;
        // Process elements in the bucket
        for (const string& neighborID : *bucket) {
            // Check if the geometry ID has been visited
            if (visited.find(neighborID) == visited.end()) {
                // Check if neighborID exists in geometrySketchMap
                auto neighborSketchIt = geometrySketchMap.find(neighborID);
                if (neighborSketchIt != geometrySketchMap.end()) {
                    visited.insert(neighborID); // Mark this ID as visited
                     GEOSGeometry* neighborGeo = neighborSketchIt->second.first;
                    double distance = jaccardDistance(ctx, queryGeo, neighborGeo);
                   // double distance = SketchJaccardDistance(ctx, queryGeo, neighborGeo);
                    
                    // cout<< "distance is "<<distance << endl;
                   // nn.push_back(make_pair(distance, make_pair(neighborID, neighborGeo))); 
                    nn.push_back(make_pair(distance, neighborID)); 
                  
                }
            }
        }
    }

// //for sketch similarity
// if (bucket != nullptr) {
//             for (const string& neighborID : *bucket) {
//                 if (visited.find(neighborID) == visited.end()) {
//                     auto neighborSketchIt = geometrySketchMap.find(neighborID);
//                     if (neighborSketchIt != geometrySketchMap.end()) {
//                         visited.insert(neighborID); // Mark this ID as visited
//                         const vector<double>& neighborSketch = neighborSketchIt->second.second;  // Using the sketch for the neighbor
//                         
//                         // Calculate sketch distance
//                         double sketchDistance = SketchJaccardDistance(sketch, neighborSketch);
// 
//                         // Store the sketch distance along with the neighbor ID
//                         nn.push_back(make_pair(sketchDistance, neighborID));
//                     }
//                 }
//             }
//         }
   // Stop timing for querying and accumulate the duration
        queryingEnd = MPI_Wtime();
        queryingTotal += queryingEnd - queryingStart;
}

        // Sort the nearest neighbors by distance
        sort(nn.begin(), nn.end(), [](const auto& a, const auto& b) { return a.first < b.first; });

        // Keep only the required number of nearest neighbors
        if (nn.size() > numNeighbors) {
            nn.resize(numNeighbors);
        }
        

        // Store the nearest neighbors for the current query
        LSHoutput.push_back(make_pair(queryID, nn));
        
         // Stop timing for querying and accumulate the duration
        queryingSortEnd = MPI_Wtime();
        queryingSortTotal += queryingSortEnd - queryingStart;

}
// //printing contents of lshoutput
// for (const auto& result : LSHoutput) {
//     const string& queryID = result.first;  // The query ID
//     const auto& neighbors = result.second;  // The vector of nearest neighbors for this query
// 
//     std::cout << "Query ID: " << queryID << std::endl;  // Print the query ID
// 
//     // Check if there are any neighbors found
//     if (neighbors.empty()) {
//         std::cout << "  No neighbors found." << std::endl;
//     } else {
//         for (const auto& neighbor : neighbors) {
//             double distance = neighbor.first;  // The distance to the neighbor
//             const string& neighborID = neighbor.second;  // The neighbor's ID
// 
//             // Print the neighbor's ID and the distance
//             std::cout << "  Neighbor ID: " << neighborID << ", Distance: " << distance << std::endl;
//         }
//     }
// 
//     std::cout << std::endl;  // Extra newline for better readability between query results
// }

// // //for sketch similarity
// // // Initialize variables to store start and end times and total durations
// // double hashingStart, hashingEnd, hashingTotal = 0.0;
// // double queryingStart, queryingEnd, queryingTotal = 0.0;
// // 
// for (const auto& item : QuerygeometrySketchMap) {
//     const string& queryID = item.first;
//     const vector<double>& querySketch = item.second.second;  // Using the sketch for the query
//     
// 
//     vector<pair<double, string>> nn;  // Nearest neighbors storing only unique id's for preventing memory leak
//     set<string> visited; // Set to keep track of visited IDs
// 
//     for (unsigned int h = 0; h < hashMaps.size(); h++) {
//     
//     hashingStart = MPI_Wtime();
//     
//     vector<int> hash = LSHHashs(const_cast<vector<double>*>(&querySketch), hashLength, seeds[h]);
//         
//         // Print the hash values for debugging
// cout << "Hash values for query ID " << queryID << ": ";
// for(int val : hash) {
//     cout << val << " ";
// }
// cout << endl;
//         
//         // Stop timing for hashing and accumulate the duration
//         hashingEnd = MPI_Wtime();
//         hashingTotal += hashingEnd - hashingStart;
// 
//         // Start timing for querying using MPI_Wtime
//         queryingStart = MPI_Wtime();
//         const list<string>* bucket = hashMaps[h]->get(hash);
//         if (bucket != nullptr) {
//             for (const string& neighborID : *bucket) {
//                 if (visited.find(neighborID) == visited.end()) {
//                     auto neighborSketchIt = geometrySketchMap.find(neighborID);
//                     if (neighborSketchIt != geometrySketchMap.end()) {
//                         visited.insert(neighborID); // Mark this ID as visited
//                         const vector<double>& neighborSketch = neighborSketchIt->second.second;  // Using the sketch for the neighbor
//                         
//                         // Calculate sketch distance
//                         double sketchDistance = SketchJaccardDistance(querySketch, neighborSketch);
// 
//                         // Store the sketch distance along with the neighbor ID
//                         nn.push_back(make_pair(sketchDistance, neighborID));
//                     }
//                 }
//             }
//         }
//         // Stop timing for querying and accumulate the duration
//         queryingEnd = MPI_Wtime();
//         queryingTotal += queryingEnd - queryingStart;
//     }
// 
//     // Sort the nearest neighbors by distance
//     sort(nn.begin(), nn.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
// 
//     // Keep only the required number of nearest neighbors
//     if (nn.size() > numNeighbors) {
//         nn.resize(numNeighbors);
//     }
// 
//     // Store the nearest neighbors for the current query
//     LSHoutput.push_back(make_pair(queryID, nn));
// }

double maxhashingTotal;
double maxqueryingTotal;
double maxqueryingSortTotal;
double maxhashmap_key_search_total;
// Use MPI_Reduce to find the maximum hashing time across all processes
MPI_Reduce(&hashingTotal, &maxhashingTotal, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
MPI_Reduce(&queryingTotal,&maxqueryingTotal, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
MPI_Reduce(&queryingSortTotal,&maxqueryingSortTotal, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
MPI_Reduce(&hashmap_key_search_total,&maxhashmap_key_search_total, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);




if (rank == 0) 
{  // Only the process with rank 0 prints the timing
    std::cout << "Max LSH Hashing Time for query geometry: " << maxhashingTotal << " seconds" << std::endl;
    std::cout << "Max time for particular lsh key search in hashmap for: " << hashmap_key_search_total << " seconds" << std::endl;
    std::cout << "Max Querying Time(only hash map search for neighbor) for query geometry: " << maxqueryingTotal << " seconds" << std::endl;
     std::cout << "Max Querying Time(only hash map search for neighbor) for query geometry after sorting: " << maxqueryingSortTotal << " seconds" << std::endl;
}


 // OVERALL time for querying 
    double endTimeQuery = MPI_Wtime();
    double totalTimeQuery = endTimeQuery - startTimeQuery;
    
 double maxTotalTimeQuery;

// Use MPI_Reduce to find the maximum hashing time across all processes
MPI_Reduce(&totalTimeQuery, &maxTotalTimeQuery, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

// The master process (rank 0) prints the maximum hashing time
if (rank == 0) 
{
    std::cout << "Maximum time for query(lsh hashing, querying hashmaps) across all processes: " << maxTotalTimeQuery << " seconds" << std::endl;
}


 
    // Clean up the GEOS context
    GEOS_finish_r(ctx);
    
 
}



//copy geometries + id's
vector<GEOSGeometry *> *copyGeometriesFromMap(const unordered_map<string, GEOSGeometry *> *original, GEOSContextHandle_t ctx) 
{
    vector<GEOSGeometry *> *copy = new vector<GEOSGeometry *>();
    for (const auto &pair : *original) {
        // Make a copy of the geometry
        GEOSGeometry *copyGeometry = GEOSGeom_clone_r(ctx, pair.second);
        copy->push_back(copyGeometry);
    }
    return copy;
}



// Function to print the contents of each HashMap in the vector
void printAllHashMaps(const std::vector<HashMap<std::string>*>& maps) {
    for (size_t i = 0; i < maps.size(); ++i) {
        std::cout << "Hash Map " << i << ":" << std::endl;
        // Call the print method on the current HashMap instance
        maps[i]->print(true); 
    }
}

// void printAllHashMaps(const std::vector<HashMap<std::string>*>& maps, std::ostream& os = std::cout) {
//     for (size_t i = 0; i < maps.size(); ++i) {
//         os << "Hash Map " << i << ":" << std::endl;
//         maps[i]->print(os, true);
//     }
// }

//generate lsh hashes of geos geometries for further querying
void LSHHashFunction(
    const unordered_map<string, GEOSGeometry*>& geos, // Changed to unordered_map
    int rank,
    const unordered_map<string, GEOSGeometry*> centeredQueryGeosWithID,
    unsigned int numNeighbor,
    int gridSize,
    int hashLength ,
   // std::vector<std::pair<std::pair<std::string, GEOSGeom_t*>, std::vector<std::pair<double, std::pair<std::string, GEOSGeom_t*>>>>>& LSHoutput,
    ResultType& LSHoutput,
    GEOSContextHandle_t ctx) 
{

 // Create a new GEOS context handler
   
    GEOSContext_setNoticeHandler_r(ctx, geosMessageHandler);
    GEOSContext_setErrorHandler_r(ctx, geosErrorHandler);
    GEOSWKTWriter* wktWriter = GEOSWKTWriter_create_r(ctx);
    
    //copy original geometries as when original geometries was used, it was resulting in segmentation fault
    unordered_map<string, GEOSGeometry*> geosCopy = copyGeometriesMap(geos, ctx);
  

   // Vector to hold centered geometries
    unordered_map<string, GEOSGeometry*> centeredGeosWithID;

  // Iterate over the unordered_map, center each geometry, and store it with its unique id
    for (const auto& pair : geos) {
    GEOSGeometry* centeredGeo = centerGeometry(ctx, pair.second);
    if (centeredGeo) {
        // Store the centered geometry with its identifier in the map
        centeredGeosWithID[pair.first] = centeredGeo;
    } else {
        std::cerr << "Error centering geometry with ID " << pair.first << std::endl;
    }
   }
   
 //  for writing output in a file
// std::string outputFilename = "cemetery_centered_polygon_" + std::to_string(rank) + ".txt";
// std::ofstream outFile("/home/asbmr/LSH_SIM_CEM/similarity_search_project/geometric-ANN-main/data/out/cem/cemetery_centered_polygon_" + std::to_string(rank) + ".txt", std::ios::out);
// 
// if (!outFile.is_open()) {
//     std::cerr << "Failed to open file for writing results." << std::endl;
//     MPI_Finalize();
//     //return 1; // or handle the error appropriately
// }
// 
// for (const auto& polygons : centeredGeosWithID) {
//     const std::string& polygonID = polygons.first;  // Query ID
//     const auto& geometry = polygons.second;  // Pair of GEOSGeom_t* and vector of doubles
// 
//     // Write the query ID
//     outFile << "polygon ID: " << polygonID << std::endl;
// 
//     // Assuming you want to print the pointer value for GEOSGeom_t*
//     //outFile << "Centered GEOS Geometry : " << geometry << std::endl;
// 
//     //outFile << std::endl << std::endl; // Blank line for separating entries
// char* PolyWKT = GEOSWKTWriter_write_r(ctx, wktWriter, geometry);
//      outFile << "Centered GEOS Geometry : " << PolyWKT << std::endl;
//      GEOSFree_r(ctx, PolyWKT); // Don't forget to free the WKT string after using it
//  
//     
// }

//outFile.close(); // Don't forget to close the file after writing

  //reading centered input geos and its id
    // for (const auto& item : centeredGeosWithID)
//      {
//     std::cout << "ID: " << item.first << ", centered geos geometry inside lshhash funct Centered Geometry: ";
//     printCoordinates(ctx, item.second); // ctx is your GEOSContextHandle_t
//     std::cout << std::endl;
//      }

 //  //copying just geometries
//   vector<GEOSGeometry*> geometries;
//   for (const auto& pair : geosCopy)
//    {
//     geometries.push_back(pair.second); // Add the geometry part of each pair to the vector
//    }


 vector<GEOSGeometry*> geometries;
for (const auto& pair : geosCopy) {
    GEOSGeometry* originalGeometry = pair.second;
    // Call the centerGeometry function on the original geometry
    GEOSGeometry* centeredGeometry = centerGeometry(ctx, originalGeometry);
    if (centeredGeometry != nullptr) {
        // Successfully centered the geometry, add it to the vector
        geometries.push_back(centeredGeometry);
         GEOSGeom_destroy_r(ctx, originalGeometry);
    } else {
        // Centering failed, handle the error accordingly
        cerr << "Warning: Failed to center geometry for ID " << pair.first << endl;
    }
}


// Start timing LSH
    double startTimeLSH = MPI_Wtime();
    
// Start construction timing  for polygon MBR, grids, sketch 
  double startTimeConstruction = MPI_Wtime();

// Start construction timing  for polygon MBR, grids, sketch 
 double startTimeMBR = MPI_Wtime();
  //Create Global MBR parallely by calling GlobalMbr_parallel function
  GEOSGeometry *global_MBR = GlobalMbr_parallel(&geometries, rank);
    
    
  //Creating grid over MBR in root only
  vector<vector<GEOSGeometry *>> global_MBR_grid = createGrid(ctx,
                                            global_MBR, gridSize); 

// OVERALL time for construction 
    double endTimeMBR = MPI_Wtime();
    double totalTimeMBR = endTimeMBR - startTimeMBR;
    
double maxTotalTimeMBR;

// Use MPI_Reduce to find the maximum hashing time across all processes
MPI_Reduce(&totalTimeMBR, &maxTotalTimeMBR, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

// Declare a map to store unique ID, geometry, and its sketch
unordered_map<string, pair<GEOSGeometry*, vector<double>>> geometrySketchMap;
unordered_map<string, pair<GEOSGeometry*, vector<double>>> trimmedGeometrySketchMap;


// The master process (rank 0) prints the maximum hashing time
if (rank == 0) {
    std::cout << "time for construction of globalMBR and grid for input: " << totalTimeMBR << " seconds" << std::endl;
}

//vector<vector<double>> localSketches;


//  //declare vector to store unique id, geometry and its sketch
 // unordered_map<string, pair<GEOSGeometry*, vector<double>>> geometrySketchMap;
 // Start construction timing  for polygon MBR, grids, sketch 
  double startTimeSket = MPI_Wtime();
 // Create sketch for each geometry in the unordered_map centeredQueryGeosWithID
   for (const auto& pair : centeredGeosWithID) 
  {
    const string& id = pair.first; // Unique ID of the geometry
    GEOSGeometry* geometry = pair.second; // Geometry pointer

    vector<double> sketchResult = sketch(ctx, &global_MBR_grid, geometry);
    if (sketchResult.empty()) 
    {
        cerr << "Warning: Sketch generated for geometry with ID " << id << " is empty." << endl;
    } else 
    {
        geometrySketchMap[id] = make_pair(geometry, sketchResult);
    }
  }
  
//   for (const auto& pair : geometrySketchMap) {
//     const std::string& id = pair.first; // Unique ID of the geometry
//     GEOSGeometry* geometry = pair.second.first; // Geometry pointer
//     const std::vector<double>& sketch = pair.second.second; // Sketch vector
// 
//     std::cout << "ID: " << id << ", Geometry Pointer: " << geometry << ", Sketch: ";
// 
//     for (const auto& value : sketch) {
//         std::cout << value << " ";
//     }
//     std::cout << std::endl; // New line after printing each geometry's sketch
// }


//Start construction timing for polygon MBR, grids, sketch 
//double startTimeSket = MPI_Wtime();

//Create sketch for each geometry in the unordered_map centeredQueryGeosWithID
// for (const auto& pair : centeredGeosWithID) {
//     const string& id = pair.first; // Unique ID of the geometry
//     GEOSGeometry* geometry = pair.second; // Geometry pointer
// 
//     vector<double> sketchResult = sketch(ctx, &global_MBR_grid, geometry);
//     if (sketchResult.empty()) {
//         cerr << "Warning: Sketch generated for geometry with ID " << id << " is empty." << endl;
//     } else {
//         geometrySketchMap[id] = make_pair(geometry, sketchResult);
//     }
// }

//Calculate local cell usage

   // Iterate over geometrySketchMap to calculate local cell usage
    int totalCells = gridSize * gridSize;
     std::vector<int> localCellUsage(totalCells, 0);
for (const auto& keyValue : geometrySketchMap) {
    const auto& sketch = keyValue.second.second; // Access the sketch vector from the pair
    for (int i = 0; i < sketch.size(); ++i) {
        if (sketch[i] > 0.0) {  // Define your own threshold for usage
            localCellUsage[i]++;
        }
    }
}

// Aggregate global cell usage
std::vector<int> globalCellUsage(totalCells, 0);
//MPI_Allreduce(localCellUsage.data(), globalCellUsage.data(), gridSize, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
MPI_Allreduce(localCellUsage.data(), globalCellUsage.data(), totalCells, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
//vector<bool> cellsToKeep(totalCells, true);

// Determine cells to keep (i.e., cells that are globally used)
vector<bool> cellsToKeep(totalCells, false);
for (int i = 0; i < globalCellUsage.size(); ++i) {
    if (globalCellUsage[i] > 0.0) {
        cellsToKeep[i] = true;
    }
}

 
// Trim sketches based on globally used cells
for (const auto& pair : geometrySketchMap) {
    const string& id = pair.first;
    GEOSGeometry* geometry = pair.second.first;
    const vector<double>& sketch = pair.second.second;
    vector<double> trimmedSketch;
    // Print the number of elements in the original sketch
   // std::cout << "ID: " << id << " - Original Sketch Size: " << sketch.size() << std::endl;

   

    for (int i = 0; i < sketch.size(); ++i) {
        if (cellsToKeep[i]) {
            trimmedSketch.push_back(sketch[i]);
        }
    }

    trimmedGeometrySketchMap[id] = make_pair(geometry, trimmedSketch);
     // Print the number of elements in the trimmed sketch
    std::cout<< " - Trimmed Sketch Size: " << trimmedSketch.size() << std::endl;

}
 
// 
// // Output the trimmed sketches if needed
// // for (const auto& pair : trimmedGeometrySketchMap) {
// //     const string& id = pair.first;
// //     GEOSGeometry* geometry = pair.second.first;
// //     const vector<double>& trimmedSketch = pair.second.second;
// //     std::cout << "ID: " << id << " Trimmed Sketch: ";
// //     for (const auto& value : trimmedSketch) {
// //         std::cout << value << " ";
// //     }
// //     std::cout << std::endl;
// //    
// // }
// 
// 
// 


// OVERALL time for construction 
    double endTimeSket = MPI_Wtime();
    double totalTimeSket = endTimeSket - startTimeSket;
    
double maxTotalTimeSket;

// Use MPI_Reduce to find the maximum hashing time across all processes
MPI_Reduce(&totalTimeSket, &maxTotalTimeSket, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

// The master process (rank 0) prints the maximum hashing time
if (rank == 0) {
    std::cout << "total time for sketch creation input after trimming: " << maxTotalTimeSket << " seconds" << std::endl;
}

 

// for writing output in a file
std::string outputSketchFilename = "sports_sketch_" + std::to_string(rank) + ".csv";  // Changed extension to .csv
std::ofstream outSKFile("/home/asbmr/LSH_SIM_CEM/similarity_search_project/geometric-ANN-main/data/out/sports/" + outputSketchFilename, std::ios::out);

if (!outSKFile.is_open()) {
    std::cerr << "Failed to open file for writing results." << std::endl;
    MPI_Finalize();
    //return 1; // or handle the error appropriately
}

for (const auto& polygons : trimmedGeometrySketchMap) {
    const std::string& polygonID = polygons.first;  // Query ID
    const auto& geometryAndketch = polygons.second;  // Pair of GEOSGeom_t* and vector of doubles

    // Write the query ID
    outSKFile << polygonID << ",";

    // Write the sketch values, comma-separated
    const std::vector<double>& sketch = geometryAndketch.second;
    for (const auto& value : sketch) {
        outSKFile << value << ",";
    }

    outSKFile << "\n"; // New line for separating entries (removed extra std::endl for blank line)
}

outSKFile.close(); // Don't forget to close the file after writing

// // OVERALL time for construction 
//     double endTimeConstruction = MPI_Wtime();
//     double totalTimeConstruction = endTimeConstruction - startTimeConstruction;
//     
// double maxTotalTimeConstruction;
// 
// // Use MPI_Reduce to find the maximum hashing time across all processes
// MPI_Reduce(&totalTimeConstruction, &maxTotalTimeConstruction, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
// 
// // The master process (rank 0) prints the maximum hashing time
// if (rank == 0) {
//     std::cout << "Maximum time for construction of globalMBR and sktech creation across all processes: " << maxTotalTimeConstruction << " seconds" << std::endl;
// }
// 
// 
// //printing sketch along with its id
// for (const auto& item : geometrySketchMap)
// {
//      const string& id = item.first; // Unique ID
//      const auto& geometryAndSketch = item.second; // Pair of geometry pointer and sketch vector
//       const vector<double>& sketchResult = geometryAndSketch.second; // Extracting the sketch vector
//  
//      cout << "ID: " << id << " - Sketch: ";
//      for (double value : sketchResult) {
//          cout << value << " ";
//      }
//      cout << endl; // End of a single sketch
// }
// 
// 
//start storing only unique id in hashmap instead of polygons for code optimization

// Start Hashing time for polygon 
  double startTimeHashing = MPI_Wtime();
  
// Initialize the vector of hashmaps to store unique IDs  
//vector<unsigned int> seeds = {100 ,360};
vector<unsigned int> seeds = {100 ,360, 550, 880, 900}; // Example seeds
//vector<unsigned int> seeds = {100 ,360, 550, 880, 900, 980, 1080, 1200};
//vector<unsigned int> seeds = {500}; // Example seeds

// int si =1;
// std::vector<HashMap<std::string>*> maps(si);  // Adjust the template parameter to std::string
// for (size_t i = 0; i < si; ++i) {
//     maps[i] = new HashMap<std::string>(5, 5, nullptr);  // Use std::string as the value type for HashMap
// }

 //maps[] = new HashMap<std::string>(5, 5, nullptr);
int nmaps =1;
  //  vector<unsigned int> seeds = {100, 200, 300, 500, 1000}; // Example seeds
   std::vector<HashMap<string>*> maps(nmaps);
for (size_t i = 0; i < nmaps; ++i) {
    maps[i] = new HashMap<string>(1000, 5, nullptr);  // Adjust the template parameter to int
}
double hashingStart, hashingEnd, hashingTotal = 0.0;
double insertStart, insertEnd, insertTotal = 0.0;

// for writing output in a file
// std::string outputHashFilename = "cem_hash_" + std::to_string(rank) + ".csv";  // Changed extension to .csv
// std::ofstream outHashFile("/home/asbmr/similarity_search_project/LSH_SIM_CEM/similarity_search_project/geometric-ANN-main/data/out/cem/" + outputHashFilename, std::ios::out);
// 
// if (!outHashFile.is_open()) {
//     std::cerr << "Failed to open file for writing LSH hashes." << std::endl;
//     MPI_Finalize();
//     return; // Handle the error as appropriate for your application
// }


for (const auto& item : trimmedGeometrySketchMap) {
    string id = item.first;  // Access the unique ID
    const vector<double>& sketch = item.second.second;  // Access the associated sketch

    if (sketch.empty()) {
        std::cerr << "Warning: Sketch for geometry ID " << id << " is empty. Skipping." << std::endl;
        continue;
    }


    for (size_t h = 0; h < nmaps; ++h) {
    hashingStart = MPI_Wtime();
       // std::vector<int> hash = LSHHashs(const_cast<std::vector<double>*>(&sketch), hashLength, seeds[h]);
        std::vector<int> hash = LSHHash(const_cast<std::vector<double>*>(&sketch), hashLength, &seeds);
        
        // Stop timing for LSH hashing and accumulate the duration
    hashingEnd = MPI_Wtime();
    hashingTotal += hashingEnd - hashingStart;
    
   //   // Write the ID and hash values to the file
//         outHashFile << id << ",";
//         for (const auto& value : hash) {
//             outHashFile << value << ",";
//         }
//         outHashFile << "\n";  // End the line for this hash


    // Start timing for insert operation using MPI_Wtime
    insertStart = MPI_Wtime();
        maps[h]->insert(hash, id);  // Insert the ID into the hash map
        
        // Stop timing for insert operation and accumulate the duration
    insertEnd = MPI_Wtime();
    insertTotal += insertEnd - insertStart;
        // cout << "Hash vector for input map " << h<< id<< ": ";
//         for (int value : hash) {
//             cout << value << " ";
//         }cout << endl;
//printAllHashMaps(maps);
    }

}

//printAllHashMaps(maps);

//for hshing and hashmap insertion 
  double maxhashingTotal;
    double maxinsertTotal;

// Use MPI_Reduce to find the maximum hashing time across all processes
MPI_Reduce(&hashingTotal, &maxhashingTotal, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
MPI_Reduce(&insertTotal, &maxinsertTotal, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

if (rank == 0)
{
std::cout << "LSH Hashing Time: " << maxhashingTotal << " seconds" << std::endl;
std::cout << "hash map insertion Time: " << insertTotal << " seconds" << std::endl;
}
//print number of buckets in a table
// for (size_t i = 0; i < maps.size(); ++i) 
// {
//     std::cout << "Table size of HashMap " << i << ": " << maps[i]->bucketCount() << std::endl;
// }

// OVERALL time for construction 
    double endTimeHashing = MPI_Wtime();
    double totalTimeHashing = endTimeHashing - startTimeHashing;
    
double maxTotalTimeHashing;

// Use MPI_Reduce to find the maximum hashing time across all processes
MPI_Reduce(&totalTimeHashing, &maxTotalTimeHashing, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

// The master process (rank 0) prints the maximum hashing time
if (rank == 0)
 {
    std::cout << "Maximum Hashing time including LSH hashing and HashMaps creation across all processes: " << maxTotalTimeHashing << " seconds" << std::endl;
    
  //printAllHashMaps(maps);
//     // Open a file to write the output of hash maps
    //  std::ofstream outFile("/home/asbmr/similarity_search_project/LSH_SIM_CEM/similarity_search_project/geometric-ANN-main/data/out/cem/maps.csv");
//     if (outFile.is_open()) {
//         std::streambuf* coutbuf = std::cout.rdbuf(); // Save old buf
//         std::cout.rdbuf(outFile.rdbuf()); // Redirect std::cout to outFile
// 
//         for (size_t i = 0; i < maps.size(); ++i) {
//             std::cout << "Hash Map " << i << ":" << std::endl;
//             maps[i]->print(true);  // Uses the original print that writes to std::cout
//         }
// 
//         std::cout.rdbuf(coutbuf); // Restore std::cout
//         outFile.close();
//     } else {
//         std::cerr << "Failed to open file for writing." << std::endl;
//     }
}
//printAllHashMaps(maps);




 //for sketch similarity
 //  queryThreadFunct(queryGeos,geometrySketchMap, hashMaps, numNeighbor, hashLength, seeds, rank, LSHoutput, global_MBR_grid);
//std::vector<std::pair<std::pair<std::string, GEOSGeom_t*>, std::vector<std::pair<double, std::pair<std::string, GEOSGeom_t*>>>>> queryResults;
 //for polygon simialrity 
 // queryThreadFunct(queryGeos, hashMaps, numNeighbor, hashLength, seeds, rank, LSHoutput, global_MBR_grid);
 
 
  
 queryThreadFunct(centeredQueryGeosWithID,trimmedGeometrySketchMap, maps, numNeighbor, hashLength, seeds, rank, LSHoutput, global_MBR_grid, cellsToKeep );

  if (LSHoutput.empty()) {
        std::cout << "LSHoutput is empty." << std::endl;
       // return;  // Exit the function early if there's nothing to print
    }
 
}


