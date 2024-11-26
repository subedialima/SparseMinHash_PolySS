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
                     ResultType& LSHoutput,
                     vector<vector<GEOSGeometry*>>& global_MBR_grid,
                      const vector<bool>& cellsToKeep) 
{

  // Create a new GEOS context handler
  GEOSContextHandle_t ctx = GEOS_init_r();
  GEOSContext_setNoticeHandler_r(ctx, geosMessageHandler);
  GEOSContext_setErrorHandler_r(ctx, geosErrorHandler);
    
 // Generate sketches for Query Geometries
  unordered_map<string, pair<GEOSGeometry*, vector<double>>> QuerygeometrySketchMap;

// Start Sketching time for query polygon 
 double startTimeQuerysketch = MPI_Wtime();
  

// Create trimmed sketch for each geometry in the unordered_map centeredQueryGeosWithID
for (const auto& pair : centeredQueryGeosWithID)
  {
    const string& id = pair.first; // Unique ID of the geometry
    GEOSGeometry* geometry = pair.second; // Geometry pointer

    vector<double> fullSketch = sketch(ctx, &global_MBR_grid, geometry); // Generate the full sketch first
    vector<double> trimmedSketch;
    for (size_t i = 0; i < fullSketch.size(); ++i) {
        if (cellsToKeep[i])
        {
            trimmedSketch.push_back(fullSketch[i]); 
            // Keep only the cells marked by cellsToKeep
        }
         
    }
    // std::cout << "ID: " << id << " - Trimmed Sketch Size: " << trimmedSketch.size() << std::endl;
    if (trimmedSketch.empty())
     {
        cerr << "Warning: Trimmed sketch generated for query geometry with ID " << id << " is empty." << endl;
     } 
    else 
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

 vector<pair<double, string>> nn;  // Nearest neighbors storing only unique id's for preventing memory leak
 
 set<string> visited; // Set to keep track of visited IDs


for (unsigned int h = 0; h < hashMaps.size(); h++) {

  hashingStart = MPI_Wtime();
  
  std::vector<int> hash = LSHHash(const_cast<std::vector<double>*>(&sketch), hashLength, &seeds);

  //Stop timing for hashing and accumulate the duration
  hashingEnd = MPI_Wtime();
  hashingTotal += hashingEnd - hashingStart;

  // Start timing for querying using MPI_Wtime
  queryingStart = MPI_Wtime();
  
  hashmap_key_search_start = MPI_Wtime();
  
    const list<string>* bucket = hashMaps[h]->get(hash);
  
    hashmap_key_search_end = MPI_Wtime();
  
     hashmap_key_search_total += hashmap_key_search_end - hashmap_key_search_start;
    if (bucket != nullptr) 
    {
        // Process elements in the bucket
        for (const string& neighborID : *bucket)
          {
            // Check if the geometry ID has been visited
            if (visited.find(neighborID) == visited.end()) 
            {
                // Check if neighborID exists in geometrySketchMap
                auto neighborSketchIt = geometrySketchMap.find(neighborID);
                if (neighborSketchIt != geometrySketchMap.end()) 
                {
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
{  
  // Only the process with rank 0 prints the timing
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

//generate lsh hashes of geos geometries for further querying
void LSHHashFunction(
    const unordered_map<string, GEOSGeometry*>& geos, // Changed to unordered_map
    int rank,
    const unordered_map<string, GEOSGeometry*> centeredQueryGeosWithID,
    unsigned int numNeighbor,
    int gridSize,
    int hashLength ,
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
    } else
    {
        std::cerr << "Error centering geometry with ID " << pair.first << std::endl;
    }
   }


 vector<GEOSGeometry*> geometries;
 for (const auto& pair : geosCopy) 
 {
    GEOSGeometry* originalGeometry = pair.second;
    // Call the centerGeometry function on the original geometry
    GEOSGeometry* centeredGeometry = centerGeometry(ctx, originalGeometry);
    if (centeredGeometry != nullptr) {
        // Successfully centered the geometry, add it to the vector
        geometries.push_back(centeredGeometry);
         GEOSGeom_destroy_r(ctx, originalGeometry);
    } else 
    {
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
  vector<vector<GEOSGeometry *>> global_MBR_grid = createGrid(ctx,global_MBR, gridSize); 

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
if (rank == 0) 
{
    std::cout << "time for construction of globalMBR and grid for input: " << totalTimeMBR << " seconds" << std::endl;
}

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
     } 
     else 
     {
        geometrySketchMap[id] = make_pair(geometry, sketchResult);
     }
  }
  

  // Iterate over geometrySketchMap to calculate local cell usage
  int totalCells = gridSize * gridSize;
  std::vector<int> localCellUsage(totalCells, 0);
  for (const auto& keyValue : geometrySketchMap) 
  {
    const auto& sketch = keyValue.second.second; // Access the sketch vector from the pair
    for (int i = 0; i < sketch.size(); ++i) 
    {
        if (sketch[i] > 0.0) 
        {  // Define your own threshold for usage
            localCellUsage[i]++;
        }
    }
  }

// Aggregate global cell usage from all local processes
std::vector<int> globalCellUsage(totalCells, 0);  
MPI_Allreduce(localCellUsage.data(), globalCellUsage.data(), totalCells, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

// Determine cells to keep (i.e., cells that are globally used)
vector<bool> cellsToKeep(totalCells, false);
for (int i = 0; i < globalCellUsage.size(); ++i)
  {
    if (globalCellUsage[i] > 0.0) {
        cellsToKeep[i] = true;
    }
  }

 
// Trim sketches based on globally used cells
for (const auto& pair : geometrySketchMap)
  {
    const string& id = pair.first;
    GEOSGeometry* geometry = pair.second.first;
    const vector<double>& sketch = pair.second.second;
    vector<double> trimmedSketch;
    for (int i = 0; i < sketch.size(); ++i) 
    {
        if (cellsToKeep[i]) {
            trimmedSketch.push_back(sketch[i]);
        }
    }

    trimmedGeometrySketchMap[id] = make_pair(geometry, trimmedSketch);
    // Print the number of elements in the trimmed sketch
    std::cout<< " - Trimmed Sketch Size: " << trimmedSketch.size() << std::endl;

}
 
// OVERALL time for construction 
double endTimeSket = MPI_Wtime();
double totalTimeSket = endTimeSket - startTimeSket;
    
double maxTotalTimeSket;

// Use MPI_Reduce to find the maximum hashing time across all processes
MPI_Reduce(&totalTimeSket, &maxTotalTimeSket, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

// The master process (rank 0) prints the maximum hashing time
if (rank == 0)
{
    std::cout << "total time for sketch creation input after trimming: " << maxTotalTimeSket << " seconds" << std::endl;
}


// Start Hashing time for polygon 
double startTimeHashing = MPI_Wtime();
  
// Initialize the vector of seeds for random number generation 
vector<unsigned int> seeds = {100 ,360, 550, 880, 900}; // Example seeds

int nmaps =1;
std::vector<HashMap<string>*> maps(nmaps);
for (size_t i = 0; i < nmaps; ++i) 
{
    maps[i] = new HashMap<string>(1000, 5, nullptr);  // Adjust the template parameter to int
}
double hashingStart, hashingEnd, hashingTotal = 0.0;
double insertStart, insertEnd, insertTotal = 0.0;


for (const auto& item : trimmedGeometrySketchMap) 
{
    string id = item.first;  // Access the unique ID
    const vector<double>& sketch = item.second.second;  // Access the associated sketch

    if (sketch.empty())
    {
        std::cerr << "Warning: Sketch for geometry ID " << id << " is empty. Skipping." << std::endl;
        continue;
    }

    for (size_t h = 0; h < nmaps; ++h) {
     hashingStart = MPI_Wtime();
      
    std::vector<int> hash = LSHHash(const_cast<std::vector<double>*>(&sketch), hashLength, &seeds);
        
    // Stop timing for LSH hashing and accumulate the duration
    hashingEnd = MPI_Wtime();
    hashingTotal += hashingEnd - hashingStart;
    
    // Start timing for insert operation using MPI_Wtime
    insertStart = MPI_Wtime();
    maps[h]->insert(hash, id);  // Insert the ID into the hash map
        
    // Stop timing for insert operation and accumulate the duration
    insertEnd = MPI_Wtime();
    insertTotal += insertEnd - insertStart;
    }
}

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
}

 //for polygon simialrity 
 queryThreadFunct(centeredQueryGeosWithID,trimmedGeometrySketchMap, maps, numNeighbor, hashLength, seeds, rank, LSHoutput, global_MBR_grid, cellsToKeep );

   if (LSHoutput.empty())
   {
        std::cout << "LSHoutput is empty." << std::endl;
       // return;  // Exit the function early if there's nothing to print
    }
 
}


