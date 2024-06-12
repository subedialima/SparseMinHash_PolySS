//#include "parse_geodata.h"
#include "query.h"
#include "thread.h"
#include <cmath>
#include <list>
#include <mpi.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <map> 
#include <set>
//#include <geos_c.h>
#include <ostream>
using namespace std;

//include <geos/export.h>

void BFquery(int rank,
             const unordered_map<string, GEOSGeometry*>& queryGeos,
             const unordered_map<string, GEOSGeometry*>& geos,
             unsigned int numNeighbors,
             vector<pair<string, vector<pair<double, string>>>>& queryResults)
{

    GEOSContextHandle_t ctx = GEOS_init_r();
    GEOSContext_setNoticeHandler_r(ctx, geosMessageHandler);
    GEOSContext_setErrorHandler_r(ctx, geosErrorHandler);

    unordered_map<string, GEOSGeometry*> centeredGeosWithID;
    // Center each geometry and store it with its identifier
    for (const auto& pair : geos)
      {
        GEOSGeometry* centeredGeo = centerGeometry(ctx, pair.second);
        if (centeredGeo) {
            centeredGeosWithID[pair.first] = centeredGeo;
        } 
        else
        {
            std::cerr << "Error centering geometry with ID " << pair.first << std::endl;
        }
      }

 // Start timing BF
    double startTimeBF = MPI_Wtime();
    
    for (const auto& queryPair : queryGeos) {
    const string& queryID = queryPair.first;
    GEOSGeometry* queryGeo = queryPair.second;
    vector<pair<double, string>> nn;  // Nearest neighbors

    for (const auto& candidatePair : centeredGeosWithID)
      {
        const string& candidateID = candidatePair.first;
        GEOSGeometry* candidateGeo = candidatePair.second;
        double distance = jaccardDistance(ctx, queryGeo, candidateGeo);

        nn.emplace_back(distance, candidateID);
      }

    // Sort the nearest neighbors by distance
    sort(nn.begin(), nn.end());

    // Keep only the top numNeighbors neighbors
    nn.resize(min(nn.size(), static_cast<size_t>(numNeighbors)));

    queryResults.push_back(make_pair(queryID, nn));
}

MPI_Barrier(MPI_COMM_WORLD) ;

    double endTimeBF = MPI_Wtime();
    double totalTimeBF = endTimeBF - startTimeBF;


// Print the timing results
    if (rank == 0) {  // Only the master process (rank 0) prints the results
        std::cout << "BF Query Processing Time: " << totalTimeBF << " seconds" << std::endl;
    }

    GEOS_finish_r(ctx);
}

