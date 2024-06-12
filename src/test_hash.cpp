#define GEOS_USE_ONLY_R_API
#include "geos_c.h"
#include "geoutil.h"
#include "util.h"
#include <fstream>
#include <random>
#include <sys/stat.h>
#include <unistd.h>
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define GRID_SIZE 64
#define HASH_LENGTH 1
#define RANDOM_SEED 42
#define REPORT_INTERVAL 100000

GEOSContextHandle_t ctx;
GEOSGeometry *MBR, *p1, *p2, *p3, *p4;
vector<vector<GEOSGeometry *>> grid;

void terminate(int parameter)
{
    GEOSGeom_destroy_r(ctx, MBR);
    GEOSGeom_destroy_r(ctx, p1);
    GEOSGeom_destroy_r(ctx, p2);
    GEOSGeom_destroy_r(ctx, p3);
    GEOSGeom_destroy_r(ctx, p4);

    for (auto r = grid.begin(); r != grid.end(); r++)
        for (auto e = r->begin(); e != r->end(); e++)
            GEOSGeom_destroy_r(ctx, *e);

    GEOS_finish_r(ctx);

    exit(0);
}

double vectorJS(vector<double> one, vector<double> two)
{
    double numerator = 0;
    double denominator = 0;
    for (unsigned int i = 0; i < one.size(); i++)
    {
        numerator += MIN(one[i], two[i]);
        denominator += MAX(one[i], two[i]);
    }

    return numerator / denominator;
}

int main()
{
    // The purpose of this program is to test the LSHHash function. For
    // two polygons with a known Jaccard Similarity, it is expected that
    // the probability of the hashes being equal will converge to their
    // Jaccard Similarity.

    // Create a new GEOS Context Handle
    //
    // While this program itself doesn't require the reentrant API, functions
    // that it calls will expect a Context Handle.
    ctx = GEOS_init_r();
    GEOSContext_setNoticeHandler_r(ctx, geosMessageHandler);
    GEOSContext_setErrorHandler_r(ctx, geosErrorHandler);

    // Set the SIGINT handler
    struct sigaction act;
    act.sa_handler = terminate;
    act.sa_flags = 0;
    sigemptyset(&act.sa_mask);
    sigaction(SIGINT, &act, NULL);

    // Too test this we need a small dataset, prehaps with four shapes.
    // Fortunatly, the output of the rest of the project is exactly that.

    string wkt;
    ifstream f("./data/out/wkt/2q.txt");
    GEOSWKTReader *parser = GEOSWKTReader_create_r(ctx);
    getline(f, wkt);
    p1 = GEOSWKTReader_read_r(ctx, parser, wkt.c_str());
    getline(f, wkt);
    p2 = GEOSWKTReader_read_r(ctx, parser, wkt.c_str());

    // Read most disimilar two
    getline(f, wkt);
    p3 = GEOSWKTReader_read_r(ctx, parser, wkt.c_str());
    getline(f, wkt);
    p4 = GEOSWKTReader_read_r(ctx, parser, wkt.c_str());
    f.close();
    GEOSWKTReader_destroy_r(ctx, parser);

    // Report inital similarities
    double p12_js = 1 - jaccardDistance(ctx, p1, p2);
    double p34_js = 1 - jaccardDistance(ctx, p3, p4);
    printf("P1-2 Jaccard Similarity: %lf\n", p12_js);
    printf("P3-4 Jaccard Similarity: %lf\n", p34_js);

    // Generate sketches of the polygons
    MBR =
        minimumBoundingRectangle(ctx, vector<GEOSGeometry *>({p1, p2, p3, p4}));
    grid = createGrid(ctx, MBR, GRID_SIZE);
    vector<double> p1Sketch = sketch(ctx, &grid, p1);
    vector<double> p2Sketch = sketch(ctx, &grid, p2);
    vector<double> p3Sketch = sketch(ctx, &grid, p3);
    vector<double> p4Sketch = sketch(ctx, &grid, p4);

    // Report the Jaccard Similarity between the vector sketchs
    double p12_sjs = vectorJS(p1Sketch, p2Sketch);
    double p34_sjs = vectorJS(p3Sketch, p4Sketch);
    printf("P1-P2 Sketch Jaccard Similarity: %lf\n", p12_sjs);
    printf("P3-P4 Sketch Jaccard Similarity: %lf\n", p34_sjs);

    // Generate the seeds for the hash function
    minstd_rand generator(RANDOM_SEED);
    vector<unsigned int> seeds(HASH_LENGTH);

    // Run the test!
    printf("This program is entering an infinite loop! Use C-c to break.\n");
    int p12Equal = 0;
    int p34Equal = 0;
    for (int count = 1; true; count++)
    {
        for (int i = 0; i < HASH_LENGTH; i++)
        {
            seeds[i] = generator();
        }

        // Generate new seeds
        vector<int> p1Hash = LSHHash(&p1Sketch, HASH_LENGTH, &seeds);
        vector<int> p2Hash = LSHHash(&p2Sketch, HASH_LENGTH, &seeds);
        if (p1Hash == p2Hash)
            p12Equal += 1;

        vector<int> p3Hash = LSHHash(&p3Sketch, HASH_LENGTH, &seeds);
        vector<int> p4Hash = LSHHash(&p4Sketch, HASH_LENGTH, &seeds);
        if (p3Hash == p4Hash)
            p34Equal += 1;

        if (count % REPORT_INTERVAL == 0)
        {
            printf(
                "%06d: p12 Equal Rate = %lf (target: %lf, difference: %lf)\n",
                count, p12Equal / (double)count, p12_sjs,
                p12_sjs - (p12Equal / (double)count));
            printf(
                "%06d: p34 Equal Rate = %lf (target: %lf, difference: %lf)\n",
                count, p34Equal / (double)count, p34_sjs,
                p34_sjs - (p34Equal / (double)count));
        }
    }

    return 0;
}
