#define GEOS_USE_ONLY_R_API
#include "geos_c.h"
#include "geoutil.h"
#include "parse_geodata.h"
#include "util.h"
#include <random>
#include <thread>
#define THREAD_COUNT 36
#define RANDOM_SEED 42
#define SAMPLE_SIZE 50
#define RAND_DURATION 10000

vector<GEOSGeometry *> g;
vector<vector<GEOSGeometry *>> grid;

void threadFunc()
{
    GEOSContextHandle_t ctx = GEOS_init_r();
    GEOSContext_setNoticeHandler_r(ctx, geosMessageHandler);
    GEOSContext_setErrorHandler_r(ctx, geosErrorHandler);

    random_device rd;

    char out[10];
    int rand_dur = rd() % RAND_DURATION;

    for (auto geo = g.begin(); geo != g.end(); geo++)
    {
        vector<double> gSketch = sketch(ctx, &grid, *geo);

        rand_dur--;
        if (rand_dur < 0)
        {
            snprintf(out, 10, "%lf\n", gSketch[rd() % SAMPLE_SIZE]);
            write(1, out, 9);
            rand_dur = rd() % RAND_DURATION;
        }
    }

    GEOS_finish_r(ctx);
}

int main()
{
    // Create a new GEOS Context Handle
    //
    // While this program itself doesn't require the reentrant API, functions
    // that it calls will expect a Context Handle.
    GEOSContextHandle_t ctx = GEOS_init_r();
    GEOSContext_setNoticeHandler_r(ctx, geosMessageHandler);
    GEOSContext_setErrorHandler_r(ctx, geosErrorHandler);

    minstd_rand generator(RANDOM_SEED);
    uniform_real_distribution dist(-1.0, 1.0);

    // g = vector<GEOSGeometry*>(SAMPLE_SIZE);
    // for (int i = 0; i < SAMPLE_SIZE; i++)
    // 	g[i] = GEOSGeom_createRectangle_r(ctx, dist(generator), dist(generator),
    // dist(generator), dist(generator));

    g = *read_wkt("./data/in/wkt/lakes_mircoscopic.txt", ctx);
    for (auto p = g.begin(); p != g.end(); p++)
        *p = centerGeometry(ctx, *p);

    GEOSGeometry *mbr = minimumBoundingRectangle(ctx, g);
    grid = createGrid(ctx, mbr, 50);

    vector<thread> threads(THREAD_COUNT);
    for (int i = 0; i < THREAD_COUNT; i++)
        threads[i] = thread(threadFunc);

    for (int i = 0; i < THREAD_COUNT; i++)
        threads[i].join();

    GEOS_finish_r(ctx);

    return 0;
}
