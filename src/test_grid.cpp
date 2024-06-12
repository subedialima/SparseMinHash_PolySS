#define GEOS_USE_ONLY_R_API
#include "geos_c.h"
#include "geoutil.h"
#include "util.h"

int main()
{
    // Create a new GEOS Context Handle
    //
    // While this program itself doesn't require the reentrant API, functions
    // that it calls will expect a Context Handle.
    GEOSContextHandle_t ctx = GEOS_init_r();
    GEOSContext_setNoticeHandler_r(ctx, geosMessageHandler);
    GEOSContext_setErrorHandler_r(ctx, geosErrorHandler);

    printf("If the grid is working, this message will be the only output of "
           "the program!\n");

    GEOSGeometry *base = GEOSGeom_createRectangle_r(ctx, -1, -1, 1, 1);
    vector<vector<GEOSGeometry *>> grid = createGrid(ctx, base, 16);

    // Without having to visualize the grid, here are some properties of the
    // grid that we can expect:
    //     1. The correct number of cells are generated.
    //     2. None of the cells overlap with each other.
    //     3. The union of all of the cells is equal to the base rectangle

    // 1. The correct number of cells are generated.
    int numCells = 0;
    for (auto r = grid.begin(); r != grid.end(); r++)
    {
        numCells += r->size();
    }
    if (numCells != 256)
        printf(
            "Returned number of cells = %d\nExpected number of cells = 256\n",
            numCells);

    // 2. None of the cells overlap with each other.
    // While it may be easier to state "none of the cells overlap with
    // their neigbors", I think that a more extensive check of all-pairs of
    // cells will provide a stronger assurance of correctness.
    //
    // There are may ways to check this. I am currently using the fillRatio
    // function since it has been causing problems with GEOS
    double fill;
    for (unsigned int row = 0; row < grid.size(); row++)
    {
        for (unsigned int col = 0; col < grid[row].size(); col++)
        {
            // Now iterate over the grid again
            for (unsigned int r = 0; r < grid.size(); r++)
            {
                for (unsigned int c = 0; c < grid[r].size(); c++)
                {
                    // Don't compare a cell to itself.
                    if (row == r && col == c)
                        continue;
                    // Find the fillRatio
                    fill = fillRatio(ctx, grid[row][col], grid[r][c]);
                    if (fill != 0.0)
                        printf("Overlap detected between grid[%u][%u] and "
                               "grid[%u][%u]: %lf\n",
                               row, col, r, c, fill);
                }
            }
        }
    }

    // 3. The union of all of the cells is equal to the base rectangle.
    // We can actually check this using the Jaccard Distance. If the union
    // if offset, or has a hole in it, etc, the Jaccard Distance will not be
    // zero.
    GEOSGeometry *temp, *u = NULL;
    for (auto r = grid.begin(); r != grid.end(); r++)
    {
        for (auto e = r->begin(); e != r->end(); e++)
        {
            if (u != NULL)
            {
                temp = u;
                u = GEOSUnion_r(ctx, temp, *e);
                GEOSGeom_destroy_r(ctx, temp);
            }
            else
                u = *e;
        }
    }

    double jd = jaccardDistance(ctx, base, u);
    if (jd != 0.0)
        printf("Reported Jaccard Distance is not 0.0: %lf\n", jd);

    GEOS_finish_r(ctx);

    return 0;
}
