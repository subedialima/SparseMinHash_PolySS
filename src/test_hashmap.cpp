#include "util.h"
#include <iostream>
#include <random>
#include <thread>

using namespace std;

void singleThreadTest(HashMap<int> *map)
{
    for (int i = 0; i < 10; i++)
        map->insert(vector({i}), i);
    map->print(true);

    cout << "******************************" << endl << endl;

    // Now that most buckets will have a value,
    // add more values to those same keys
    random_device rd;
    mt19937 generator(rd());
    uniform_int_distribution<> distrib(0, 9);

    for (int i = 0; i < 20; i++)
        map->insert(vector({distrib(generator)}), i);
    map->print(true);
}

void threadFunc(HashMap<int> *map, int id)
{
    // Randomly generate keys and always insert them with id as value
    random_device rd;
    mt19937 generator(rd());
    uniform_int_distribution<> dist(0, 9);

    for (int i = 0; i < 20; i++)
    {
        vector<int> key = vector<int>({dist(generator), dist(generator)});
        if (!map->tryInsert(key, id))
            printf("Failed to insert (<%d, %d>,  %d), bucket unavailable\n",
                   key[0], key[1], id);
    }
}

void multiThreadTest(HashMap<int> *map)
{
    vector<thread> threads = vector<thread>(10);
    for (int i = 0; i < 10; i++)
    {
        threads[i] = thread(threadFunc, map, i);
    }

    for (int i = 0; i < 10; i++)
    {
        threads[i].join();
    }

    printf("\n");
    map->print(true);

    // Sample 10 random keys from the table
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dist(0, 9);

    for (int i = 0; i < 10; i++)
    {
        vector<int> key = vector<int>({dist(gen), dist(gen)});
        const list<int> *vals = map->get(key);

        if (vals == nullptr)
        {
            printf("Cannot fetch <%d, %d>: No values returned\n", key[0],
                   key[1]);
        }
        else
        {
            printf("Fetched <%d, %d>: ", key[0], key[1]);
            for (auto v : *vals)
                cout << v << " ";
            cout << "\n";
        }
    }
}

int main()
{
    // Single threaded tests
    HashMap<int> test(11, 1);
    singleThreadTest(&test);
    printf("Hash Map Size: %d\n", test.size());

    // Multi-threaded tests
    // Also increase the hash length
    HashMap<int> test2(17, 2);
    multiThreadTest(&test2);
    printf("Hash Map Size: %d\n", test2.size());

    return 0;
}
