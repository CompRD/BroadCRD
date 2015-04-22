#include "seednet/SeedClient.h"

int main(int argc, char** argv)
{
    SeedClient sc("localhost");
    for (int x = 0; x < 1000; x++)
    {
        seed_idx seed;
        seed = sc.get_next_seed();
        printf("got seed: %d\n", seed);
        if (seed == -1)
            break;
        sc.set_seed_status(seed, SEED_DONE);
        printf("set seed to done: %d\n", seed);
    }
    sc.shutdown();
    printf("shutdown sever\n");
    return 0;
}
