#include <iostream>
#include <chrono>
#include <random>

#include "MutoSiddon.h"

int main () {
    int nSimu = 10000;
    std::default_random_engine generator{0};
    std::uniform_real_distribution<double> distribution(-500, 500.0);

    VoxelGrid grid;
    grid.x_min = grid.y_min = grid.z_min= -100.0;
    grid.dx = grid.dy = grid.dz = 2.0;
    grid.nx = grid.ny = grid.nz = 100;

    MutoSiddon siddon (grid);

    Vector3 p1;
    Vector3 p2;

    // start timing
    auto ts = std::chrono::steady_clock::now();

    std::vector<VoxelData> ans;

    for (int i=0; i<nSimu; ++i) {
        p1(0) = distribution(generator);
        p1(1) = distribution(generator);
        p1(2) = distribution(generator);
        p2(0) = distribution(generator);
        p2(1) = distribution(generator);
        p2(2) = distribution(generator);
        ans = siddon.getVoxelPath(p1, p2);

        if (ans.size() != 0) {
            std::cout << "total points: " << ans.size() << std::endl;

            std::cout << "In points: " << p1(0) << " " << p1(1) << " " << p1(2) << std::endl;
            std::cout << "Out points: " << p2(0) << " " << p2(1) << " " << p2(2) << std::endl;
            
            for (VoxelData i : ans) {
                std::cout << i.x << " " << i.y << " " << i.z << " : " << i.length << std::endl;
            }
        }
        
    }

    // stopwatch
    auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - ts);
    // std::cout << nSimu << " simulations took " << time_span.count() << " seconds." << std::endl;



    return 0;
}