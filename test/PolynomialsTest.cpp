#include <iostream>
#include <chrono>
#include <random>

#include "MutoTypes.h"

#include <unsupported/Eigen/Polynomials>

int main (int argc, char** argv) {
    int nSimu = 1000000;
    std::random_device rd;
    std::default_random_engine generator{rd()};
    std::uniform_real_distribution<double> distribution(-1, 1.0);

    // start timing
    auto ts = std::chrono::steady_clock::now();

    PolynomialSolver<double, 3> solver;

    for (int i=0; i<nSimu; ++i) {
        Vector4d poly;
        std::vector<double> roots;
        for (int j=0; j<4; ++j) { poly[j] = distribution(generator); }

        solver.compute(poly);
        solver.realRoots(roots);

        // for (int j=0; j<4; ++j) { std::cout << poly[j] << ", "; }
        // std::cout<<std::endl;
        // for (int j=0; j<3; ++j) { std::cout << "root " << j << "=" << roots[j] << " "; }
        // std::cout<<std::endl<<std::endl;
    }

    // stopwatch
    auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - ts);
    std::cout << "solve polynomials took: " << time_span.count() << " seconds. " << std::endl;

    return 0;
}