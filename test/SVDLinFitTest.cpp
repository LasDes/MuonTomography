#include <iostream>
#include <chrono>
#include <random>

#include <Eigen/SVD>

#include "MutoTypes.h"

int main () {
    int nSimu = 10;
    std::default_random_engine generator{0};
    std::uniform_real_distribution<double> uniform(-1, 1);
    std::normal_distribution<double> normal(0, 0.01);


    // start timing
    auto ts = std::chrono::steady_clock::now();

    Vector3d center (0,0,0);

    for (int i=0; i<nSimu; ++i) {
        // generate random data
        double a = uniform(generator);
        double b = uniform(generator);
        double c = uniform(generator);

        std::cout << "True values are:" << a << " " << b << " " << c << std::endl;

        MatrixXd data(2,  3);
        for (int j=0; j<2; j++) {
            data(j, 0) = center(0) + a * (j-1) * (10.0);// + normal(generator);
            data(j, 1) = center(1) + b * (j-1) * (10.0);// + normal(generator);
            data(j, 2) = center(2) + c * (j-1) * (10.0);// + normal(generator);
        }
        // std::cout << center << std::endl;
        // std::cout << data << std::endl;

        // fit to data 
        double xm = data.col(0).mean();
        double ym = data.col(1).mean();
        double zm = data.col(2).mean();

        for (int j=0; j<2; j++) {
            data(j, 0) -= xm;
            data(j, 1) -= ym;
            data(j, 2) -= zm;
        }
        JacobiSVD<MatrixXd> svd(data, ComputeThinU | ComputeThinV);
        std::cout << "SVD values are:" << svd.matrixV().col(0).transpose() << std::endl;
        // std::cout << "Compare to abc:" << svd.matrixV().col(0)(0) / a << " " 
        //                                << svd.matrixV().col(0)(1) / b << " " 
        //                                << svd.matrixV().col(0)(2) / c << " "  << std::endl;
    }

    // stopwatch
    auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - ts);
    std::cout << nSimu << " simulations took " << time_span.count() << " seconds. " << std::endl;



    return 0;
}