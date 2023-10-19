#include <mpi.h>
#include <cstdio>
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <random>

double norm(std::vector <double> v){
    double sum = 0;
    for (auto elem : v){
        sum += elem * elem;
    }
    return sqrt(sum);
}

double vector_mult(std::vector <double> a, std::vector <double> b){
    assert(a.size() == b.size());
    double res = 0;
    for (int i = 0; i < a.size(); i++){
        res += a[i] * b[i];
    }
    return res;
}

void mult(std::vector <std::vector <double> > A, std::vector <double> x, int raw1, int raw2, std::vector <double>& res){
    double sum;
    for (int i = raw1; i <= raw2; i++){
        sum = 0;
        for (int j = 0; j < x.size(); j++){
            sum += A[i][j] * x[j];
        }
        res[i] = sum;
    }
}

std::vector <double> mult_digit_vector(double a, std::vector <double> b) {
    for (int i = 0; i < b.size(); i++){
        b[i] *= a;
    }
    return b;
}

std::vector <double> vector_diff(std::vector <double> a, std::vector <double> b) {
    assert(a.size() == b.size());
    for (int i = 0; i < a.size(); i++){
        a[i] -= b[i];
    }
    return a;
}

void printMatrix(const std::vector <std::vector <double> >& mt){
    for (int i = 0; i < mt.size(); i++){
        for (int j = 0; j < mt[0].size(); j++){
            std::cout << mt[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

void fillMatrixWithRandom(std::vector <std::vector <double> >& mt) {
    std::random_device rd;
    std::mt19937 gen(rd());
    for (int i = 0; i < mt.size(); i++){
        for (int j = 0; j < mt[0].size(); j++){
            mt[i][j] = (1.0 + gen() % 100) / 10.f;
        }
    }
}

std::vector <double> solve(std::vector <std::vector <double> > A, double eps){
    int n = A.size(), m = A[0].size();
    std::vector <double> x(m, 0.1);
    //std::vector <double> b(m, n+1);
    std::vector <double> b = {9, 8, 7};
    std::vector <double> res(m, 0);
    std::vector <double> Ay(m, 0);

    double tau, error;
    int iter = 1;
    while(1){
        mult(A, x, 0, n-1, res);
        auto y = vector_diff(res, b);
        mult(A, y, 0, n-1, Ay);
        tau = vector_mult(y, Ay) / vector_mult(Ay, Ay);
        x = vector_diff(x, mult_digit_vector(tau, y));

        std::cout << "tau: " << tau << std::endl;

        mult(A, x, 0, n-1, res);

        error = norm(vector_diff(res, b)) / norm(b);
        std::cout << iter << " " << error << std::endl;
        iter++;

        if (error < eps) break;
    }

    return x;
}

int main(int argc, char** argv) {
//    // Initialize the MPI environment
//    MPI_Init(nullptr, nullptr);
//
//    // Get the number of processes
//    int world_size;
//    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
//
//    // Get the rank of the process
//    int world_rank;
//    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
//
//    // Get the name of the processor
//    char processor_name[MPI_MAX_PROCESSOR_NAME];
//    int name_len;
//    MPI_Get_processor_name(processor_name, &name_len);
//
//    // Print off a hello world message
//    printf("Hello world from processor %s, rank %d out of %d processors\n",
//           processor_name, world_rank, world_size);
//
//    // Finalize the MPI environment.
//    MPI_Finalize();

    std::vector <std::vector <double> > a(3, std::vector <double> (3));
    for (int i = 0; i < a.size(); i++){
        for (int j = 0; j < a[0].size(); j++){
            std::cin >> a[i][j];
        }
    }
    //fillMatrixWithRandom(a);
    printMatrix(a);
    auto res = solve(a, 0.00001);
    for (double re : res){
        std::cout << re << " ";
    }
    std::cout << std::endl;


}

/*
2.0 1.0 1.0 1.0
1.0 2.0 1.0 1.0
1.0 1.0 2.0 1.0
1.0 1.0 1.0 2.0

1 2 3
1 2 1
3 1 3

2.0 1.0 1.0
1.0 2.0 1.0
1.0 1.0 2.0
 */