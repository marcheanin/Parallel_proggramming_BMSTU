#include <mpi.h>
#include <cstdio>
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <random>
#include <chrono>

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

int* get_interval(int proc, int size, const int* interval)
{
    int* range = new int[2];
    int interval_size = (interval[1] - interval[0]) / size;

    range[0] = interval[0] + interval_size * proc;
    range[1] = interval[0] + interval_size * (proc + 1);
    range[1] = range[1] == interval[1] - 1 ? interval[1] : range[1];
    return range;
}

void multiprocess_mult(std::vector <std::vector <double> > A, std::vector <double> x, std::vector <double>& res, int count, int rank) {
    if (rank == 0){
        int chunk_size = A.size() / count, bonus = A.size() - chunk_size * count;
        int to_thread = 1;
        for (int start = 0, end = chunk_size; start < A.size(); start = end, end = start + chunk_size){
            if (bonus) {
                end++;
                bonus--;
            }
            std::cout << start << " " << end - 1 << " " << to_thread << std::endl;
            int interval[2] = {start, end - 1};
            MPI_Send(&interval, 2, MPI_INT, to_thread, 0, MPI_COMM_WORLD);
            to_thread++;
            for (int i = 1; i < count; i++){
                int* recieved_res;
                MPI_Recv(&recieved_res, A[0].size(), MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (int i = 0; i < res.size(); i++){
                    res[i] += recieved_res[i];
                }
            }
        }
    }
    else{
        int interval[2];
        MPI_Recv(&interval, 2, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        std::vector <double> part_res;
        mult(A, x, interval[0], interval[1], part_res);
        int* send_res = new int[A[0].size()];
        for (int i = interval[0]; i <= interval[1]; i++){
            MPI_Send(send_res, A[0].size(), MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
    }
}

std::vector <double> solve(std::vector <std::vector <double> > A, double eps, int col_ranks, int rank){
    uint n = A.size(), m = A[0].size();
    std::vector <double> x(m, 5);
    //std::vector <double> b(m, n+1);
    std::vector <double> b (n);
    for (double & i : b){
        i = n+1;
    }
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

    MPI_Init(&argc, &argv);

    int rank = -1, col_ranks = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &col_ranks);

    auto go = std::chrono::high_resolution_clock::now();

//    std::vector <std::vector <double> > a(10, std::vector <double> (10));
//    for (int i = 0; i < a.size(); i++){
//        for (int j = 0; j < a[0].size(); j++){
//            if (i == j) a[i][j] = 2.0;
//            else a[i][j] = 1.0;
//        }
//    }
//
//    printMatrix(a);
//    auto res = solve(a, 0.00001, col_ranks, rank);
//    for (double re : res){
//        std::cout << re << " ";
//    }
//    std::cout << std::endl;
//
//    auto finish = std::chrono::high_resolution_clock::now();
//
//    auto duration = duration_cast<std::chrono::microseconds>(finish - go);
//    std::cout << "Time: "<< duration.count() << " ms" << std::endl;

    std::vector <std::vector <double> > a(3, std::vector <double> (3));
    for (int i = 0; i < a.size(); i++){
        for (int j = 0; j < a[0].size(); j++){
            if (i == j) a[i][j] = 2.0;
            else a[i][j] = 1.0;
        }
    }
    std::vector <double> x = {1, 2, 3};

    std::vector <double> res;
    multiprocess_mult(a, x, res, col_ranks, rank);

    for (double re : res){
       std::cout << re << " ";
    }
    MPI_Finalize();

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