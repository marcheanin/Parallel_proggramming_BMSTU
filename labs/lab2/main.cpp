#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include "mpich/mpi.h"
#include "mpich/mpi_proto.h"
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


int main(int argc, char **argv) {
    int matrix_size = 500;
    //std::cin >> matrix_size;

    std::vector <std::vector <double> > A(matrix_size, std::vector <double> (matrix_size));
    for (int i = 0; i < A.size(); i++){
        for (int j = 0; j < A[0].size(); j++){
            if (i == j) A[i][j] = 2.0;
            else A[i][j] = 1.0;
        }
    }
    std::vector <double> x (matrix_size, 5);
    int rank = -1, col_ranks = -1;

    uint n = A.size(), m = A[0].size();
    //std::vector <double> b(m, n+1);
    std::vector <double> res(n, 0.0);
    std::vector <double> res2(n, 0.0);
    std::vector <double> Ay(n, 0.0);
    std::vector <double> b (n);
    for (double & i : b){
        i = n+1;
    }

    int chunk_size, bonus;
    std::vector <std::pair <int, int> > intervals;
    double eps = 0.00001;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &col_ranks);

    auto go = std::chrono::high_resolution_clock::now();

    chunk_size = A.size() / col_ranks, bonus = A.size() - chunk_size * col_ranks;
    for (int start = 0, end = chunk_size; start < A.size(); start = end, end = start + chunk_size) {
        if (bonus) {
            end++;
            bonus--;
        }
        intervals.emplace_back(start, end-1);
    }
    std::vector <double> send_res (n);
    std::vector <double> collect_res (n);
    std::vector <double> y (n);
    for (int k = 0; k < 30; k++) {
        if (rank == 0) {
            double tau, error;

            mult(A, x, intervals[rank].first, intervals[rank].second, res);
            send_res = x;
            for (int i = 1; i < col_ranks; i++) {
                MPI_Send(&send_res[0], send_res.size(), MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
            }

            for (int i = 1; i < col_ranks; i++) {
                MPI_Recv(&collect_res[0], n, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (int j = intervals[i].first; j <= intervals[i].second; j++) {
                    res[j] = collect_res[j];
                }
            }
            y = vector_diff(res, b);

            mult(A, y, intervals[rank].first, intervals[rank].second, Ay);
            send_res = x;

            for (int i = 1; i < col_ranks; i++) {
                MPI_Send(&send_res[0], send_res.size(), MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
            }
            for (int i = 1; i < col_ranks; i++) {
                MPI_Recv(&collect_res[0], n, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (int j = intervals[i].first; j <= intervals[i].second; j++) {
                    Ay[j] = collect_res[j];
                }
            }

            tau = vector_mult(y, Ay) / vector_mult(Ay, Ay);

            x = vector_diff(x, mult_digit_vector(tau, y));

            mult(A, x, intervals[rank].first, intervals[rank].second, res2);
            send_res = x;
            for (int i = 1; i < col_ranks; i++) {
                MPI_Send(&send_res[0], send_res.size(), MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
            }
            for (int i = 1; i < col_ranks; i++) {
                MPI_Recv(&collect_res[0], n, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                for (int j = intervals[i].first; j <= intervals[i].second; j++) {
                    res2[j] = collect_res[j];
                }
            }

//           for (auto i : x){
//               std::cout << i << " ";
//           }
//            std::cout << std::endl;

            error = norm(vector_diff(res2, b)) / norm(b);

            if (error < eps) {
                std::cout << "Result: " << std::endl;
//                for (double re : x){                  ВЫВОД УБРАН, ТАК КАК ОН ЛОМАЕТ ПРОЦЕССЫ НА БОЛЬШИХ МАТРИЦАХ
//                    std::cout << re << " ";
//                }
                std::cout << std::endl;
                fflush(stdout);

                MPI_Finalize();
                break;
            }
        }
    }
    for (int k = 0; k < 30; k++) {
        if (rank != 0) {
            MPI_Recv(&send_res[0], n, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            mult(A, send_res, intervals[rank].first, intervals[rank].second, collect_res);
            MPI_Send(&collect_res[0], n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    }
    auto finish = std::chrono::high_resolution_clock::now();

    auto duration = duration_cast<std::chrono::microseconds>(finish - go);
    std::cout << "Time: "<< duration.count() << " ms" << std::endl;

    std::cout << "Result: " << std::endl;
    if (rank == 0)
        for (double re : x){
            std::cout << re << " ";
        }
}
/* 1: 9.2
 * 2: 8.2
 * 3: 7.8
 * 4: 7.6
 * 5: 7.5
 */