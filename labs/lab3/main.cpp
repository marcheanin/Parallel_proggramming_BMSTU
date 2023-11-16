#include <omp.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <random>
#include <chrono>
#include <ratio>

double norm(const std::vector <double>& v){
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
    int i, j;
    #pragma omp parallel for shared(A, x, res, raw1, raw2) private(i, j)
    for (i = raw1; i <= raw2; i++){
        res[i] = 0;
        for (j = 0; j < x.size(); j++){
            res[i] += A[i][j] * x[j];
        }
    }
}

void default_mult(std::vector <std::vector <double> > A, std::vector <double> x, int raw1, int raw2, std::vector <double>& res){
    double sum;
    int i, j;
    for (i = raw1; i <= raw2; i++){
        sum = 0;
        for (j = 0; j < x.size(); j++){
            sum += A[i][j] * x[j];
        }
        res[i] = sum;
    }
}

std::vector <double> mult_digit_vector(double a, std::vector <double> b) {
    int i = 0;
    #pragma omp parallel for shared(a, b) private(i)
    for (int i = 0; i < b.size(); i++){
        b[i] *= a;
    }
    return b;
}

std::vector <double> mult_digit_vector_default(double a, std::vector <double> b) {
    int i = 0;
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
    std::vector <double> x(m, 10);
    //std::vector <double> b(m, n+1);
    std::vector <double> b (n);
    for (double & i : b){
        i = n+1;
    }
    std::vector <double> res(m, 0);
    std::vector <double> Ay(m, 0);

    double tau, error;
    int iter = 1;
    while(true){
        mult(A, x, 0, n-1, res);

        auto y = vector_diff(res, b);
        mult(A, y, 0, n-1, Ay);
        tau = vector_mult(y, Ay) / vector_mult(Ay, Ay);
        x = vector_diff(x, mult_digit_vector(tau, y));

        std::cout << "tau: " << tau << std::endl;

        mult(A, x, 0, n-1, res);

        assert(norm(b) != 0);

        error = norm(vector_diff(res, b)) / norm(b);
        std::cout << iter << " " << error << std::endl;
        iter++;

        if (error < eps) break;
    }

    return x;
}

std::vector <double> solve_without_parallel(std::vector <std::vector <double> > A, double eps){
    int n = A.size(), m = A[0].size();
    std::vector <double> x(m, 10);
    //std::vector <double> b(m, n+1);
    std::vector <double> b (n);
    for (double & i : b){
        i = n+1;
    }
    std::vector <double> res(m, 0);
    std::vector <double> Ay(m, 0);

    double tau, error;
    int iter = 1;
    while(true){
        default_mult(A, x, 0, n-1, res);

        auto y = vector_diff(res, b);
        default_mult(A, y, 0, n-1, Ay);
        tau = vector_mult(y, Ay) / vector_mult(Ay, Ay);
        x = vector_diff(x, mult_digit_vector_default(tau, y));

        std::cout << "tau: " << tau << std::endl;

        default_mult(A, x, 0, n-1, res);

        assert(norm(b) != 0);

        error = norm(vector_diff(res, b)) / norm(b);
        std::cout << iter << " " << error << std::endl;
        iter++;

        if (error < eps) break;
    }

    return x;
}

int main(int argc, char** argv) {
    int matrix_sixe = 10000;
    std::vector <std::vector <double> > a(matrix_sixe, std::vector <double> (matrix_sixe));
    for (int i = 0; i < a.size(); i++){
        for (int j = 0; j < a[0].size(); j++){
            if (i == j) a[i][j] = 2.0;
            else a[i][j] = 1.0;
        }
    }

    auto go = std::chrono::high_resolution_clock::now();
    auto res = solve(a, 0.00001);
    auto finish = std::chrono::high_resolution_clock::now();

    auto duration =  std::chrono::duration_cast <std::chrono::microseconds> (finish - go);
    std::cout << duration.count() << " ms" << std::endl;

    go = std::chrono::high_resolution_clock::now();
    auto res1 = solve_without_parallel(a, 0.00001);
    finish = std::chrono::high_resolution_clock::now();

    duration =  std::chrono::duration_cast <std::chrono::microseconds> (finish - go);
    std::cout << duration.count() << " ms" << std::endl;

}