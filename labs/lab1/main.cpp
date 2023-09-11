#include <iostream>
#include <random>
#include <thread>
#include <chrono>
#include <vector>
#include <ctime>
#include <mutex>

std::mutex mtx;

void mult(const std::vector <std::vector <int> >& a,
          const std::vector <std::vector <int> >& b,
          std::vector <std::vector <int> >& c,
          int x1, int y1, int x2, int y2){

    int t = 0;
    for (int i = x1; i <= x2; i++){
        for (int j = y1; j <= y2; j++){
            for (int k = 0; k < a.size(); k++){
                t += a[i][k] * b[k][j];
            }
            c[i][j] = t;
            t = 0;
        }
    }

}

void fillMatrixWithRandom(std::vector <std::vector <int> >& mt) {
    std::random_device rd;
    std::mt19937 gen(rd());
    for (int i = 0; i < mt.size(); i++){
        for (int j = 0; j < mt[0].size(); j++){
            mt[i][j] = 1 + gen() % 100;
        }
    }
}

void printMatrix(const std::vector <std::vector <int> >& mt){
    for (int i = 0; i < mt.size(); i++){
        for (int j = 0; j < mt[0].size(); j++){
            std::cout << mt[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

bool compareMatrix(const std::vector <std::vector <int> >& a, const std::vector <std::vector <int> >& b) {
    for (int i = 0; i < a.size(); i++){
        for (int j = 0; j < a[0].size(); j++){
            if (a[i][j] != b[i][j])
                return false;
        }
    }
    return true;
}

int main() {

    const int matrix_size = 1000;

    std::vector <std::vector <int> > a(matrix_size, std::vector <int> (matrix_size));
    std::vector <std::vector <int> > b(matrix_size, std::vector <int> (matrix_size));
    std::vector <std::vector <int> > c(matrix_size, std::vector <int> (matrix_size));
    std::vector <std::vector <int> > d(matrix_size, std::vector <int> (matrix_size));

    fillMatrixWithRandom(a);
    fillMatrixWithRandom(b);

    //printMatrix(a);
    //printMatrix(b);

    auto go = std::chrono::high_resolution_clock::now();
    mult(a, b, c, 0, 0, a.size()-1, a.size()-1);
    auto finish = std::chrono::high_resolution_clock::now();

    auto duration = duration_cast<std::chrono::microseconds>(finish - go);
    std::cout << "Mult without threads: "<< duration.count() << " ms" << std::endl;
    //printMatrix(c);
    std::cout << std::endl;

    std::vector <std::thread> thvec;
    int max_th_count = 10;

    for (int th_num = 1; th_num <= max_th_count; th_num++) {
        go = std::chrono::high_resolution_clock::now();
        int chunk_size = matrix_size / th_num, bonus = matrix_size - chunk_size * th_num;

        for (int start = 0, end = chunk_size; start < matrix_size; start = end, end = start + chunk_size){
            if (bonus) {
                end++;
                bonus--;
            }
            thvec.emplace_back(mult, std::ref(a), std::ref(b), std::ref(d), start, 0,
                              end - 1, matrix_size - 1);
        }
        for (auto &t : thvec) {
            t.join();
        }
        finish = std::chrono::high_resolution_clock::now();
        std::cout << "Equal? " << compareMatrix(c, d) << std::endl;
        thvec.clear();
        for (auto &v: d) {
            std::fill(v.begin(), v.end(), 0);
        }

        duration = duration_cast <std::chrono::microseconds> (finish - go);
        std::cout << "Mult with " << th_num <<  " threads: "<< duration.count() << " ms" << std::endl;

    }


    //printMatrix(c);


    return 0;
}
