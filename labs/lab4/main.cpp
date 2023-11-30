#include <iostream>
#include <thread>
#include <mutex>
#include <vector>
#include <random>
#include <condition_variable>

const int num_of_philosophers = 5;

std::mutex print_lock;

std::atomic <bool> meal_flag {true};

//class Waiter{
//    std::mutex mutex;
//    std::condition_variable cv;
//};

class Fork{
public:
    std::mutex mutex;
    std::atomic <bool> locked = false;
    void take_fork() {
        while(locked){};
        locked = true;
    }
    void free(){
        locked = false;
    }
};

class Table{
public:
    std::atomic <bool> ready {false};
    std::array <Fork, num_of_philosophers> forks;
};

class Philosopher {
private:
    Fork &lfork;
    Fork &rfork;
    int num{};
    Table &table;
    std::thread phthread;
    std::mt19937 rd { std::random_device{}() };
    std::string status;


public:
    int num_finished_meals = 0;

    Philosopher(int _num, Table &_table, Fork &l, Fork &r) :
        num(_num),
        table(_table),
        lfork(l),
        rfork(r),
        phthread(&Philosopher::process, this){}

    ~Philosopher(){
        std::cout << "Killing philosopher " << std::to_string(num) << std::endl;
        phthread.join();
    }

    void process(){
        while(!table.ready){}
        while(table.ready && meal_flag) {
            think();
            if (!meal_flag) return;
            eat();
        }
    }

    void output(const std::string& text) const {
        std::lock_guard<std::mutex> cout_lock(print_lock);
        std::cout << "philosopher № " << num << text << std::endl;
    }

    void eat() {
        //while(lfork.locked);
        //std::lock_guard<std::mutex> left_fork_lock(lfork.mutex, std::adopt_lock);
        lfork.take_fork();
        //lfork.locked = true;

        status = "tacked left fork";

        thread_local std::uniform_int_distribution<> forking(1, 3);
        std::this_thread::sleep_for(std::chrono::milliseconds(forking(rd) * 100));

        //while(rfork.locked);
        //std::lock_guard<std::mutex> right_fork_lock(rfork.mutex, std::adopt_lock);
        rfork.take_fork();
        //lfork.locked = true;

        //output(" is eating");
        status = "is eating";

        thread_local std::uniform_int_distribution <> dist (1, 3);
        std::this_thread::sleep_for(std::chrono::milliseconds(dist(rd) * 100));

        //output(" finished eating.");
        status = "finished eating";
        lfork.free();
        rfork.free();
        num_finished_meals++;
    }

    void think() {
        thread_local std::uniform_int_distribution <> dist (1, 3);
        std::this_thread::sleep_for(std::chrono::milliseconds(dist(rd) * 100));

        //output(" is thinking.");
        status = "is thinking";
    }

    std::string get_status(){
        return std::to_string(num) + ": " + status + " ";
    }

};


int main() {
    std::cout << "Start dinner for " << num_of_philosophers << " guests" << std::endl;

    Table table;
    std::array <Philosopher, num_of_philosophers> philosophers {
            {
                    {1, table, table.forks[0], table.forks[1]},
                    {2, table, table.forks[1], table.forks[2]},
                    {3, table, table.forks[2], table.forks[3]},
                    {4, table, table.forks[3], table.forks[4]},
                    {5, table, table.forks[0], table.forks[4]}
            }};

    table.ready = true;
    //std::this_thread::sleep_for(std::chrono::seconds(15));

    for (int t = 0; t < 100; t++){
        for (auto & philosopher : philosophers){
            std::cout << philosopher.get_status();
        }
        std::cout << std::endl;
//        for (auto & fork : table.forks) {
//            std::cout << fork.mutex.try_lock() << " ";
//        }
//        std::cout << std::endl;
        std::this_thread::sleep_for(std::chrono::milliseconds(250));
    }

    meal_flag = false;
    std::cout << "Dinner done." << std::endl;
    for (auto & philosopher : philosophers){
        std::cout << philosopher.num_finished_meals << std::endl;
    }

    return 0;
}
