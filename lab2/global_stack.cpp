#include <cstdio>
#include <vector>
#include <cmath>
#include <mutex>
#include <cassert>
#include <thread>
#include <array>
#include <iostream>
#include <sstream>


struct Task_t
{
    double start = 0;
    double end = 0;
    double fn_start = 0;
    double fn_end = 0;
    double area = 0;
};

double 
Func( double x)
{
    return std::sin(1/x);
}


const double kEpsilon = 1e-7;
const size_t kMaxLocalStackSize  = 1;
const size_t kMaxGlobalStackSize = 1000;


class IntegralSolver {
public:
    IntegralSolver(double start, double end) {
        double fn_start = Func( start);
        double fn_end   = Func( end);
        double area     = (fn_end + fn_start) * (end - start) / 2;

        gStack.push_back( Task_t{ start, end, fn_start, fn_end, area});
    }

    void start(int threadCount) {
        kThreadsCount = threadCount;
        std::vector<std::thread> threads;

        auto start_time = std::chrono::steady_clock::now();
        for (int i = 0; i < threadCount; ++i) {
            threads.emplace_back([this]() { this->process(); });
        }

        for (auto& th : threads) th.join();

        auto stop_time = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds> (stop_time - start_time).count();
        std::cout << "Elapsed: " << elapsed / 1000. << " sec \n";
    }

    double getResult() const { return gResult; }

private:
    std::mutex gStackMutex;
    std::mutex gHasTaskMutex;
    std::vector<Task_t> gStack;

    std::mutex gResultMutex;
    double gResult = 0;
    size_t gActiveThreads = 0;

    size_t kThreadsCount = 0;
  
    void process() {
        Task_t task;

        for ( ;; )
        {
            {
                gHasTaskMutex.lock();

                std::lock_guard<std::mutex> scoped_lock{ gStackMutex};

                assert( !gStack.empty() );

                {
                    task = gStack.back();
                    gStack.pop_back();
                }

                if ( !gStack.empty() )
                {
                    gHasTaskMutex.unlock();
                }

                if ( task.start > task.end)
                {
                    break;
                }

                gActiveThreads++;
            }

            {
                double result = 0;
                std::vector<Task_t> stack;

                for ( ;; )
                {
                    double mid = (task.start + task.end) / 2;
                    double fn_mid = Func( mid);

                    double area_start_mid = (fn_mid + task.fn_start) * (mid - task.start) / 2;
                    double area_mid_end   = (task.fn_end + fn_mid) * (task.end - mid) / 2;
                    
                    double new_area = area_start_mid + area_mid_end;

                    if ( std::abs( new_area - task.area) >= kEpsilon * std::abs( new_area) )
                    {
                        stack.push_back( Task_t{ task.start, mid, task.fn_start, fn_mid, area_start_mid});
                        task.start = mid;
                        task.fn_start = fn_mid;
                        task.area = area_mid_end;
                    } else
                    {
                        result += new_area;

                        if ( stack.empty() )
                        {
                            break;
                        }

                        {
                            task = stack.back();
                            stack.pop_back();
                        }
                    }

                    // Move tasks to the global stack
                    if ( stack.size() >= kMaxLocalStackSize )
                    {
                        std::lock_guard<std::mutex> scoped_lock{ gStackMutex};

                        if ( gStack.empty() )
                        {
                            size_t n_elements = std::min( stack.size(), kMaxGlobalStackSize);
                            auto from = stack.end() - n_elements;
                            auto to   = stack.end();
                            gStack.insert( gStack.end(), make_move_iterator( from), make_move_iterator( to));
                            stack.erase( from, to);

                            gHasTaskMutex.unlock();
                        }
                    }
                }

                {
                    std::lock_guard<std::mutex> scoped_lock{ gResultMutex};
                    gResult += result;
                }
            }

            {
                std::lock_guard<std::mutex> scoped_lock{ gStackMutex};
                
                gActiveThreads--;

                if ( gActiveThreads == 0 
                    && gStack.empty() )
                {
                    for ( size_t i = 0; i != kThreadsCount; ++i )
                    {
                        gStack.push_back( Task_t{ 1, 0, 0, 0, 0});
                    }

                    gHasTaskMutex.unlock();
                }
            }
        }
    }
};


int 
main()
{
    const int threadsCount = 4;
    IntegralSolver solver(0.001, 1.0);
    solver.start(threadsCount); 

    std::cout << "Integral result is: " << solver.getResult() << std::endl;

    return 0;
}

