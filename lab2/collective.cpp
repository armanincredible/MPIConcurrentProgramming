#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <functional>
#include <cmath>
#include <queue>

struct Task {
    double start;
    double end;
    uint32_t steps;
};

class IntegralSolver {
public:
    IntegralSolver(double a, double b, uint32_t totalSteps, uint32_t tasks) {

        for (int i = 0; i < tasks; ++i) {
            double taskStart = a + (b - a) / tasks * i;
            double taskEnd = taskStart + (b - a) / tasks;
            queue.push({taskStart, taskEnd, totalSteps / tasks});
        }
    }

    void start(int threadCount) {
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

    double getResult() const { return result; }

private:
    double result = 0.0;
    std::queue<Task> queue;
    std::mutex mutex;
  
    void process() {
        while (true) {
            mutex.lock();
            if (queue.empty()) {
                mutex.unlock();
                return; 
            }
            Task task = queue.front();
            queue.pop();
            mutex.unlock();

            double localSum = 0.0;
            double stepSize = (task.end - task.start) / task.steps;
            for (int i = 0; i < task.steps; ++i) {
                double x = task.start + (i + 0.5) * stepSize;
                localSum += std::sin(1/x) * stepSize;
            }

            mutex.lock();
            result += localSum;
            mutex.unlock();
        }
    }
};

int main(void)
{
    IntegralSolver solver(0.001, 1.0, 1410065400, 100);  
    solver.start(4); 

    std::cout << "Integral result is: " << solver.getResult() << std::endl;
}
