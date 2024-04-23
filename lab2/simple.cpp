
#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <functional>
#include <cmath>

struct Task {
    double start;
    double end;
    double result;
    uint32_t steps;
};

class IntegralSolver {
public:
    IntegralSolver(double a, double b, uint32_t totalSteps, int tasksNumber):threadCount_(tasksNumber) {
        tasks_ = new Task[tasksNumber];
        for (int i = 0; i < tasksNumber; ++i) {
            tasks_[i].start = a + i * (b - a) / tasksNumber;
            tasks_[i].end = a + (i + 1) * (b - a) / tasksNumber;
            tasks_[i].steps = totalSteps / tasksNumber;
        }
    }

    void start(void) {
        std::vector<std::thread> threads;
        auto start_time = std::chrono::steady_clock::now();

        for (int i = 0; i < threadCount_; ++i) {
            threads.emplace_back([this, i]() { this->process(&this->tasks_[i]); });
        }

        for (auto& th : threads){
            th.join();
        }

        auto stop_time = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds> (stop_time - start_time).count();
        std::cout << "Elapsed: " << elapsed / 1000. << " sec \n";

        for (int i = 0; i < threadCount_; ++i)
        {
            result += tasks_[i].result;
        }
    }

    double getResult() const { return result; }

private:
    double result = 0.0;
    Task* tasks_;
    int threadCount_ = 0;
  
    void process(Task* data) {
        double range = data->end - data->start;
        double step = range / data->steps;
        double integral = 0.0;
        for(int i = 0; i < data->steps; ++i) {
            double x = data->start + (i + 0.5) * step;
            integral += sin(1/x);
        }
        data->result = integral * step;
    }
};

int main(void)
{
    const int threadsCount = 4;
    IntegralSolver solver(0.001, 1.0, 100000000, threadsCount);  
    solver.start();

    std::cout << "Integral result is: " << solver.getResult() << std::endl;
}
