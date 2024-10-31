#include <iostream>
#include <random>
#include <chrono>
#include <omp.h>

#define MAX_STEPS 100000

int main() {
    int a, b, x, x_tmp, N, i, hits_b = 0, stepCount, sum_stepCount = 0, NUM_THREADS;
    double p;
    
    std::cout << "Введите начальную и конечную точку" << std::endl;
    std::cin >> a >> b;
    std::cout << "Введите исходную точку и количество частиц" << std::endl;
    std::cin >> x >> N;
    std::cout << "Введите вероятность перехода" << std::endl;
    std::cin >> p;
    std::cout << "Введите количество нитей" << std::endl;
    std::cin >> NUM_THREADS;

    int glob_time = time(NULL);

    #pragma omp parallel num_threads(NUM_THREADS) default(none) shared(N, a, b, x, p, glob_time) reduction(+:hits_b) reduction(+:sum_stepCount)
    {
        int my_id = omp_get_thread_num();
        std::mt19937 my_gen(glob_time * my_id);
        std::uniform_real_distribution<double> dis(0.0, 1.0);
        
        #pragma omp for private(x_tmp, stepCount)
        for (i = 0; i < N; ++i) {
            x_tmp = x;
            stepCount = 0;
            while (x_tmp != a && x_tmp != b && stepCount <= MAX_STEPS) {
                stepCount++;
                x_tmp += dis(my_gen) > p ? 1 : -1;
            }
            if (x_tmp == b) hits_b++;
            sum_stepCount += stepCount;
        }
    }

    std::cout << "Вероятность достижения b: " << (double) hits_b / (double) N << std::endl;
    std::cout << "Среднее время жизни частицы: " << sum_stepCount / N << std::endl;

    return 0;
}
