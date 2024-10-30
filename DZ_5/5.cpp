#include <iostream>
#include <random>
#include <chrono>
#include <omp.h>

#define MAX_STEPS 10000

int main() {
    int a, b, x, x_tmp, N, i, hits_b = 0, stepCount, NUM_THREADS;
    double p;
    
    std::cout << "Введите начальную и конечную точку" << std::endl;
    std::cin >> a >> b;
    std::cout << "Введите исходную точку и количество частиц" << std::endl;
    std::cin >> x >> N;
    std::cout << "Введите вероятность перехода" << std::endl;
    std::cin >> p;
    std::cout << "Введите количество нитей" << std::endl;
    std::cin >> NUM_THREADS;

    auto start = std::chrono::system_clock::now();

    #pragma omp parallel num_threads(NUM_THREADS) default(none) shared(N, a, b, x, p, NUM_THREADS) reduction(+:hits_b)
    {
        #pragma omp for private(x_tmp, stepCount)
        for (i = 0; i < N; ++i) {
            std::random_device generator;
            std::uniform_real_distribution<double> distribution(0.0, 1.0);
            x_tmp = x;
            stepCount = 0;
            while (x_tmp != a && x_tmp != b && stepCount <= MAX_STEPS) {
                stepCount++;
                if (distribution(generator) > p)
                    x_tmp += 1;
                else
                    x_tmp -= 1;
            }
            if (x_tmp == b)
                hits_b++;
        }
    }

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;

    std::cout << "Probability of reaching b: " << (double) hits_b / (double) N << std::endl;
    std::cout << "Basic time: " << elapsed.count() << " ms" << std::endl;
    std::cout << "Time per dot: " << elapsed.count() / N << " ms" << std::endl;

    delete[] go_random;
    return 0;
}
