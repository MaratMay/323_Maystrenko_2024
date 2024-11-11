#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <string>

int main(int argc, char** argv) {
    if (argc != 2) return 1;

    int N = 1024;
    int n_iter = strtol(argv[1], NULL, 10);

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        std::cout << "Размер сетки: " << N << "x" << N << std::endl;
        std::cout << "Количество итераций: " << n_iter << std::endl;
    }

    // Создаем ленточную подобласть 💕
    int local_N = N / size;
    std::vector<double> local_grid(local_N * N);
    std::vector<double> new_local_grid(local_N * N);
    std::vector<double> first_row(N - 2);
    std::vector<double> last_row(N - 2);

    // Инициализируем сетку случайными значениями 🌟
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    
    for (int i = 0; i < local_N * N; ++i) {
        local_grid[i] = dis(gen);
        new_local_grid[i] = local_grid[i];
    }

    double start_time = MPI_Wtime();
    for (int iter = 0; iter < n_iter; iter++) {
        // Обмениваемся границами с соседями 💌
        if (rank > 0) { 
            MPI_Sendrecv(local_grid.data() + N, N - 2, MPI_DOUBLE, rank - 1, 0, first_row.data(), N - 2, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (rank < size - 1) {
            MPI_Sendrecv(local_grid.data() + (local_N - 1) * N + 1, N - 2, MPI_DOUBLE, rank + 1, 0, last_row.data(), N - 2, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // Реализуем алгоритм 💀
        for (int i=1; i < local_N - 1; i++)
            for (int j=1; j < N - 1; j++)
                new_local_grid[i * N + j] = 0.25 * (
                    local_grid[(i - 1) * N + j] +
                    local_grid[(i + 1) * N + j] +
                    local_grid[i * N + (j - 1)] +
                    local_grid[i * N + (j + 1)]);

        if (rank > 0)
            for (int i=1; i < N-1; i++)
                new_local_grid[i] = 0.25 * (
                    local_grid[i - 1] + 
                    local_grid[i + 1] + 
                    first_row[i - 1] + 
                    first_row[i + 1]);

        if (rank < size - 1)
            for (int i=1; i < N-1; i++)
                new_local_grid[(local_N - 1) * N + i] = 0.25 * (
                    local_grid[(local_N - 1) * N + i - 1] +
                    local_grid[(local_N - 1) * N + i + 1] + 
                    last_row[i - 1] + 
                    last_row[i + 1]
                );
    }

    double local_norm = 0.0;
    for (int i=0; i < size; i++) local_norm += (local_grid[i] - new_local_grid[i]) * (local_grid[i] - new_local_grid[i]);

    local_grid.swap(new_local_grid);

    double global_norm;
    MPI_Reduce(&local_norm, &global_norm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    double end_time = MPI_Wtime();
    double elapsed_time = end_time - start_time;
    double max_time;
    MPI_Reduce(&elapsed_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        std::cout << "Итоговая норма: " << std::sqrt(global_norm) << std::endl;
        std::cout << "Время выполнения: " << max_time << " секунд" << std::endl;
    }

    MPI_Finalize();
    return 0;
}
