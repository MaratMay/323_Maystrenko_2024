#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>

#define N 1024
#define N_ITER 1000

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Создаем ленточную подобласть 💕
    int local_N = N / size;
    std::vector<double> local_grid(local_N * N);
    std::vector<double> new_local_grid(local_N * N);

    // Инициализируем сетку случайными значениями 🌟
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    for (int i=0; i < local_N * N; i++) local_grid[i] = dis(gen);

    for (int iter = 0; iter < N_ITER; iter++) {
        // Обмениваемся границами с соседями 💌
        if (rank > 0) { 
            MPI_Send(&local_grid[0], N, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&local_grid[local_N * N], N, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (rank < size - 1) {
            MPI_Send(&local_grid[(local_N - 1) * N], N, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&local_grid[-N], N, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        for (int i=1; i < local_N - 1; i++)
            for (int j=1; j < N - 1; j++)
                new_local_grid[i * N + j] = 0.25 * (
                    local_grid[(i - 1) * N + j] +
                    local_grid[(i + 1) * N + j] +
                    local_grid[i * N + (j - 1)] +
                    local_grid[i * N + (j + 1)]);
        std::swap(local_grid, new_local_grid);
    }

    double local_norm = 0.0;
    for (int i=1; i < local_N - 1; i++) {
        for (int j=1; j < N - 1; j++) {
            double diff = local_grid[i*N + j] - new_local_grid[i*N + j];
            local_norm += diff * diff;
        }
    }

    double global_norm;
    MPI_Reduce(&local_norm, &global_norm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (!rank) std::cout << "Итоговая норма: " << std::sqrt(global_norm) << std::endl;

    MPI_Finalize();
    return 0;
}
