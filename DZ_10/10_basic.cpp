#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>

#define N 128 // Сторона куба
#define ITER 100 // Количество итераций метода Якоби

void go_random(std::vector<double>& grid, int nx, int ny, int nz) {
    for (int i = 0; i < nx * ny * nz; ++i) {
        grid[i] = static_cast<double>(rand()) / RAND_MAX;
    }
}

double go_norm(const std::vector<double>& grid, const std::vector<double>& newGrid, int nx, int ny, int nz) {
    double localNorm = 0.0;
    for (int i = 0; i < nx * ny * nz; i++) {
        double diff = newGrid[i] - grid[i];
        localNorm += diff * diff;
    }
    return sqrt(localNorm);
}

MPI_Datatype createBoundaryType(int nx, int ny, int nz, int direction) {
    MPI_Datatype boundaryType;
    switch(direction) {
        case 0:
            MPI_Type_vector(ny * nz, 1, nx, MPI_DOUBLE, &boundaryType);
            break;
        case 1:
            MPI_Type_vector(nz, nx, nx * ny, MPI_DOUBLE, &boundaryType);
            break;
        case 2:
            MPI_Type_contiguous(nx * ny, MPI_DOUBLE, &boundaryType); //срез элементов подряд + эффективность
            break;
    }
    MPI_Type_commit(&boundaryType);
    return boundaryType;
}

void jacobiIteration(std::vector<double>& grid, std::vector<double>& newGrid, int nx, int ny, int nz) {
    for (int k = 1; k < nz - 1; k++)
        for (int j = 1; j < ny - 1; j++)
            for (int i = 1; i < nx - 1; i++) {
                int idx = i + j * nx + k * nx * ny;                
                newGrid[idx] = (
                    (2 * grid[idx] - grid[idx-1] - grid[idx+1]) / (nx * nx) +
                    (2 * grid[idx] - grid[idx-nx] - grid[idx+nx]) / (ny * ny) +
                    (2 * grid[idx] - grid[idx-nx*ny] - grid[idx+nx*ny]) / (nz * nz)
                );
            }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Создание декартовой топологии
    int dims[3] = {0, 0, 0};
    MPI_Dims_create(size, 3, dims);
    int periods[3] = {0, 0, 0};
    MPI_Comm cartComm;
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 1, &cartComm);

    // Размер локальной области
    int nx = N / dims[0], ny = N / dims[1], nz = N / dims[2];

    // Инициализация сетки
    std::vector<double> grid(nx * ny * nz);
    std::vector<double> newGrid(nx * ny * nz);
    go_random(grid, nx, ny, nz);

    // Создание типов для обмена
    MPI_Datatype xBoundary = createBoundaryType(nx, ny, nz, 0);
    MPI_Datatype yBoundary = createBoundaryType(nx, ny, nz, 1);
    MPI_Datatype zBoundary = createBoundaryType(nx, ny, nz, 2);

    int left, right, up, down, front, back;
    MPI_Cart_shift(cartComm, 0, 1, &left, &right);
    MPI_Cart_shift(cartComm, 1, 1, &down, &up);
    MPI_Cart_shift(cartComm, 2, 1, &front, &back);

    // Итерации метода Якоби
    double globalNorm;

    double start_time = MPI_Wtime();
    for (int iter = 0; iter < ITER; ++iter) {
        // Обмен границ
        MPI_Request requests[6];
        MPI_Isend(grid.data(), 1, xBoundary, left, 0, cartComm, &requests[0]);
        MPI_Irecv(grid.data(), 1, xBoundary, right, 0, cartComm, &requests[1]);
        MPI_Isend(grid.data(), 1, yBoundary, down, 0, cartComm, &requests[2]);
        MPI_Irecv(grid.data(), 1, yBoundary, up, 0, cartComm, &requests[3]);
        MPI_Isend(grid.data(), 1, zBoundary, front, 0, cartComm, &requests[4]);
        MPI_Irecv(grid.data(), 1, zBoundary, back, 0, cartComm, &requests[5]);

        MPI_Waitall(6, requests, MPI_STATUSES_IGNORE);

        // Выполнение итерации
        jacobiIteration(grid, newGrid, nx, ny, nz);

        // Вычисление нормы на последней итерации
        if (iter == ITER - 1) {
            double localNorm = go_norm(grid, newGrid, nx, ny, nz);
            MPI_Reduce(&localNorm, &globalNorm, 1, MPI_DOUBLE, MPI_SUM, 0, cartComm);
        }

        // Обновление сетки
        grid.swap(newGrid);
    }
    double end_time = MPI_Wtime();
    double elapsed_time = end_time - start_time;
    double max_time;
    MPI_Reduce(&elapsed_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "Глобальная норма: " << globalNorm << std::endl;
        std::cout << "Время выполнения: " << max_time << " секунд" << std::endl;
    }

    // Освобождение производных типов
    MPI_Type_free(&xBoundary);
    MPI_Type_free(&yBoundary);
    MPI_Type_free(&zBoundary);

    MPI_Comm_free(&cartComm);
    MPI_Finalize();
    return 0;
}