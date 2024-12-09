#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <string>

#define N 64

void printMatrix(const std::vector<int>& matrix) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            std::cout << matrix[i * N + j] << " ";
        }
        std::cout << std::endl;
    }
}

std::vector<int> generateRandomMatrix(int key) {
    std::mt19937 generator(key);
    std::uniform_int_distribution<> distribution(1, 10);
    std::vector<int> result(N*N);
    for (int i = 0; i < N*N; i++) result[i] = distribution(generator);
    return result;
}

std::vector<int> multiply(std::vector<int> &A, std::vector<int> &B) {
    std::vector<int> result(N * N, 0);
    for (int i=0; i < N; i++)
        for (int j=0; j < N; j++)
            for (int k=0; k < N; k++)
                result[i*N + j] += A[i*N + k] * B[k*N + j];
    return result;
}

void add_to_me(std::vector<int> &A, std::vector<int> B) {
    for (int i=0; i < N*N; i++) A[i] += B[i];
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int p = sqrt(size);

    if (rank == 0) {
        std::cout << "Работаем с матрицами размера " << N * p << " на " << N * p << std::endl;
        std::cout << "Каждому процессу (" << size << " штук) выделяется матрица размера " << N << " на " << N << std::endl;
    }

    auto A = generateRandomMatrix(rank);
    auto B = generateRandomMatrix(rank + size);
    std::vector<int> C(N*N, 0), tmp_1(N*N), tmp_2(N*N);

    double start_time = MPI_Wtime();

    int dims[] = {p, p};
    int periods[] = {0, 0};
    int reorder = 0;
    MPI_Comm Square_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &Square_comm);

    int my_coords[2];
    MPI_Cart_coords(Square_comm, rank, 2, my_coords);

    int left, right, up, down, request_iter=0;
    MPI_Request requests[p * 2];

    for(int i=1; i < p; i++) {
        MPI_Cart_shift(Square_comm, 0, i, &left, &right);
        if (left != MPI_PROC_NULL) MPI_Isend(B.data(), N * N, MPI_INT, left, 0, Square_comm, &requests[request_iter++]);
        if (right != MPI_PROC_NULL) MPI_Isend(B.data(), N * N, MPI_INT, right, 0, Square_comm, &requests[request_iter++]);

        MPI_Cart_shift(Square_comm, 1, i, &up, &down);
        if (up != MPI_PROC_NULL) MPI_Isend(A.data(), N * N, MPI_INT, up, 0, Square_comm, &requests[request_iter++]);
        if (down != MPI_PROC_NULL) MPI_Isend(A.data(), N * N, MPI_INT, down, 0, Square_comm, &requests[request_iter++]);
    }

    int src_A;
    int src_B;

    for(int i=0; i < p; i++) {
        src_A = my_coords[0] * p + i;
        src_B = my_coords[1] + i * p;

        if (src_A == rank) tmp_1 = A; 
        else MPI_Recv(tmp_1.data(), N * N, MPI_INT, src_A, MPI_ANY_TAG, Square_comm, MPI_STATUS_IGNORE);

        if (src_B == rank) tmp_2 = B;
        else MPI_Recv(tmp_2.data(), N * N, MPI_INT, src_B, MPI_ANY_TAG, Square_comm, MPI_STATUS_IGNORE);
        
        add_to_me(C, multiply(tmp_1, tmp_2));
    }

    MPI_Waitall(request_iter, requests, MPI_STATUSES_IGNORE);
    double end_time = MPI_Wtime();
    double elapsed_time = end_time - start_time;
    double max_time;

    MPI_Barrier(Square_comm);
    MPI_Reduce(&elapsed_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, Square_comm);

    std::vector<int> main_A, main_B, main_C;

    if (rank == 0) {
        std::cout << "Время выполнения: " << max_time << " секунд" << std::endl;
        std::cout << "0 A:------------------" << std::endl;
        printMatrix(A);

        for (int i=1; i < size; i++) {
            MPI_Recv(tmp_1.data(), N * N, MPI_INT, i, MPI_ANY_TAG, Square_comm, MPI_STATUS_IGNORE);
            std::cout << i << " A:------------------" << std::endl;
            printMatrix(tmp_1); }

        std::cout << "0 B:------------------" << std::endl;
        printMatrix(B);

        for (int i=1; i < size; i++) {
            MPI_Recv(tmp_1.data(), N * N, MPI_INT, i, MPI_ANY_TAG, Square_comm, MPI_STATUS_IGNORE);
            std::cout << i << " B:------------------" << std::endl;
            printMatrix(tmp_1);}

        std::cout << "0 C:------------------" << std::endl;
        printMatrix(C);

        for (int i=1; i < size; i++) {
            MPI_Recv(tmp_1.data(), N * N, MPI_INT, i, MPI_ANY_TAG, Square_comm, MPI_STATUS_IGNORE);
            std::cout << i << " C:------------------" << std::endl;
            printMatrix(tmp_1); }

    } else {
        MPI_Send(A.data(), N * N, MPI_INT, 0, 0, Square_comm);
        MPI_Send(B.data(), N * N, MPI_INT, 0, 0, Square_comm);
        MPI_Send(C.data(), N * N, MPI_INT, 0, 0, Square_comm);
    }

    MPI_Comm_free(&Square_comm);
    MPI_Finalize();
    return 0;
}