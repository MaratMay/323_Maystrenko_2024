#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <string>

#define N 1024
#define MAX_ITER 10000

void fillRandomMainField(std::vector<char>& main_field, int& n_alive_now, double fill_percentage = 0.2) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    n_alive_now = 0;
    for (int i = 0; i < N * N; ++i) {
        if (dis(gen) < fill_percentage) {
            main_field[i] = 1;
            n_alive_now++;
        } else {
            main_field[i] = 0;
        }
    }
}


int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n_iter = 0, n_alive_now = 0, n_alive_past = -1, n_alive_now_per_thread;

    std::vector<char> main_field(N * N, 0);
    if (rank == 0) {
        fillRandomMainField(main_field, n_alive_now);
        std::cout << "Поле " << N << " x " << N << " случайно заполнено! Начальное количество живых клеток: " << n_alive_now << std::endl;
    }

    int local_N = N / size;
    char tmp_sum;
    std::vector<char> local_grid(local_N * N);
    std::vector<char> new_local_grid(local_N * N);
    std::vector<char> next_first_row(N);
    std::vector<char> prev_last_row(N);

    MPI_Scatter(main_field.data(), N * local_N, MPI_CHAR, local_grid.data(), N * local_N, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Request request_recv_1, request_recv_2, request_send_1, request_send_2;
    MPI_Bcast(&n_alive_now, 1, MPI_INT, 0, MPI_COMM_WORLD);

    double start_time = MPI_Wtime();
    while (n_alive_now != n_alive_past && n_iter < MAX_ITER) {
        n_iter++;
        n_alive_now_per_thread = 0;

        if (rank > 0) {
            MPI_Isend(local_grid.data(), N, MPI_CHAR, rank - 1, 0, MPI_COMM_WORLD, &request_send_1);
            MPI_Irecv(prev_last_row.data(), N, MPI_CHAR, rank - 1, 0, MPI_COMM_WORLD, &request_recv_1);
        }
        if (rank < size - 1) {
            MPI_Isend(local_grid.data() + (local_N - 1) * N, N, MPI_CHAR, rank + 1, 0, MPI_COMM_WORLD, &request_send_2);
            MPI_Irecv(next_first_row.data(), N, MPI_CHAR, rank + 1, 0, MPI_COMM_WORLD, &request_recv_2);
        }

        if (rank > 0) MPI_Wait(&request_recv_1, MPI_STATUS_IGNORE);
        if (rank < size - 1) MPI_Wait(&request_recv_2, MPI_STATUS_IGNORE);

        for (int i=1; i < local_N - 1; i++) {
            for (int j=0; j < N; j++) {

                tmp_sum = local_grid[(i - 1) * N + j] + local_grid[(i + 1) * N + j];
                if (j > 0) tmp_sum += local_grid[(i - 1) * N + j - 1] + local_grid[i * N + j - 1] + local_grid[(i + 1) * N + j - 1];
                if (j < N - 1) tmp_sum += local_grid[(i - 1) * N + j + 1] + local_grid[i * N + j + 1] + local_grid[(i + 1) * N + j + 1];
                
                if ((tmp_sum > 3) || (tmp_sum < 2) || ((tmp_sum == 2) && (local_grid[i * N + j] == 0))) new_local_grid[i * N + j] = 0;
                else {new_local_grid[i * N + j] = 1; n_alive_now_per_thread++;}
            }
        }

        for (int j=0; j < N; j++) {

            tmp_sum = local_grid[N + j];
            if (rank > 0) tmp_sum += prev_last_row[j];

            if (j > 0) {
                tmp_sum += local_grid[j - 1] + local_grid[N + j - 1];
                if (rank > 0) tmp_sum += prev_last_row[j - 1];
            }

            if (j < N - 1) {
                tmp_sum += local_grid[j + 1] + local_grid[N + j + 1];
                if (rank > 0) tmp_sum += prev_last_row[j + 1];
            }
            if ((tmp_sum > 3) || (tmp_sum < 2) || ((tmp_sum == 2) && (local_grid[j] == 0))) new_local_grid[j] = 0;
            else {new_local_grid[j] = 1; n_alive_now_per_thread++;}

            tmp_sum = local_grid[(local_N - 2) * N + j];
            if (rank > 0) tmp_sum += next_first_row[j];

            if (j > 0) {
                tmp_sum += local_grid[(local_N - 2) * N + j - 1] + local_grid[(local_N - 1) * N + j - 1];
                if (rank > 0) tmp_sum += next_first_row[j - 1];
            }

            if (j < N - 1) {
                tmp_sum += local_grid[(local_N - 2) * N + j + 1] + local_grid[(local_N - 1) * N + j + 1];
                if (rank > 0) tmp_sum += next_first_row[j + 1];
            }
            if ((tmp_sum > 3) || (tmp_sum < 2) || ((tmp_sum == 2) && (local_grid[(local_N - 1) * N + j] == 0))) new_local_grid[(local_N - 1) * N + j] = 0;
            else {new_local_grid[(local_N - 1) * N + j] = 1;  n_alive_now_per_thread++;}
        }

        MPI_Barrier(MPI_COMM_WORLD);
        local_grid.swap(new_local_grid);

        n_alive_past = n_alive_now;
        n_alive_now = 0;
        MPI_Reduce(&n_alive_now_per_thread, &n_alive_now, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Bcast(&n_alive_now, 1, MPI_INT, 0, MPI_COMM_WORLD);


    }

    double end_time = MPI_Wtime();
    double elapsed_time = end_time - start_time;
    double max_time;
    MPI_Reduce(&elapsed_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        std::cout << "K-ый шаг: " << n_iter << std::endl;
        std::cout << "Время выполнения: " << max_time << " секунд" << std::endl;
    }

    MPI_Finalize();
    return 0;
}
