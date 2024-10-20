
#include <iostream>
#include <random>
#include <sys/time.h>
#include <immintrin.h> 
#include <string.h>
#define N 640

double get_time() 
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return t.tv_sec + t.tv_usec / 1000.0;
}

void fillMatrix(float *matrix) {
    std::default_random_engine generator;
    std::uniform_real_distribution<float> distribution(-1.0f, 1.0f);
    for (int i = 0; i < N * N; i++) matrix[i] = distribution(generator);
}

void Multiply_irl(const float *A, const float *B, float *C) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            float sum = 0.0f;
            for (int k = 0; k < N; k++) {
                sum += A[i * N + k] * B[k * N + j];
            }
            C[i * N + j] = sum;
        }
    }
}

void Multiply_vec(const float *A, const float *B, float *C) {
    __m128 aVec, bVec, cVec, tmpVec;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) { //зафиксировали элемент матрицы
            aVec = _mm_broadcast_ss(A + i * N + j);
            for (int k = 0; k < N; k += 4) { //начинаем его умножать, но уже ВЕКТОРНО!
                bVec = _mm_loadu_ps(&B[j * N + k]);
                tmpVec = _mm_mul_ps(aVec, bVec);

                cVec = _mm_loadu_ps(&C[i * N + k]);
                cVec = _mm_add_ps(cVec, tmpVec);
                _mm_storeu_ps(&C[i * N + k], cVec);
            }
        }
    }
}

int main() {
    float A[N * N], B[N * N], C[N * N], D[N * N];
    memset(C, 0, N * N);
    memset(D, 0, N * N);
    fillMatrix(A);
    fillMatrix(B);

    double time_irl = get_time();
    Multiply_irl(A, B, C);
    time_irl = get_time() - time_irl;

    double time_vec = get_time();
    Multiply_vec(A, B, D);
    time_vec = get_time() - time_vec;

    std::cout << "Basic Multiply time: " << time_irl << " ms" << std::endl;
    std::cout << "Vectorized Multiply time: " << time_vec << " ms" << std::endl;

    return 0;
}
