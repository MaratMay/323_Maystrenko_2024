#include <iostream>
#include <random>
#include <chrono>
#include <immintrin.h> 
#include <string.h>
#define N 2048

using namespace std::chrono;

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
    float *A = (float *) malloc(sizeof(float) * N * N);
    float *B = (float *) malloc(sizeof(float) * N * N);
    float *C = (float *) malloc(sizeof(float) * N * N);
    float *D = (float *) malloc(sizeof(float) * N * N);
    memset(C, 0, N * N);
    memset(D, 0, N * N);
    fillMatrix(A);
    fillMatrix(B);

    high_resolution_clock::time_point t1 = high_resolution_clock::now();    
    Multiply_irl(A, B, C);
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double, std::milli> time_span = duration_cast<duration<double, std::milli>>(t2 - t1);
    std::cout << "Basic Multiply time: " << time_span.count() << " ms" << std::endl;

    t1 = high_resolution_clock::now();    
    Multiply_vec(A, B, D);
    t2 = high_resolution_clock::now();
    time_span = duration_cast<duration<double, std::milli>>(t2 - t1);
    std::cout << "Vectorized Multiply time: " << time_span.count() << " ms" << std::endl;

    //свободу попугаям !
    free(A);
    free(B);
    free(C);
    free(D);

    return 0;
}
