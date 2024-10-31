#include <iostream>
#include <vector>
#include <algorithm>
#include <omp.h>
#include <cstdlib>
#include <ctime>
#include <iomanip>

void merge(std::vector<int>& arr, int left, int mid, int right) 
{
    //std::cout << "left= " << left << " right= " << right << std::endl;
    std::vector<int> tmp(right - left + 1);
    int i=left, j=mid+1, k=0;

    while (i <= mid && j <= right) {
        if (arr[i] <= arr[j]) tmp[k++] = arr[i++];
        else tmp[k++] = arr[j++];
    }

    while (i <= mid) tmp[k++] = arr[i++];
    while (j <= right) tmp[k++] = arr[j++];

    for (i = 0; i < k; i++) arr[left + i] = tmp[i];
}

void parallelMergeSort(std::vector<int>& arr, int left, int right) 
{
    if (left < right) {
        int mid = left + (right - left) / 2;

        #pragma omp task shared(arr) if(right - left > 1000) //если меньше, то не создаём task
        parallelMergeSort(arr, left, mid);

        #pragma omp task shared(arr) if(right - left > 1000)
        parallelMergeSort(arr, mid + 1, right);

        #pragma omp taskwait //ждёт детишек
        merge(arr, left, mid, right);
    }
}

int compare(const void* a, const void* b) 
{
    return (*(int*)a - *(int*)b);
}

int main() {
    int N, p;
    std::cout << "Введите размер массива (N) и количество потоков (p): ";
    std::cin >> N >> p;

    std::vector<int> arr(N);
    std::vector<int> arr_copy(N);
    srand(time(0));

    for (int i = 0; i < N; i++) {
        arr[i] = rand() % 1000;
        arr_copy[i] = arr[i];
    }

    double start_time = omp_get_wtime();
    #pragma omp parallel num_threads(p)
    {
        #pragma omp single //запуск рекурсии
        parallelMergeSort(arr, 0, N - 1);
    }
    double end_time = omp_get_wtime();
    double parallel_time = end_time - start_time;

    start_time = omp_get_wtime();
    qsort(arr_copy.data(), N, sizeof(int), compare);
    end_time = omp_get_wtime();
    double qsort_time = end_time - start_time;

    bool parallel_sorted = std::is_sorted(arr.begin(), arr.end());
    std::cout << "Параллельная сортировка " << (parallel_sorted ? "прошла успешно" : "не удалась") << "! " << (parallel_sorted ? "😍👍" : "😱❌") << std::endl;

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Время параллельной сортировки: " << parallel_time << " секунд 🚀" << std::endl;
    std::cout << "Время qsort(): " << qsort_time << " секунд 🐢" << std::endl;
    std::cout << "Параллельная сортировка быстрее в " << qsort_time / parallel_time << " раз! 🚀💨" << std::endl;
    
    return 0;
}
