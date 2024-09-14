#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <sys/time.h>

struct arguments 
{
    double h;
    double start;
    double n_per_thread;
    double ans;
};

void* func(void* arg) {
    double x, n_per_thread, sum, h;
    struct arguments* my_args = (struct arguments*)(arg);

    h = my_args->h;
    x = my_args->start;
    n_per_thread = my_args->n_per_thread;
    sum = 4.0 / (1 + x * x);
    
    for (int i=1; i < n_per_thread; i++) {
        x += h;
        sum += 4.0 / (1 + x * x);
    }
    my_args->ans = sum;
    pthread_exit(NULL);
}

double get_time() 
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return t.tv_sec + t.tv_usec / 1000000.0;
}

int main(int args, char **argv)
{
    double a=0, b=1, ans=0, h, n_per_thread;

    int n = atof(argv[1]);
    int n_threads = atoi(argv[2]);

    h = ((double)(b - a)) / n;
    n_per_thread = (double) n / n_threads;

    pthread_t my_threads[n_threads];
    struct arguments my_params[n_threads];
    
    double start = get_time();

    for(int i=0; i < n_threads; i++) {
        my_params[i].h = h;
        my_params[i].start = a + i * ((double)(b - a) / n_threads);
        my_params[i].n_per_thread = n_per_thread;
        my_params[i].ans = 0;
        pthread_create(&my_threads[i], NULL, &func, &my_params[i]);
    }
    void* ans_from_thread;

    for (int i = 0; i < n_threads; i++) 
    {
        pthread_join(my_threads[i], NULL);
        ans += my_params[i].ans;
    }

    ans *= h;
    printf("%f\n", ans);

    double end = get_time();
    double time_used = end - start;

    printf("Time used = %f\n", time_used);

    return 0;
}
