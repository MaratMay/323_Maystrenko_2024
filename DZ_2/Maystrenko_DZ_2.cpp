#include <iostream>
#include "pthread.h"
#include <queue>
#include <cassert>

#define MESSAGES_PER_READER 2
#define MESSAGES_PER_WRITER 3
#define QUE_CAPACITY 10

template <typename T>
class MyConcurrentQueue {
    std::queue<T> queue;
    int queue_limit = QUE_CAPACITY;
    pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_cond_t cp = PTHREAD_COND_INITIALIZER;
    pthread_cond_t cg = PTHREAD_COND_INITIALIZER;
    public:  
        MyConcurrentQueue() {
            pthread_mutex_init(&mutex, nullptr);
            pthread_cond_init(&cp, NULL);
            pthread_cond_init(&cg, NULL);
        }

        ~MyConcurrentQueue() {
            pthread_mutex_destroy(&mutex);
            pthread_cond_destroy(&cp);
            pthread_cond_destroy(&cg);
        }

        void put(const T& item) {
            pthread_mutex_lock(&mutex);
            while (queue.size() == queue_limit) pthread_cond_wait(&cp, &mutex);
            queue.push(item);
            std::cout << "Write: " << item << std::endl;
            pthread_cond_signal(&cg);
            pthread_mutex_unlock(&mutex);
        }
        
        T get() {
            pthread_mutex_lock(&mutex);
            while (queue.empty()) pthread_cond_wait(&cg, &mutex);
            T item = queue.front();
            std::cout << "Read: " << item << std::endl;
            queue.pop();
            pthread_cond_signal(&cp);
            pthread_mutex_unlock(&mutex);
            return item;
        }
    };

void* writer(void* arg) {
    MyConcurrentQueue<int>* queue = static_cast<MyConcurrentQueue<int>*>(arg);
    for(int j=0; j < MESSAGES_PER_WRITER; j++) queue -> put(random() % 1000);
    return nullptr;
}

void* reader(void* arg) {
    MyConcurrentQueue<int>* queue = static_cast<MyConcurrentQueue<int>*>(arg);
    for(int j=0; j < MESSAGES_PER_READER; j++) queue -> get();
    return nullptr;
}

int main() {
    MyConcurrentQueue<int> my_queue;
    int num_readers, num_writers;
    std::cout << "Введите количество писателей и читателей (например, 3 3): ";
    std::cin >> num_writers >> num_readers;
    pthread_t readers[num_readers];
    pthread_t writers[num_writers];

    for (int i=0; i < num_writers; i++) pthread_create(&writers[i], nullptr, writer, &my_queue);

    for (int i=0; i < num_readers; i++) pthread_create(&readers[i], nullptr, reader, &my_queue);

    for (int i=0; i < num_writers; i++) pthread_join(writers[i], nullptr);

    for (int i=0; i < num_readers; i++) pthread_join(readers[i], nullptr);

    return 0;
}
