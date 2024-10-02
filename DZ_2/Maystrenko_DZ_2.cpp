#include <iostream>
#include <pthread.h>
#include <queue>
#include <unistd.h>

#define NO_MORE_WRITERS 4
#define NO_MORE_READERS 3
#define WRITING 2
#define READING 1

template <typename T>
class MyConcurrentQueue {
    std::queue<T> queue;
    int num_writers, num_readers;
    char doing = WRITING;
    size_t maxSize;
    pthread_mutex_t mutex_write;
    pthread_mutex_t mutex_read;
    pthread_cond_t go_write;
    pthread_cond_t go_read;
    public:
        MyConcurrentQueue(size_t _maxSize, int a, int b) : maxSize(_maxSize), num_writers(a), num_readers(b) {
            pthread_mutex_init(&mutex_write, nullptr);
            pthread_mutex_init(&mutex_read, nullptr);
            pthread_cond_init(&go_write, nullptr);
            pthread_cond_init(&go_read, nullptr);
        }

        ~MyConcurrentQueue() {
            pthread_mutex_destroy(&mutex_write);
            pthread_mutex_destroy(&mutex_read);
            pthread_cond_destroy(&go_read);
            pthread_cond_destroy(&go_write);
        }

        void put(const T& item) {
            pthread_mutex_lock(&mutex_write);
            if (num_readers == 0) doing = NO_MORE_READERS;
            if (doing == NO_MORE_READERS) {
                if (queue.size() >= maxSize) {
                    std::cout << "We're out of space, sorry >_<" << std::endl;
                } else {
                    queue.push(item);
                    std::cout << "Write: " << item << std::endl;
                } pthread_mutex_unlock(&mutex_write);
            } else {
                while (doing == READING) pthread_cond_wait(&go_write, &mutex_write);
                queue.push(item);
                std::cout << "Write: " << item << std::endl;
                num_writers -= 1;
                doing = (num_writers) ? READING : NO_MORE_WRITERS;
                pthread_cond_signal(&go_read);
                pthread_mutex_unlock(&mutex_read);
            }
        }

        T get() {
            pthread_mutex_lock(&mutex_read);
            T item = -1;
            if (num_writers == 0) doing = NO_MORE_WRITERS;
            if (doing == NO_MORE_WRITERS) {
                if (queue.size() == 0) {
                    std::cout << "We're out of elements, sorry >_<" << std::endl;
                } else {
                    T item = queue.front();
                    queue.pop();
                    std::cout << "Read: " << item << std::endl;
                } 
                pthread_mutex_unlock(&mutex_read);
            } else {
                while (doing == WRITING) pthread_cond_wait(&go_read, &mutex_read);
                T item = queue.front();
                queue.pop();
                std::cout << "Read: " << item << std::endl;
                num_readers -= 1;
                doing = (num_readers) ? WRITING : NO_MORE_READERS;
                pthread_cond_signal(&go_write);
                pthread_mutex_unlock(&mutex_write);
            } return item;
        }
    };

void* writer(void* arg) {
    MyConcurrentQueue<int>* queue = static_cast<MyConcurrentQueue<int>*>(arg);
    queue-> put(random() % 1000);
    return nullptr;
}

void* reader(void* arg) {
    MyConcurrentQueue<int>* queue = static_cast<MyConcurrentQueue<int>*>(arg);
    queue-> get();
    return nullptr;
}

int main() {

    int num_readers, num_writers;
    std::cout << "Введите количество писателей и читателей (например, 3 3): ";
    std::cin >> num_readers >> num_writers;

    MyConcurrentQueue<int> queue(5, num_writers, num_readers);

    pthread_t readers[num_readers];
    pthread_t writers[num_writers];

    for (int i=0; i < num_writers; i++) pthread_create(&writers[i], nullptr, writer, &queue);
    for (int i=0; i < num_readers; i++) pthread_create(&readers[i], nullptr, reader, &queue);
    for (int i=0; i < num_writers; i++) pthread_join(writers[i], nullptr);
    for (int i=0; i < num_readers; i++) pthread_join(readers[i], nullptr);
    return 0;
}
