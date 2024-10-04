#include <iostream>
#include <pthread.h>
#include <queue>
#include <unistd.h>

#define NO_MORE_WRITERS 4
#define NO_MORE_READERS 3
#define WRITING 2
#define READING 1

#define MESSAGES_PER_READER 2
#define MESSAGES_PER_WRITER 3
#define QUE_CAPACITY 10

template <typename T>
class MyConcurrentQueue {
    std::queue<T> queue;
    int num_writers, num_readers;
    int doing = WRITING;
    size_t maxSize;
    pthread_mutex_t mutex;
    public:
        MyConcurrentQueue(size_t _maxSize, int a, int b) : maxSize(_maxSize), num_writers(a * MESSAGES_PER_WRITER), num_readers(b * MESSAGES_PER_READER) {
            pthread_mutex_init(&mutex, nullptr);
        }

        ~MyConcurrentQueue() {
            pthread_mutex_destroy(&mutex);
        }

        void put(const T& item) {
            pthread_mutex_lock(&mutex);
            if (num_readers == 0) doing = NO_MORE_READERS;
            while(doing == READING) {
                pthread_mutex_unlock(&mutex);
                pthread_mutex_lock(&mutex);
            }
            if (doing == NO_MORE_READERS) {
                if (queue.size() >= maxSize) {
                    std::cout << "We're out of space, sorry >_<" << std::endl;
                } else {
                    queue.push(item);
                    std::cout << "Write: " << item << std::endl;
                } pthread_mutex_unlock(&mutex);
            } else {
                queue.push(item);
                std::cout << "Write: " << item << std::endl;
                num_writers -= 1;
                doing = (num_writers) ? READING : NO_MORE_WRITERS;
                pthread_mutex_unlock(&mutex);
            }
        }

        T get() {
            pthread_mutex_lock(&mutex);
            if (num_writers == 0) doing = NO_MORE_WRITERS;
            while(doing == WRITING) {
                pthread_mutex_unlock(&mutex);
                pthread_mutex_lock(&mutex);
            }
            T item = -1;
            if (doing == NO_MORE_WRITERS) {
                if (queue.size() == 0) {
                    std::cout << "We're out of elements, sorry >_<" << std::endl;
                } else {
                    T item = queue.front();
                    queue.pop();
                    std::cout << "Read: " << item << std::endl;
                } pthread_mutex_unlock(&mutex);
            } else {
                T item = queue.front();
                queue.pop();
                std::cout << "Read: " << item << std::endl;
                num_readers -= 1;
                doing = (num_readers) ? WRITING : NO_MORE_READERS;
                pthread_mutex_unlock(&mutex);
            } return item;
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

    int num_readers, num_writers;
    std::cout << "Введите количество писателей и читателей (например, 3 3): ";
    std::cin >> num_writers >> num_readers;

    MyConcurrentQueue<int> queue(QUE_CAPACITY, num_writers, num_readers);

    pthread_t readers[num_readers];
    pthread_t writers[num_writers];

    for (int i=0; i < num_writers; i++) pthread_create(&writers[i], nullptr, writer, &queue);

    for (int i=0; i < num_readers; i++) pthread_create(&readers[i], nullptr, reader, &queue);

    for (int i=0; i < num_writers; i++) pthread_join(writers[i], nullptr);

    for (int i=0; i < num_readers; i++) pthread_join(readers[i], nullptr);

    return 0;
}
