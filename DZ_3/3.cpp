#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <papi.h>

struct CSRGraph {
    int n_vertex;
    std::vector<int> row_ptr;
    std::vector<int> col_ind;
    std::vector<float> weights;
};

int findMaxWeightVertex(const CSRGraph& graph) {
    int max_weight_vertex = -1;
    float max_weight = -1, sum_weight;

    for (int i=0; i < graph.n_vertex; i++) {
        sum_weight = 0;
        for (int j=graph.row_ptr[i]; j < graph.row_ptr[i+1]; j++) {
            if (graph.col_ind[j] % 2 == 0) sum_weight += graph.weights[j];
        }
        if (sum_weight > max_weight) {
            max_weight = sum_weight;
            max_weight_vertex = i;
        }
    }
    return max_weight_vertex;
}

int findMaxRankVertex(const CSRGraph& graph) {
    int max_rank_vertex = -1;
    float max_rank = -1, rank, W_vert;

    for (int i=0; i < graph.n_vertex; i++) {
        rank = 0;
        for (int j=graph.row_ptr[i]; j < graph.row_ptr[i+1]; j++) {
            W_vert = 0;
            for (int k=graph.row_ptr[graph.col_ind[j]]; k < graph.row_ptr[graph.col_ind[j]+1]; k++) {
                W_vert += graph.weights[k] * (graph.row_ptr[graph.col_ind[j]+1] - graph.row_ptr[graph.col_ind[j]]);
            }
            rank += graph.weights[j] * W_vert;
        }
        if (rank > max_rank) {
            max_rank = rank;
            max_rank_vertex = i;
        }
    }
    return max_rank_vertex;
}

int main() {
    CSRGraph graph;

    int n_vertex, n_friends_vertex, tmp;
    std::cout << "Введите число вершин графа" << std::endl;
    std::cin >> n_vertex;

    std::cout << "Введите " << n_vertex + 1 << " чисел для вектора row_ptr" << std::endl;
    for(int i=0; i < n_vertex + 1; i++) {
      std::cin >> tmp;
      graph.row_ptr.push_back(tmp);
    }

    n_friends_vertex = tmp;

    std::cout << "Введите " << n_friends_vertex << " чисел для вектора col_ind" << std::endl;
    for(int i=0; i < n_friends_vertex; i++) {
      std::cin >> tmp;
      graph.col_ind.push_back(tmp);
    }

    float float_tmp;

    std::cout << "Введите " << n_friends_vertex << " чисел для вектора weights" << std::endl;
    for(int i=0; i < n_friends_vertex; i++) {
      std::cin >> float_tmp;
      graph.weights.push_back(float_tmp);
    }

    graph.n_vertex = n_vertex;

    int Eventset = PAPI_NULL, code, retval = PAPI_NULL;
    long long values[3];
                
    PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT;
    PAPI_create_eventset(&Eventset);

    PAPI_add_event(Eventset, PAPI_L1_TCM);
    PAPI_add_event(Eventset, PAPI_L2_TCM);
    char event_name[] = "perf::PERF_COUNT_HW_CPU_CYCLES";
    PAPI_event_name_to_code(event_name, &code);
        
    PAPI_add_event(Eventset, code);
        
    PAPI_start(Eventset);
    std::cout << "Результат по первому алгоритму: " << findMaxWeightVertex(graph) + 1 << "-ая вершина графа" << std::endl;
    PAPI_stop(Eventset, values);
    std::cout << "L1 Data Misses: " << values[0] << std::endl;
    std::cout << "L2 Data Misses: " << values[1] << std::endl;
    std::cout << "PERF_COUNT_HW_CPU_CYCLES: " << values[2] << std::endl;

    PAPI_reset(Eventset);
    PAPI_start(Eventset); 
    std::cout << "Результат по второму алгоритму: " << findMaxRankVertex(graph) + 1 << "-ая вершина графа" << std::endl;
    PAPI_stop(Eventset, values);
    std::cout << "L1 Data Misses: " << values[0] << std::endl;
    std::cout << "L2 Data Misses: " << values[1] << std::endl;
    std::cout << "PERF_COUNT_HW_CPU_CYCLES: " << values[2] << std::endl;

    PAPI_destroy_eventset(&Eventset);
    PAPI_shutdown();
    return 0;
}
