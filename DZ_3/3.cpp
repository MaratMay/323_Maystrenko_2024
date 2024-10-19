#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <papi.h>

class CSR_graph {
    int row_count; //number of vertices in graph
    unsigned int col_count; //number of edges in graph
    
    std::vector<unsigned int> row_ptr;
    std::vector<int> col_ids;
    std::vector<double> vals;

public:

    void read(const char* filename) {
        FILE *graph_file = fopen(filename, "rb");
        fread(reinterpret_cast<char*>(&row_count), sizeof(int), 1, graph_file);
        fread(reinterpret_cast<char*>(&col_count), sizeof(unsigned int), 1, graph_file);

        std::cout << "Row_count = " << row_count << ", col_count = " << col_count << std::endl;
        
        row_ptr.resize(row_count + 1);
        col_ids.resize(col_count);
        vals.resize(col_count);
         
        fread(reinterpret_cast<char*>(row_ptr.data()), sizeof(unsigned int), row_count + 1, graph_file);
        fread(reinterpret_cast<char*>(col_ids.data()), sizeof(int), col_count, graph_file);
        fread(reinterpret_cast<char*>(vals.data()), sizeof(double), col_count, graph_file);
        fclose(graph_file);
    }

    void print_vertex(int idx) {
        for (int col = row_ptr[idx]; col < row_ptr[idx + 1]; col++) {
            std::cout << col_ids[col] << " " << vals[col] <<std::endl;
        }
        std::cout << std::endl;
    }

    void reset() {
        row_count = 0;
        col_count = 0;
        row_ptr.clear();
        col_ids.clear();
        vals.clear();
    }
    int findMaxWeightVertex() {
        int max_weight_vertex = -1;
        float max_weight = -1, sum_weight;

        for (int i=0; i < row_count; i++) {
            sum_weight = 0;
            for (int j=row_ptr[i]; j < row_ptr[i+1]; j++) {
                if (col_ids[j] % 2 == 0) sum_weight += vals[j];
            }
            if (sum_weight > max_weight) {
                max_weight = sum_weight;
                max_weight_vertex = i;
            }
        }
        return max_weight_vertex;
    }

    int findMaxRankVertex() {
        int max_rank_vertex = -1;
        float max_rank = -1, rank, W_vert;

        for (int i=0; i < row_count; i++) {
            rank = 0;
            for (int j=row_ptr[i]; j < row_ptr[i+1]; j++) {
                W_vert = 0;
                for (int k=row_ptr[col_ids[j]]; k < row_ptr[col_ids[j]+1]; k++) {
                    W_vert += vals[k] * (row_ptr[col_ids[j]+1] - row_ptr[col_ids[j]]);
                }
                rank += vals[j] * W_vert;
            }
            if (rank > max_rank) {
                max_rank = rank;
                max_rank_vertex = i;
            }
        }
        return max_rank_vertex;
    }
}; 

#define N_TESTS 5

int main () {
    const char* filenames[N_TESTS];
    filenames[0] = "synt";
    filenames[1] = "road_graph";
    filenames[2] = "stanford";
    filenames[3] = "youtube";
    filenames[4] = "syn_rmat";

    /* https://drive.google.com/file/d/183OMIj56zhqN12Aui1kxv76_Ia0vTPIF/view?usp=sharing архив с тестами, 
        распаковать командой tar -xzf 
    */

    int Eventset = PAPI_NULL, code, retval = PAPI_NULL;
    long long values[3];
                
    PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT;
    PAPI_create_eventset(&Eventset);

    PAPI_add_event(Eventset, PAPI_L1_TCM);
    PAPI_add_event(Eventset, PAPI_L2_TCM);
    char event_name[] = "perf::PERF_COUNT_HW_CACHE_REFERENCES";
    PAPI_event_name_to_code(event_name, &code);
        
    PAPI_add_event(Eventset, code);

    for (int n_test = 0; n_test < N_TESTS; n_test++) {
        CSR_graph a;
        a.read(filenames[n_test]);
        PAPI_start(Eventset);
        std::cout << "Результат по первому алгоритму: " << a.findMaxWeightVertex() + 1 << "-ая вершина графа" << std::endl;
        PAPI_stop(Eventset, values);
        std::cout << "L1 Data Misses: " << values[0] << std::endl;
        std::cout << "L2 Data Misses: " << values[1] << std::endl;
        std::cout << "PERF_COUNT_HW_CACHE_REFERENCES: " << values[2] << std::endl;

        PAPI_reset(Eventset);
        PAPI_start(Eventset); 
        std::cout << "Результат по второму алгоритму: " << a.findMaxRankVertex() + 1 << "-ая вершина графа" << std::endl;
        PAPI_stop(Eventset, values);
        std::cout << "L1 Data Misses: " << values[0] << std::endl;
        std::cout << "L2 Data Misses: " << values[1] << std::endl;
        std::cout << "PERF_COUNT_HW_CACHE_REFERENCES: " << values[2] << std::endl;
    }
    PAPI_destroy_eventset(&Eventset);
    PAPI_shutdown();
    return 0;
}
