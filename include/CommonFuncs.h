#ifndef COMMON_FUNCS_H
#define COMMON_FUNCS_H
#include "headers.h"
template <typename _Ty>
static inline void log_info(_Ty val) {
    std::cout << val << std::endl;
}

/// Log information
template <typename _Ty>
static inline void log_info(const std::string title, _Ty val) {
    std::cout << title << ": " << val << std::endl;
}

/// Math, pow2
static inline double pow2(const double t) {
    return t * t;
}

/// Math, log2
static inline double log2(const size_t n) {
    return log(n) / log(2);
}

/// Math, logcnk
static inline double logcnk(const size_t n, size_t k) {
    k = k < n - k ? k : n - k;
    double res = 0;
    for (auto i = 1; i <= k; i++) res += log(double(n - k + i) / i);
    return res;
}

void print_vector(std::vector<uint32_t> v) {
    for(auto node_id: v) {
        std::cout << node_id << ' ';
    }
    std::cout << '\n';
    return;
}

void print_vector(std::vector<size_t> v) {
    for(auto node_id: v) {
        std::cout << node_id << ' ';
    }
    std::cout << '\n';
    return;
}

void print_vector(VecDouble v) {
    for(auto x: v) {
        std::cout << x << ' ';
    }
    std::cout << '\n';
    return;
}

void print_vector(std::vector<UVWEdge> v) {
    for(auto& tuple: v) {
        std::cout << tuple.first << ' ' << tuple.second.first << ' ' << tuple.second.second << std::endl;
    }
    std::cout << '\n';
    return;
}

void print_vector(std::vector<uint32_t> v, uint32_t max_len) {
    if (max_len > v.size()) max_len = v.size();
    for(auto i=0; i<max_len; i++) {
        std::cout << v[i] << ' ';
    }
    std::cout << '\n';
    return;
}


static size_t gen_rand_node_by_weight_LT(const EdgeList& edges) {
    const double weight = dsfmt_gv_genrand_open_close();
    size_t min_idx = 0, max_idx = edges.size() - 1;
    if (weight < edges.front().second) return 0;
    if (weight > edges[max_idx].second) return edges.size();
    while (max_idx > min_idx) {
        const size_t mid_idx = (min_idx + max_idx) / 2;
        const auto mid_weight = edges[mid_idx].second;
        if (weight <= mid_weight) max_idx = mid_idx;
        else min_idx = mid_idx+1;
    }
    return max_idx;
}

template <typename T>
double comp_mean(std::vector<T> v) {
    return (double)std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}

double variance(VecLargeNum a, int n) {
    // Compute mean (average of elements)
    int sum = 0;
    for (int i = 0; i < n; i++)
        sum += a[i];
    double mean = (double)sum /
                  (double)n;

    // Compute sum squared
    // differences with mean.
    double sqDiff = 0;
    for (int i = 0; i < n; i++)
        sqDiff += (a[i] - mean) *
                  (a[i] - mean);
    return sqDiff / n;
}

std::vector<UVWEdge> read_UVWEdges(string filename) {
    uint32_t u, v;
    double w;
    ifstream myfile(filename);
    std::vector<UVWEdge> res;
    while (myfile >> u >> v >> w) {
        res.push_back(make_pair(u, Edge(v, w)));
    }
    myfile.close();
    return res;
}

void write_UVWEdges(std::vector<UVWEdge>& tuple_list, string filename) {
    ofstream outfile(filename);
    for (auto& tuple : tuple_list) {
        outfile << tuple.first << '\t' << tuple.second.first << '\t' << tuple.second.second << '\n';
    }
    outfile.close();
    log_info("Finish saving the edges to " + filename);
    return;
}

void ExitMessage(string msg) {
    cout << msg << endl;
    exit(1);
}

bool make_dir(string dir_path) {
    bool ret = false;
    DIR* cur_dir = nullptr;
    cur_dir = opendir(dir_path.c_str());
    if (cur_dir == NULL)
    {
        log_info("--- This path does not exist, creating ---");
        ret = (mkdir(dir_path.c_str(), 0777) == 0);
        if (!ret)   log_info("--- mkdir FAILED ---");
        else log_info("--- mkdir success ---");
    }
    else {
        ret = true;
        log_info("--- This dir exists ---");
        closedir(cur_dir);
    }
    return ret;
}

void save_vector(NodeList v, string out_file) {
    ofstream myfile(out_file);
    for (auto& i : v) {
        myfile << i << '\n';
    }
    return;
    
}

bool check_file_exist(string file_path) {
    ifstream ifile;
    ifile.open(file_path);
    if(ifile) {
        ifile.close();
        return true;
    } else {
        return false;
    }
}

inline void make_min_heap(VecLargeNum& vec)
{
	// Min heap
	const auto size = vec.size();
	if (2 <= size)
	{
		for (auto hole = (size + 1) / 2; hole--;)
		{
			const auto val = vec[hole];
			size_t i, child;
			for (i = hole; i * 2 + 1 < size; i = child)
			{
				// Find smaller child
				child = i * 2 + 2;
				if (child == size || vec[child - 1] < vec[child])
				{
					// One child only or the left child is smaller than the right one
					--child;
				}

				// Percolate one level
				if (vec[child] < val)
				{
					vec[i] = vec[child];
				}
				else
				{
					break;
				}
			}
			vec[i] = val;
		}
	}
}

/// Replace the value for the first element and down-heap this element.
inline void min_heap_replace_min_value(VecLargeNum& vec, const size_t& val)
{
	// Increase the value of the first element
	const auto size = vec.size();
	size_t i, child;
	for (i = 0; i * 2 + 1 < size; i = child)
	{
		// Find smaller child
		child = i * 2 + 2;
		if (child == size || vec[child - 1] < vec[child])
		{
			// One child only or the left child is smaller than the right one
			--child;
		}

		// Percolate one level
		if (vec[child] < val)
		{
			vec[i] = vec[child];
		}
		else
		{
			break;
		}
	}
	vec[i] = val;
}

VecLargeNum max_cover_by_heap(VecVecLargeNum& vec_u, VecVecLargeNum& vec_v, uint32_t target_size) {
    size_t num_u = vec_u.size();
    size_t num_v = vec_v.size();
    VecLargeNum res;
    // Use a heap to store the <node, coverage> pair
    std::priority_queue<PairIntSizet, std::vector<PairIntSizet>, CompareBySecond> heap;
    VecLargeNum coverage(num_u, 0);       // coverage[v] is the RR sets covered by v
    int cnt = 0;
    for(uint32_t i=0; i<num_u; i++) {
        // store coverage
        size_t node_i_cov = vec_u[i].size();
        PairIntSizet tmp(make_pair(i, node_i_cov));
        heap.push(tmp);
        coverage[i] = node_i_cov;
        
    }

    VecBool RRsets_mark(vec_v.size(), false);
    VecBool node_mark(num_u, false);   

    uint32_t max_idx;
    size_t cov_num = 0;
    while (res.size() < target_size) {
        PairIntSizet top = heap.top();
        heap.pop();
        // Lazy Update
        if (top.second > coverage[top.first]) {
            // Update coverage of top
            top.second = coverage[top.first];
            heap.push(top);
            continue;
        }
        max_idx = top.first;
        VecLargeNum& e = vec_u[max_idx];     // e: the RR sets covered by the node max_idx
        cov_num += coverage[max_idx];
        res.push_back(max_idx);
        node_mark[max_idx] = true;

        
        // After selecting one node, we need to remove the covered RR sets from the coverage graph
        // e: the RR sets covered by the node max_idx
        for (uint32_t j=0; j<e.size(); j++) {
            if (RRsets_mark[e[j]]) continue;        // If the RR set has been removed
            VecLargeNum node_list = vec_v[e[j]];
            
            
            for (uint32_t l=0; l<node_list.size(); ++l){
                if (!node_mark[node_list[l]]) {
                    coverage[node_list[l]]--;
                }
            }
            RRsets_mark[e[j]] = true;
        }
    }
    return res;
}

double comp_cov(VecVecLargeNum& vec_u, VecVecLargeNum& vec_v, VecLargeNum vec_seed)
{
    ASSERT(vec_seed.size() > 0);
    std::vector<bool> vecBoolVst = std::vector<bool>(vec_v.size());
    std::vector<bool> vecBoolSeed(vec_u.size());
    for (auto seed : vec_seed) vecBoolSeed[seed] = true;
    for (auto seed : vec_seed)
    {
        for (auto node : vec_u[seed])
        {
            vecBoolVst[node] = true;
        }
    }
    return 1.0 * std::count(vecBoolVst.begin(), vecBoolVst.end(), true);
}

// size_t hash_RRset(VecuInt32& RRset) {
//     size_t seed=0;
//     if (RRset.size() == 1) {
//         boost::hash_combine(seed, RRset[0]);
//     }
//     if (RRset.size() == 2) {
//         boost::hash_combine(seed, std::min(RRset[0], RRset[1]));
//         boost::hash_combine(seed, std::max(RRset[0], RRset[1]));
//     }
//     // boost::hash<std::string> string_hash;
//     return seed;
// }
#endif