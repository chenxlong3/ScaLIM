#ifndef GRAPH_BASE_H
#define GRAPH_BASE_H
#include "headers.h"
#include "CommonFuncs.h"
#include "IOcontroller.h"
typedef long long ll;
class GraphBase {
public:
    long m=0, n=0;
    Graph _G;
    Graph _inv_G;
    NodeList _indegree;
    VecDouble _weighted_in_degree;
    VecDouble _weighted_out_degree;
    string folder;
    string _prob_mode = "NONE";
    // vector<vector<double>> _Pr;
    // vector<vector<double>> _inv_Pr;
    GraphBase() {
        log_info("Empty Initialization GraphBase");
        return;
    }

    GraphBase(const string& folder_name, string f_name = "edgelist_ic.txt", bool format_g=true) {
        this->folder = folder_name;
        const string attr_filename = folder_name + "attr.txt";
        const string filename = folder_name + f_name;
        
        log_info("File:");
        cout << filename << endl;
        ifstream file(filename);
        ifstream attr_file(attr_filename);
        long u, v;
        double w;
        string s;
        while (attr_file >> s) {
            if (s.substr(0, 2) == "n=") {
                this->n = atoi(s.substr(2).c_str());
            } else if (s.substr(0, 2) == "m=") {
                this->m = atoi(s.substr(2).c_str());
            }
        }
        log_info("Number of Nodes", this->n);
        log_info("Number of Edges", this->m);
        attr_file.close();
        
        this->_G.resize(this->n);
        this->_inv_G.resize(this->n);
        this->_indegree.resize(this->n, 0);
    
        while(file >> u >> v >> w) {
            this->_G[u].push_back(Edge(v, w));
            this->_inv_G[v].push_back(Edge(u, w));
            this->_indegree[v]++;
            // Compute weighted degree
            // this->_weighted_in_degree[v] += w;
            // this->_weighted_out_degree[u] += w;
        }
        file.close();
        std::cout << "Finish reading the graph file." << std::endl;
        if (format_g) {
            log_info("--- Start formatting graph ---");
            TIO::save_graph_struct(folder_name, this->_G, false);
            // TIO::save_graph_struct(filename, this->_inv_G, true);
        }
    }

    GraphBase(const string& folder_name, string f_name, string prob_mode) {
        this->folder = folder_name;
        const string attr_filename = folder_name + "attr.txt";
        const string filename = folder_name + f_name;
        string post_fix = f_name.substr(f_name.size() - 3);
        // const string filename = "/root/research_proj/dyn_IM/my_proj/data/DBLP.txt";
        log_info("File:");
        cout << filename << endl;
        ifstream file(filename);
        ifstream attr_file(attr_filename);
        long u, v;
        double w;
        string s;
        while (attr_file >> s) {
            if (s.substr(0, 2) == "n=") {
                this->n = atoi(s.substr(2).c_str());
            } else if (s.substr(0, 2) == "m=") {
                this->m = atoi(s.substr(2).c_str());
            }
        }
        log_info("Number of Nodes", this->n);
        log_info("Number of Edges", this->m);
        attr_file.close();
        
        this->_G.resize(this->n);
        this->_inv_G.resize(this->n);
        
        this->_weighted_in_degree.resize(this->n);
        this->_weighted_out_degree.resize(this->n);
        std::fill(this->_indegree.begin(), this->_indegree.end(), 0);
        std::fill(this->_weighted_in_degree.begin(), this->_weighted_in_degree.end(), 0.0);
        std::fill(this->_weighted_out_degree.begin(), this->_weighted_out_degree.end(), 0.0);
        if (post_fix == "txt")
        {
            while(file >> u >> v >> w) {
                this->_G[u].push_back(Edge(v, w));
                // this->_inv_G[v].push_back(Edge(u, w));
                // this->_indegree[v]++;
                // Compute weighted degree
                // this->_weighted_in_degree[v] += w;
                // this->_weighted_out_degree[u] += w;
            }
            file.close();
        }
        else {
            TIO::load_graph_struct(folder_name, this->_G, false);
            // TIO::load_graph_struct(filename, this->_inv_G, true);
        }
        this->comp_indeg();
        this->set_probability(prob_mode);
        
        std::cout << "Finish reading the graph file." << std::endl;
    }

    void set_probability(string prob_mode = "NONE") {
        this->_prob_mode = prob_mode;
        if (prob_mode == "NONE") {
            log_info("--- Probability: Set by File ---");
        } else if (prob_mode == "UNI_RAND") {
            log_info("--- Probability: Uniformly Random ---");

            for (uint32_t i = 0; i < this->_G.size(); i++) {
                auto& nbrs = this->_G[i];
                for (auto& node : nbrs) {
                    double prob = dsfmt_gv_genrand_open_close();
                    node.second = prob;
                }
            }
        } else if (prob_mode == "UNI_001") {
            log_info("--- Probability Uniformly selected from {0.1, 0.01, 0.001}");
            VecDouble p_list = {0.1, 0.01, 0.001};
            for (uint32_t i = 0; i < this->_G.size(); i++) {
                auto& nbrs = this->_G[i];  // out-neighbors of node i
                for (auto& node : nbrs) {
                    double prob = p_list[dsfmt_gv_genrand_uint32_range(3)];
                    node.second = prob;
                }
            }
        } else if (prob_mode == "WC") {
            log_info("--- Probability: WC ---");
            for (uint32_t i = 0; i < this->_G.size(); i++) {
                auto& nbrs = this->_G[i];
                for (auto& node : nbrs) {
                    double prob = 1.0 / this->_indegree[node.first];  // For every nbr, the link is 1/nbr's indeg
                    node.second = prob;
                }
            }
            // for (uint32_t i = 0; i < this->_inv_G.size(); i++) {
            //     auto& nbrs = this->_inv_G[i];
            //     for (auto& node : nbrs) {
            //         node.second = double(1.0 / nbrs.size());  // for every inverse nbr, the edge probability is 1 / nbrs.size()
            //     }
            // }
        }
        else if (prob_mode == "001")
        {
            log_info("--- Probability: Uniformly Random ---");

            for (uint32_t i = 0; i < this->_G.size(); i++) {
                auto& nbrs = this->_G[i];
                for (auto& node : nbrs) {
                    double prob = 0.01;
                    node.second = prob;
                }
            }
        }
        else if (prob_mode == "0001")
        {
            log_info("--- Probability: Uniformly Random ---");

            for (uint32_t i = 0; i < this->_G.size(); i++) {
                auto& nbrs = this->_G[i];
                for (auto& node : nbrs) {
                    double prob = 0.001;
                    node.second = prob;
                }
            }
        }
        this->get_inv_G_and_stats();
        log_info("--- Finish setting probability distribution ---");
        return;
    }

    void comp_indeg() {
        this->_indegree.resize(this->n, 0);
        for (uint32_t i = 0; i < this->_G.size(); i++) {
            auto& nbrs = this->_G[i];  // out-neighbors of node i
            for (auto& node : nbrs) {
                this->_indegree[node.first]++;
            }
        }
        return;
    }

    void get_inv_G_and_stats() {
        // Refresh inv_G
        Graph().swap(this->_inv_G);
        this->_inv_G.resize(this->n);
        
        for (uint32_t i = 0; i < this->_G.size(); i++) {
            auto& nbrs = this->_G[i];  // out-neighbors of node i
            for (auto& node : nbrs) {
                double prob = node.second;
                this->_inv_G[node.first].push_back(Edge(i, prob));
                this->_weighted_in_degree[node.first] += prob;
                this->_weighted_out_degree[i] += prob;
            }
        }
        return;
    }

    void format_graph(const string& folder_name, string f_name = "edgelist_ic.txt") {
        this->folder = folder_name;
        const string attr_filename = folder_name + "attr.txt";
        const string filename = folder_name + f_name;
        log_info("File:");
        cout << filename << endl;
        ifstream file(filename);
        ifstream attr_file(attr_filename);
        long u, v;
        double w;
        string s;
        while (attr_file >> s) {
            if (s.substr(0, 2) == "n=") {
                this->n = atoi(s.substr(2).c_str());
            } else if (s.substr(0, 2) == "m=") {
                this->m = atoi(s.substr(2).c_str());
            }
        }
        log_info("Number of Nodes", this->n);
        log_info("Number of Edges", this->m);
        attr_file.close();
        
        this->_G.resize(this->n);
    
        while(file >> u >> v >> w) {
            this->_G[u].push_back(Edge(v, w));
        }
        file.close();
        std::cout << "Finish reading the graph file." << std::endl;
    
        log_info("--- Start formatting graph ---");
        TIO::save_graph_struct(folder_name, this->_G, false);
        return;
    }
};

static double inf_evaluate(const Graph& graph, const NodeList& vec_seed, const CascadeModel cascade_model, const uint32_t num_mc=10000, bool display=true) {
    /*
    Input: graph: Adjacency list
    vec_seed: vector of the seed nodes
    cascade_model: the information diffusion model we use
    num_mc: number of monte carlo simulations
    */
    if(display) std::cout << ">>>Evaluate influence...\n";
    std::queue<uint32_t> q;
    uint32_t num_V = graph.size();
    double inf = (double)vec_seed.size();       // k is already added here
    std::vector<uint32_t> vec_active;
    std::deque<bool> activated(num_V, false);
    
    // For LT specifically
    std::vector<uint32_t> last_time_vis(num_V, 0);     // Last time the node was visited
    std::vector<double> vec_node_threshold(num_V);
    std::vector<double> vec_activated_weight(num_V, 0.0);
    for(auto seed_id: vec_seed) activated[seed_id] = true;
    for(uint32_t i=0; i<num_mc; i++) {
        for(auto seed_id: vec_seed) {
            q.push(seed_id);
        }
        
        // IC
        if(cascade_model == IC) {
            while(!q.empty()) {
                uint32_t node_id = q.front();
                q.pop();
                // Traverse out neighbors
                for(auto& nbr: graph[node_id]) {
                    if(activated[nbr.first]) continue;
                    double rand_num = dsfmt_gv_genrand_open_close();
                    if(rand_num <= nbr.second) {
                        activated[nbr.first] = true;
                        vec_active.push_back(nbr.first);
                        q.push(nbr.first);
                    }
                }
            }
        }
        else if(cascade_model == LT) {
            while (!q.empty()) {
                uint32_t node_id = q.front();
                q.pop();
                for(auto& nbr: graph[node_id]) {
                    if(activated[nbr.first]) continue;
                    // Have entered a new iteration
                    if(last_time_vis[nbr.first] < i+1) {
                        last_time_vis[nbr.first] = i+1;
                        vec_node_threshold[nbr.first] = dsfmt_gv_genrand_open_close();
                        vec_activated_weight[nbr.first] = 0.0;
                    }
                    vec_activated_weight[nbr.first] += nbr.second;
                    if(vec_activated_weight[nbr.first] >= vec_node_threshold[nbr.first]) {
                        activated[nbr.first] = true;
                        vec_active.push_back(nbr.first);
                        q.push(nbr.first);
                    }
                }
            }
            
        }
        inf += (double)vec_active.size() / num_mc;
        // Reset
        // print_vector(vec_active);
        for(auto node_id: vec_active) activated[node_id] = false;
        vec_active.clear();
    }
    return inf;
}
#endif