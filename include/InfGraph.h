#ifndef INF_GRAPH_H
#define INF_GRAPH_H
#include "GraphBase.h"
#include "Timer.h"
#include "ResultInfo.h"

// All arguments passed to func
class Argument{
public:
    uint32_t rand_seed = 0;
    string folder_name = "";
    string graph_file = "";
    string probability_mode = "";
    
    CascadeModel model = IC;
    string str_model = "";
    
    string seed_mode = "";
    string method="";

    uint32_t k_seed=0, k_edges = 0, beta = 1;
    double epsilon = -0.6;
    double delta=0;
    
    bool fast_truncated = true;
    bool pruned = false;
    uint32_t num_cand_edges = 0;
    size_t num_samples = 0;     // num_mc or num of RR sets
    
    Argument() {

    };

    Argument(uint32_t _k_seed, uint32_t k_links, uint32_t _num_cand_edges, double eps, double _delta, CascadeModel _model, uint32_t _beta) {
        this->k_seed = _k_seed;
        this->k_edges = k_links;
        this->num_cand_edges = _num_cand_edges;
        this->epsilon = eps;
        this->delta = _delta;
        this->model = _model;
        this->beta = _beta;
    };
    void check_arguments_eligible() {
        if (this->folder_name == "")
            ExitMessage("argument dataset missing");
        if (this->k_seed == 0)
            ExitMessage("argument k_seed missing");
        if (this->k_edges == 0)
            ExitMessage("argument k_edges missing");
        if (this->epsilon < -0.5)
            ExitMessage("argument epsilon missing");
        if (this->delta == 0)
            ExitMessage("argument delta missing");

        if (this->rand_seed == 0)
            log_info("argument rand_seed missing, 0 as default");
        if (this->str_model == "")
            log_info("argument model missing, IC as default");
        if (this->method == "") {
            log_info("argument method missing, ScaLIM as default");
            this->method = "ScaLIM";
        }
        if (this->seed_mode == "") {
            log_info("argument seed_mode missing, IM as default");
            this->seed_mode = "IM";
        }
        if (this->probability_mode == "") {
            log_info("argument probability_mode missing, NONE as default");
            this->probability_mode == "NONE";
        }
        if (this->graph_file == "") {
            log_info("argument graph_file missing, vec.graph as default");
            this->graph_file = "vec.graph";
        }
        if (this->num_samples > 0) {
            this->epsilon = 0;
        }
        return;
    }
};

class InfGraph: public GraphBase {
    public:
        VecVecInt32 _RRsets;      // _RRsets[i] is the vector storing nodes of the i-th RRset
        VecVecLargeNum _node_to_RRsets;  // _node_to_RRsets[i] is the vector storing the indices of RR sets containing node i
        VecVecInt32 _RRsets_val;      // _RRsets[i] is the vector storing nodes of the i-th RRset
        VecVecLargeNum _node_to_RRsets_val;  // _node_to_RRsets[i] is the vector storing the indices of RR sets containing node i
        VecVecInt32 _RRsets_delta;       // Marginal RR sets, all RR sets not covered by S
        VecVecLargeNum _node_to_RRsets_delta;  // node to RR sets not covered by S

        CascadeModel _casc_model = IC;
        NodeList _seed_set_IM;
        NodeList _seed_set_to_augment;
        // VecLargeNum _node_sampled_times; // For Debug
        std::vector<uint32_t> _RRset_traversed_edges;
        VecBool _vec_bool_vis;
        VecBool _vec_bool_seed;
        
        size_t _cur_RRsets_num=0;
        size_t _cur_RRsets_num_val=0;
        LL _num_traversed_edges=0, _tmp_num_traversed_edges=0;
        VecSetInt _SET_RRsets;
        VecLargeNum _marginal_cov;
        VecBool _vec_bool_RRset_cov_by_S;           // Bool, indicating whether a RR set is covered by S
        VecBool _vec_bool_RRset_val_cov_by_S;           // Bool, indicating whether a RR set is covered by S

        VecBool _vec_bool_cand_edges_selected;
        
        string _seed_mode = "RAND";

        std::vector<UVWEdge> _vec_cand_edges;       // vector of all candidate edges
        std::vector<UVWEdge> _vec_selected_edges;   // vector of selected edges
        std::vector<UVWEdge> _vec_selected_edges_MCGreedy;   // vector of selected edges by degree
        std::vector<std::priority_queue<PairSizetDouble, std::vector<PairSizetDouble>, CompareBySecondDoubleForSizet>> _vec_node_prob_heap;
        std::vector<std::vector<PairSizetDouble>> _vec_node_prob_vec;

        // vector mapping from candidate edge index to target node
        VecVecLargeNum _vec_node_to_cand_edge_idx;
        
        string _cur_working_folder = "";
        string _cand_edges_filename = "";
        string _seed_filename = "";
        Argument args = Argument();
        bool _truncated_or_not = true;

        InfGraph(): GraphBase() {
            return;
        }
        InfGraph(string folder_name, string f_name="edgelist_ic.txt", bool format_g=true): GraphBase(folder_name, f_name, true) {
            this->_node_to_RRsets.resize(this->n);
            this->_vec_bool_vis.resize(this->n, false);
            this->_vec_bool_seed.resize(this->n, false);
            return;
        }

        InfGraph(string folder_name, string f_name, string prob_mode): GraphBase(folder_name, f_name, prob_mode) {
            this->_node_to_RRsets.resize(this->n);
            this->_node_to_RRsets_val.resize(this->n);
            this->_node_to_RRsets_delta.resize(this->n);
            this->_vec_bool_vis.resize(this->n, false);
            this->_vec_bool_seed.resize(this->n, false);
            this->_vec_node_to_cand_edge_idx.resize(this->n);
            return;
        }

        void clean_RRsets_InfGraph() {
            VecVecInt32().swap(this->_RRsets); VecVecLargeNum().swap(this->_node_to_RRsets);
            VecVecInt32().swap(this->_RRsets_val); VecVecLargeNum().swap(this->_node_to_RRsets_val);
            VecVecInt32().swap(this->_RRsets_delta); VecVecLargeNum().swap(this->_node_to_RRsets_delta);
            VecBool().swap(this->_vec_bool_RRset_cov_by_S); VecBool().swap(this->_vec_bool_RRset_val_cov_by_S);
            this->_cur_RRsets_num=0; this->_cur_RRsets_num_val=0;
            this->_num_traversed_edges=0;
            this->_node_to_RRsets.resize(this->n);
            this->_node_to_RRsets_val.resize(this->n);
            this->_node_to_RRsets_delta.resize(this->n);
            
            return;
        }

        void set_casc_model(CascadeModel model) {
            this->_casc_model = model;
            return;
        }

        void set_args(Argument arguments) {
            this->args = arguments;
            
            this->_seed_mode = this->args.seed_mode;
            this->_cand_edges_filename = this->folder + "candidate_edges_" + to_string(this->args.num_cand_edges) + "_" + this->_seed_mode + "_" + to_string(this->args.k_seed) + this->_prob_mode + ".txt";
            this->_seed_filename = this->folder + "seeds_" + this->_seed_mode + "_" + to_string(this->args.k_seed) + this->_prob_mode + ".txt";
            
            log_info("--- Arguments ---");
            log_info("k_seed", this->args.k_seed);
            log_info("k_edges", this->args.k_edges);
            log_info("num_cand_edges", this->args.num_cand_edges);
            log_info("epsilon", this->args.epsilon);
            log_info("delta", this->args.delta);
            log_info("seed_mode", this->_seed_mode);
            log_info("random seed", this->args.rand_seed);
            log_info("beta", this->args.beta);
            log_info("num_samples", this->args.num_samples);
            log_info("pruned", this->args.pruned);
            return;
        }

        void set_seed(string seed_mode = "RAND", bool from_file=false) {
            ASSERT(this->_cur_working_folder != "");
            if (from_file)
            {
                string seed_filename = this->_seed_filename;
                log_info("--- Initializing seed nodes from file: " + this->_seed_filename + "---");
                ifstream seed_file(seed_filename);
                uint32_t u;

                while (seed_file >> u) {
                    this->_seed_set_to_augment.push_back(u);
                    this->_vec_bool_seed[u] = true;
                }
                seed_file.close();
                print_vector(this->_seed_set_to_augment);
                return;
            }
            
            if (seed_mode == "RAND") {
                log_info("--- Select the seed set randomly ---");
                for (uint32_t i = 0; i < this->args.k_seed; i++) {
                    uint32_t rand_node = dsfmt_gv_genrand_uint32_range(this->n);
                    // log_info(rand_node);
                    this->_seed_set_to_augment.push_back(rand_node);
                    // log_info(rand_node);
                    this->_vec_bool_seed[rand_node] = true;
                }
            } else if (seed_mode == "IM") {
                if (this->_seed_set_IM.size() == 0)
                {
                    ExitMessage("Please run an IM algorithm first.");
                }
                log_info("--- Select the seed set by IMM ---");
                this->_seed_set_to_augment = NodeList(this->_seed_set_IM.begin(), this->_seed_set_IM.end());
                for (auto& node : this->_seed_set_IM) {
                    this->_vec_bool_seed[node] = true;
                }
            }
            else if (seed_mode == "OUTDEG") {
                log_info("--- Select the seed set by OUTDEG ---");
                std::priority_queue<PairIntSizet, std::vector<PairIntSizet>, CompareBySecond> heap;
                for (size_t i = 0; i < this->n; i++) {
                    heap.push(make_pair(i, this->_G[i].size()));
                }
                for (size_t i = 0; i < this->args.k_seed; i++) {
                    this->_seed_set_to_augment.push_back(heap.top().first);
                    heap.pop();
                }
            }
            print_vector(this->_seed_set_to_augment);
            save_vector(this->_seed_set_to_augment, this->_seed_filename);
            return;
        }

        void init_RRsets_cov_by_S(string RRset_mode="select") {
            if (RRset_mode == "select") {
                ASSERT(this->_cur_RRsets_num > 0);
                VecBool().swap(this->_vec_bool_RRset_cov_by_S);
                this->_vec_bool_RRset_cov_by_S.resize(this->_cur_RRsets_num, false);
                for (auto seed : this->_seed_set_to_augment) {
                    for (auto RRset : this->_node_to_RRsets[seed]) {
                        this->_vec_bool_RRset_cov_by_S[RRset] = true;
                    }
                }
            } else {
                ASSERT(this->_cur_RRsets_num_val > 0);
                VecBool().swap(this->_vec_bool_RRset_val_cov_by_S);
                this->_vec_bool_RRset_val_cov_by_S.resize(this->_cur_RRsets_num_val, false);
                for (auto seed : this->_seed_set_to_augment) {
                    for (auto RRset : this->_node_to_RRsets_val[seed]) {
                        this->_vec_bool_RRset_val_cov_by_S[RRset] = true;
                    }
                }
            }

            return;
        }

        double comp_inf_cov_by_S() {
            ASSERT(this->_vec_bool_RRset_cov_by_S.size() > 0);
            return (1.0 * this->n / this->_cur_RRsets_num) * this->_vec_bool_RRset_cov_by_S.count();
        }

        void init_seed_cov() {
            // this->init_RRsets_cov_by_S();
            if (this->_vec_bool_RRset_cov_by_S.size() == 0) {
                std::cout << "ERROR: Please first run init_RRsets_cov_by_S()" << std::endl;
                exit(1);
            }

            VecLargeNum().swap(this->_marginal_cov);
            std::vector<UVWEdge>().swap(this->_vec_selected_edges);

            this->_marginal_cov.resize(this->n, 0);
            // log_info("The total number of RR sets covered by S is", std::count(this->_vec_bool_RRset_cov_by_S.begin(), this->_vec_bool_RRset_cov_by_S.end(), true));
            for (size_t i = 0; i < this->_RRsets.size(); i++) {
                if (!this->_vec_bool_RRset_cov_by_S[i]) {
                    for (auto x : this->_RRsets[i]) {
                        this->_marginal_cov[x]++;
                    }
                }
            }
            return;
        }
        
        void add_edge(UVWEdge& edge) {
            uint32_t u = edge.first, v = edge.second.first;
            double w = edge.second.second;
            this->_G[u].push_back(Edge(v, w));
            return;
        }

        void pop_edge(UVWEdge& edge) {
            uint32_t u = edge.first;
            this->_G[u].pop_back();
            return;
        }

        void add_edges(std::vector<UVWEdge>& edge_tuple, uint32_t st, uint32_t end) {
            uint32_t u, v;
            double w;
            for (uint32_t i = st; i < end; i++) {
                UVWEdge& e_tuple = edge_tuple[i];
                u = e_tuple.first;
                v = e_tuple.second.first;
                w = e_tuple.second.second;
                this->m++;
                this->_G[u].push_back(Edge(v, w));
                this->_inv_G[v].push_back(Edge(u, w));
                this->_indegree[v]++;
                this->_weighted_in_degree[v] += w;
                this->_weighted_out_degree[u] += w;
            }
            return;
        }

        bool build_one_RRset(uint32_t v, const size_t RRset_idx, const bool flag_count_edges=false, bool fast_trunc_or_not=true) {
            NodeList vec_vis_nodes;         // The list of visited nodes
            size_t num_visit_node = 0, cur_idx = 0;     // Number of activated nodes; Current traversal index
            uint32_t num_trav_edges = 0;
            // For v, the target node
            this->_node_to_RRsets[v].push_back(RRset_idx);
            vec_vis_nodes.push_back(v);
            num_visit_node++;
            this->_vec_bool_vis[v] = true;
            bool end_flag = false;
            bool cov_by_S = false;
            if (this->_vec_bool_seed[v])  // If the target node is a seed node and trunc
            {
                cov_by_S = true;
                this->_vec_bool_RRset_cov_by_S.push_back(true);
                if (fast_trunc_or_not) {
                    end_flag = true;
                }
            }
            while (cur_idx < num_visit_node && !end_flag) {         // when there exist new nodes and we have not encountered any seed node
                const auto expand = vec_vis_nodes[cur_idx++];
                num_trav_edges += this->_inv_G[expand].size();
                if (this->_casc_model == IC) {
                    // if the current node is one-hop activated by the seed set
                    for(auto& nbr: this->_inv_G[expand]) {
                        const auto nbr_id = nbr.first;
                        if (this->_vec_bool_vis[nbr_id]) {
                            // if(num_visit_node == 1) std::cout << "Skip...";
                            continue;
                        }
                        // if the current node is a seed, then skip
                        // if (this->_vec_bool_seed[nbr_id] && this->args.pruned) {
                        //     continue;
                        // }
                        const auto rand_double = dsfmt_gv_genrand_open_close();
                        if (rand_double <= nbr.second) {
                            // The node nbr is activated
                            vec_vis_nodes.push_back(nbr_id);
                            num_visit_node++;
                            this->_vec_bool_vis[nbr_id] = true;
                            this->_node_to_RRsets[nbr_id].push_back(RRset_idx);

                            if (this->_vec_bool_seed[nbr_id]) {
                                cov_by_S = true;
                                this->_vec_bool_RRset_cov_by_S.push_back(true);
                                if (fast_trunc_or_not) {
                                    end_flag = true;
                                    break;
                                }
                            }
                        }
                    }
                }
                else if (this->_casc_model == LT) {
                    if (this->_inv_G[expand].size() == 0) break;
                    num_trav_edges += this->_inv_G[expand].size();
                    const auto next_nbr_idx = gen_rand_node_by_weight_LT(this->_inv_G[expand]);
                    if (next_nbr_idx >= this->_inv_G[expand].size()) break; 
                    const auto nbr_id = this->_inv_G[expand][next_nbr_idx].first;
                    if (this->_vec_bool_vis[nbr_id]) break;

                    vec_vis_nodes.push_back(nbr_id);
                    num_visit_node++;
                    this->_vec_bool_vis[nbr_id] = true;
                    this->_node_to_RRsets[nbr_id].push_back(RRset_idx);
                }
            }
            
            ASSERT(num_visit_node == (size_t)vec_vis_nodes.size());
            
            for (auto i=0; i<num_visit_node; i++) this->_vec_bool_vis[vec_vis_nodes[i]] = false;
            this->_RRsets.push_back(vec_vis_nodes);
            if (!cov_by_S) this->_vec_bool_RRset_cov_by_S.push_back(false);
            
            // NodeList().swap(vec_vis_nodes);
            if (flag_count_edges) {
                this->_num_traversed_edges += num_trav_edges;
            }
            return cov_by_S;
        }

        bool build_one_RRset_val(uint32_t v, const size_t RRset_idx, const bool flag_count_edges=false, bool fast_trunc_or_not=true) {
            NodeList vec_vis_nodes;         // The list of visited nodes
            size_t num_visit_node = 0, cur_idx = 0;     // Number of activated nodes; Current traversal index
            uint32_t num_trav_edges = 0;
            // For v, the target node
            this->_node_to_RRsets_val[v].push_back(RRset_idx);
            vec_vis_nodes.push_back(v);
            num_visit_node++;
            this->_vec_bool_vis[v] = true;
            bool end_flag = false;
            bool cov_by_S = false;
            if (this->_vec_bool_seed[v])  // If the target node is a seed node and
            {
                cov_by_S = true;
                this->_vec_bool_RRset_val_cov_by_S.push_back(true);
                if (fast_trunc_or_not) {
                    end_flag = true;
                }
            }

            while (cur_idx < num_visit_node && !end_flag) {         // when there exist new nodes and we have not encountered any seed node
                const auto expand = vec_vis_nodes[cur_idx++];
                num_trav_edges += this->_inv_G[expand].size();
                if (this->_casc_model == IC) {
                    for(auto& nbr: this->_inv_G[expand]) {
                        const auto nbr_id = nbr.first;
                        if (this->_vec_bool_vis[nbr_id]) {
                            // if(num_visit_node == 1) std::cout << "Skip...";
                            continue;
                        }
                        const auto rand_double = dsfmt_gv_genrand_open_close();
                        if (rand_double <= nbr.second) {
                            // The node nbr is activated
                            vec_vis_nodes.push_back(nbr_id);
                            num_visit_node++;
                            this->_vec_bool_vis[nbr_id] = true;
                            this->_node_to_RRsets_val[nbr_id].push_back(RRset_idx);
                            if (this->_vec_bool_seed[nbr_id]) {
                                cov_by_S = true;
                                this->_vec_bool_RRset_val_cov_by_S.push_back(true);
                                if (fast_trunc_or_not) {
                                    end_flag = true;
                                    break;
                                }
                            }
                        }
                    }
                }
                else if (this->_casc_model == LT) {
                    if (this->_inv_G[expand].size() == 0) break;
                    num_trav_edges += this->_inv_G[expand].size();
                    const auto next_nbr_idx = gen_rand_node_by_weight_LT(this->_inv_G[expand]);
                    if (next_nbr_idx >= this->_inv_G[expand].size()) break; 
                    const auto nbr_id = this->_inv_G[expand][next_nbr_idx].first;
                    if (this->_vec_bool_vis[nbr_id]) break;

                    vec_vis_nodes.push_back(nbr_id);
                    num_visit_node++;
                    this->_vec_bool_vis[nbr_id] = true;
                    this->_node_to_RRsets_val[nbr_id].push_back(RRset_idx);
                }
            }
            
            ASSERT(num_visit_node == (size_t)vec_vis_nodes.size());
            
            for (auto i=0; i<num_visit_node; i++) this->_vec_bool_vis[vec_vis_nodes[i]] = false;
            this->_RRsets_val.push_back(vec_vis_nodes);
            if (!cov_by_S) this->_vec_bool_RRset_val_cov_by_S.push_back(false);
            // NodeList().swap(vec_vis_nodes);
            
            if (flag_count_edges) {
                this->_num_traversed_edges += num_trav_edges;
            }
            return cov_by_S;
        }

        void build_RRsets(size_t num_RRsets, bool fast_truncated=false, string RRsets_mode="select") {
            if (num_RRsets > SIZE_MAX) {
                std::cout << "Error: The number of RR sets is too large." << std::endl;
                exit(1);   
            }
            if (RRsets_mode == "select") {
                const auto prev_size = this->_cur_RRsets_num;
                this->_cur_RRsets_num = std::max(this->_cur_RRsets_num, num_RRsets);
                std::cout << "current RR sets number: " << this->_cur_RRsets_num << std::endl;
                for (auto i = prev_size; i < num_RRsets; i++) {
                    uint32_t rand_node = dsfmt_gv_genrand_uint32_range(this->n);
                    // this->_node_sampled_times[rand_node]++;
                    this->build_one_RRset(rand_node, i, false, fast_truncated);
                }
                ASSERT(this->_cur_RRsets_num == this->_RRsets.size());
            } else {
                const auto prev_size = this->_cur_RRsets_num_val;
                this->_cur_RRsets_num_val = std::max(this->_cur_RRsets_num_val, num_RRsets);
                std::cout << "current RR sets number: " << this->_cur_RRsets_num_val << std::endl;
                for (auto i=prev_size; i<num_RRsets; i++) {
                    uint32_t rand_node = dsfmt_gv_genrand_uint32_range(this->n);
                    // this->_node_sampled_times[rand_node]++;
                    this->build_one_RRset_val(rand_node, i, false, fast_truncated);
                }
                ASSERT(this->_cur_RRsets_num_val == this->_RRsets_val.size());
            }
            return;
        }
        
        double comp_inf_by_cov(const NodeList& vec_seed) {
            VecBool vec_bool_vis(this->_cur_RRsets_num);
            for (auto seed: vec_seed) {
                for(auto RRset: this->_node_to_RRsets[seed]) {
                    vec_bool_vis[RRset] = true;
                }
            }
            return 1.0 * vec_bool_vis.count() * this->n / this->_cur_RRsets_num;
        }

        double edge_selection_on_edge_RRsets() {
            VecDouble live_prob_per_RRset(this->_RRsets_delta.size(), 1.0);    // store the live probability of every RR set
            VecLargeNum vec_selected_edges_idx;
            // Use a heap to store the <node, prob coverage> pair
            std::priority_queue<PairSizetDouble, std::vector<PairSizetDouble>, CompareBySecondDoubleForSizet> heap;
            VecDouble coverage(this->_vec_cand_edges.size(), 0.0);       // coverage[v] is the expected number of RR sets covered by v

            for(size_t i=0; i<this->_vec_cand_edges.size(); i++) {
                // store coverage
                // if the node is not covered 
                double prob_for_i = this->_vec_cand_edges[i].second.second;
                uint32_t cur_node = this->_vec_cand_edges[i].second.first;
                double edge_i_cov = prob_for_i * this->_node_to_RRsets_delta[cur_node].size();
                PairSizetDouble tmp(make_pair(i, edge_i_cov));
                heap.push(tmp);
                coverage[i] = edge_i_cov;
            }
            log_info("Finish initializing marginal coverage, start selection");

            this->_vec_selected_edges.clear();
            size_t top_edge_idx;
            uint32_t max_node_idx;
            double cov_num = 0.0;
            while (this->_vec_selected_edges.size() < this->args.k_edges) {
                PairSizetDouble top = heap.top();
                heap.pop();

                // Lazy Update
                if (top.second > coverage[top.first]) {
                    // Update coverage of top
                    top.second = coverage[top.first];
                    heap.push(top);
                    continue;
                }
                top_edge_idx = top.first;
                cov_num += coverage[top_edge_idx];
                double origin_cov = coverage[top_edge_idx];
                
                // selection
                max_node_idx = this->_vec_cand_edges[top_edge_idx].second.first;   // the index of the edge's target node
                double top_prob = this->_vec_cand_edges[top_edge_idx].second.second;      // the probability of the target edge
                this->_vec_selected_edges.push_back(this->_vec_cand_edges[top_edge_idx]);
                vec_selected_edges_idx.push_back(top_edge_idx);

                // After selecting one node, we need to reduce marginal gain of every node
                // e: the RR sets covered by the node max_idx
                VecLargeNum& RRsets_cov_by_max_node = this->_node_to_RRsets_delta[max_node_idx];
                for (uint32_t j=0; j<RRsets_cov_by_max_node.size(); j++) {
                    // Process one RR set covered by the node max_node_idx
                    size_t RRset_idx = RRsets_cov_by_max_node[j];
                    NodeList& node_list = this->_RRsets_delta[RRset_idx];
                    for (uint32_t& node: node_list){
                        for (auto& edge: this->_vec_node_prob_vec[node])
                        {
                            size_t e_idx = edge.first;
                            if (e_idx == top_edge_idx)
                            {
                                continue;
                            }
                            double cur_prob = edge.second;
                            coverage[e_idx] -= (live_prob_per_RRset[RRset_idx] * top_prob * cur_prob);
                        }
                    }
                    // update RR set live prob
                    live_prob_per_RRset[RRset_idx] *= (1 - top_prob);
                }
                
                coverage[top_edge_idx] = 0.0;
                //update edges pointing to the same node
                // if (!this->_vec_node_prob_heap[max_idx].empty()) {
                //     // for (uint32_t j=0; j<RRsets_cov_by_max_node.size(); j++) {
                //     //     size_t RRset_idx = RRsets_cov_by_max_node[j];
                //     //     coverage[max_idx] += live_prob_per_RRset[RRset_idx];
                //     // }
                    
                //     double new_prob = this->_vec_node_prob_heap[max_idx].top().second;
                //     coverage[max_idx] = origin_cov / top_prob * new_prob * (1-top_prob);
                //     // coverage[max_idx] *= new_prob;
                //     heap.push(make_pair(max_idx, coverage[max_idx]));
                // }
                
            }
            log_info("cov_num", cov_num);
            return cov_num;
        }
        double edge_selection_greedy() {
            VecDouble live_prob_per_RRset(this->_RRsets_delta.size(), 1.0);    // store the live probability of every RR set
            VecLargeNum vec_selected_edges_idx;
            // Use a heap to store the <node, prob coverage> pair
            std::priority_queue<PairIntDouble, std::vector<PairIntDouble>, CompareBySecondDouble> heap;
            VecDouble coverage(this->n, 0.0);       // coverage[v] is the expected number of RR sets covered by v
            // Timer timer("Test");
            
            for(uint32_t i=0; i<this->n; i++) {
                // store coverage
                // if the node is not covered 
                if (this->_vec_node_prob_heap[i].empty())
                    continue;
                double prob_for_i = this->_vec_node_prob_heap[i].top().second;
                double node_i_cov = prob_for_i * this->_node_to_RRsets_delta[i].size();
                PairIntDouble tmp(make_pair(i, node_i_cov));
                heap.push(tmp);
                coverage[i] = node_i_cov;
            }
            log_info("Finish initializing marginal coverage, start selection");
            
            this->_vec_selected_edges.clear();
            uint32_t max_idx;
            double cov_num = 0.0;
            double update_time = 0.0;
            while (this->_vec_selected_edges.size() < this->args.k_edges) {
                PairIntDouble top = heap.top();
                heap.pop();

                // Lazy Update
                if (top.second > coverage[top.first]) {
                    // Update coverage of top
                    top.second = coverage[top.first];
                    heap.push(top);
                    continue;
                }
                max_idx = top.first;
                cov_num += coverage[max_idx];
                double origin_cov = coverage[max_idx];
                
                // selection
                size_t top_edge_idx = this->_vec_node_prob_heap[max_idx].top().first;   // the index of the target edge
                double top_prob = this->_vec_node_prob_heap[max_idx].top().second;      // the probability of the target edge

                this->_vec_node_prob_heap[max_idx].pop();       // pop the edge for node max_idx
                this->_vec_selected_edges.push_back(this->_vec_cand_edges[top_edge_idx]);
                vec_selected_edges_idx.push_back(top_edge_idx);

                // After selecting one node, we need to reduce marginal gain of every node
                // e: the RR sets covered by the node max_idx
                VecLargeNum& RRsets_cov_by_max_node = this->_node_to_RRsets_delta[max_idx];
                for (uint32_t j=0; j<RRsets_cov_by_max_node.size(); j++) {
                    size_t RRset_idx = RRsets_cov_by_max_node[j];
                    NodeList node_list = this->_RRsets_delta[RRset_idx];
                    // timer.get_operation_time();
                    for (uint32_t& node: node_list){
                        if (node != max_idx && !this->_vec_node_prob_heap[node].empty()) {
                            double cur_node_prob = this->_vec_node_prob_heap[node].top().second;
                            coverage[node] -= (live_prob_per_RRset[RRset_idx] * top_prob * cur_node_prob);
                        }
                    }

                    // update RR set live prob
                    live_prob_per_RRset[RRset_idx] *= (1 - top_prob);
                }
                
                //update edges pointing to the same node
                coverage[max_idx] = 0.0;
                if (!this->_vec_node_prob_heap[max_idx].empty()) {
                    double new_prob = this->_vec_node_prob_heap[max_idx].top().second;
                    
                    coverage[max_idx] = origin_cov / top_prob * new_prob * (1-top_prob);
                    heap.push(make_pair(max_idx, coverage[max_idx]));
                }
            }
            
            // Post-process
            for (auto& idx : vec_selected_edges_idx) {
                auto& uvw_edge = this->_vec_cand_edges[idx];
                uint32_t targ_node = uvw_edge.second.first;
                double p_e = uvw_edge.second.second;
                this->_vec_node_prob_heap[targ_node].push(make_pair(idx, p_e));
            }

            log_info("cov_num", cov_num);
            return cov_num;
        }

        double edge_selection_by_max_cover(VecVecLargeNum& vec_u, VecVecLargeNum& vec_v, VecLargeNum& res) {
            size_t num_u = vec_u.size();
            size_t num_v = vec_v.size();
            uint32_t target_size = this->args.k_edges;

            // Use a heap to store the <node, coverage> pair
            std::priority_queue<PairSizet, std::vector<PairSizet>, CompareBySecond> heap;
            VecLargeNum coverage(num_u, 0);       // coverage[v] is the RR sets covered by v
            int cnt = 0;
            for(uint32_t i=0; i<num_u; i++) {
                // store coverage
                size_t node_i_cov = vec_u[i].size();
                PairSizet tmp(make_pair(i, node_i_cov));
                heap.push(tmp);
                coverage[i] = node_i_cov;
                
            }

            VecBool RRsets_mark(vec_v.size(), false);
            VecBool node_mark(num_u, false);   

            uint32_t max_idx;
            size_t cov_num = 0;
            while (res.size() < target_size) {
                PairSizet top = heap.top();
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
            
            // Update selected edges
            std::vector<UVWEdge>().swap(this->_vec_selected_edges);
            for (auto i : res) {
                this->_vec_selected_edges.push_back(this->_vec_cand_edges[i]);
            }
            return (double)cov_num;
        }

        // Link Recommendation
        void generate_candidate_edges(size_t num) {
            ASSERT(this->_seed_set_to_augment.size() > 0);
            vector<UVWEdge> res;
            VecSetInt added_nodes(this->args.k_seed);
            VecBool vis_by_v(this->n, false);
            for (uint32_t i = 0; i < this->_seed_set_to_augment.size(); i++) {
                for (auto& node : this->_G[this->_seed_set_to_augment[i]]) {
                    added_nodes[i].insert(node.first);
                }
            }

            for (size_t i = 0; i < num; i++) {
                uint32_t rand_idx = dsfmt_gv_genrand_uint32_range(this->_seed_set_to_augment.size());
                uint32_t rand_node_S = this->_seed_set_to_augment[rand_idx];
                
                uint32_t rand_node = dsfmt_gv_genrand_uint32_range(this->n);
                int cnt = 0;
                while (this->_vec_bool_seed[rand_node] || added_nodes[rand_idx].count(rand_node) != 0)
                {
                    rand_idx = dsfmt_gv_genrand_uint32_range(this->_seed_set_to_augment.size());
                    rand_node_S = this->_seed_set_to_augment[rand_idx];
                    rand_node = dsfmt_gv_genrand_uint32_range(this->n);
                    cnt++;
                    if (cnt > 1000) {
                        log_info("--- FAILED to generate candidate edges in 1000 iterations ---");
                        return;
                    }
                }
                
                double w = 0.0;
                if (this->_casc_model == IC) {
                    uint32_t targ_indeg = this->_inv_G[rand_node].size(), st_outdeg = this->_G[rand_node_S].size();
                    double avg_w_in, avg_w_out;
                    if (targ_indeg == 0)
                        avg_w_in = dsfmt_gv_genrand_open_close();
                    if (st_outdeg == 0) avg_w_out = dsfmt_gv_genrand_open_close();
                    if (targ_indeg != 0 && st_outdeg != 0) {
                        avg_w_in = this->_weighted_in_degree[rand_node] / this->_inv_G[rand_node].size();
                        avg_w_out = this->_weighted_out_degree[rand_node_S] / this->_G[rand_node_S].size();
                    }

                    w = (avg_w_in + avg_w_out) / 2;
                } else
                    w = std::max(1 - this->_weighted_in_degree[rand_node], 0.0);
                
                res.push_back(make_pair(rand_node_S, Edge(rand_node, w)));
                added_nodes[rand_idx].insert(rand_node);
            }
            log_info("--- Start writing candidate edges ---");
            this->write_cand_edges(res);
            this->_vec_cand_edges = res;
            return;
        }

        void generate_candidate_edges() {
            ASSERT(this->_seed_set_to_augment.size() > 0);
            vector<UVWEdge> res;
            VecSetInt added_nodes(this->args.k_seed);

            for (uint32_t i = 0; i < this->_seed_set_to_augment.size(); i++) {
                added_nodes[i].insert(this->_seed_set_to_augment[i]);
                for (auto& node : this->_G[this->_seed_set_to_augment[i]]) {
                    added_nodes[i].insert(node.first);
                }
            }
            
            for (size_t i = 0; i < this->n; i++) {
                for (size_t j=0; j < this->args.k_seed; j++) {
                    uint32_t seed = this->_seed_set_to_augment[j];
                    double w = 0.0;
                    if (this->_casc_model == IC) {
                        uint32_t targ_indeg = this->_inv_G[i].size(), st_outdeg = this->_G[seed].size();
                        double avg_w_in, avg_w_out;
                        
                        if (targ_indeg == 0)
                            avg_w_in = dsfmt_gv_genrand_open_close();
                        if (st_outdeg == 0) avg_w_out = dsfmt_gv_genrand_open_close();
                        if (targ_indeg != 0 && st_outdeg != 0) {
                            avg_w_in = this->_weighted_in_degree[i] / targ_indeg;
                            avg_w_out = this->_weighted_out_degree[seed] / st_outdeg;
                        }
                        w = (avg_w_in + avg_w_out) / 2;
                        if (w > 1.0) {
                            ExitMessage("probability over 1, please check");
                        }

                    } else
                        w = std::max(1 - this->_weighted_in_degree[i], 0.0);
                    if (this->args.probability_mode == "0001" ) {
                        w = 0.001;
                    } else if (this->args.probability_mode == "001") {
                        w = 0.01;
                    }
                    if (w > 1.0) {
                        log_info("Target's weighted indeg", this->_weighted_in_degree[i]);
                        log_info("Target's indeg", this->_inv_G[i].size());
                        log_info("Seed's weighted outdeg", this->_weighted_out_degree[seed]);
                        log_info("Seed's outdeg", this->_G[seed].size());
                        ASSERT(w <= 1.0);
                    }
                    if (!this->_vec_bool_seed[i] && added_nodes[j].count(i) == 0)
                    {
                        
                        res.push_back(make_pair(seed, Edge(i, w)));
                        added_nodes[j].insert(i);
                    }
                }
            }
            
            log_info("--- Writing candidate edges to " + this->_cand_edges_filename+".vec" + " ---");
            // this->write_cand_edges(res);
            this->_vec_cand_edges = res;
            TIO::save_file(this->_cand_edges_filename+".vec", this->_vec_cand_edges);
            return;
        }

        void write_cand_edges(vector<UVWEdge> edge_tuple) {
            ASSERT(this->_cur_working_folder != "");
            ofstream save_file(this->_cand_edges_filename);
            uint32_t st, end;
            double w;
            for(auto& tup : edge_tuple) {
                st = tup.first; end = tup.second.first; w = tup.second.second;
                save_file << st << '\t' << end << '\t' << w << '\n';
            }
            save_file.close();
            log_info("--- Finish saving the candidate edges. ---");
            return;
        }

        void read_cand_edges() {
            if (this->_vec_cand_edges.size() > 0) {
                log_info("--- Candidate edges already loaded ---");
                return;
            }
            // std::vector<UVWEdge>().swap(this->_vec_cand_edges);
            ASSERT(this->_vec_cand_edges.size() == 0);
            
            string filename = this->_cand_edges_filename + ".vec";
            if (check_file_exist(filename))
            {
                log_info("Loading Candidate edges from", filename);
                TIO::load_file(filename, this->_vec_cand_edges);
                this->args.num_cand_edges = this->_vec_cand_edges.size();
                
            }
            else {
                filename = this->_cand_edges_filename;
                ifstream myfile(filename);
                uint32_t u, v;
                double w;
                while (myfile >> u >> v >> w) {
                    this->_vec_cand_edges.push_back(make_pair(u, Edge(v, w)));
                }
                myfile.close();
                TIO::save_file(this->_cand_edges_filename + ".vec", this->_vec_cand_edges);
                this->args.num_cand_edges = this->_vec_cand_edges.size();
            }
            
            return;
        }

        // update the coverage
        void update_RRsets_cov_by_S(std::vector<UVWEdge>& selected_edges) {
            ASSERT(this->_RRsets.size() > 0);
            for (auto& e_tuple: selected_edges)
            {
                // update coverage after selecting
                uint32_t target_node = e_tuple.second.first;
                double w = e_tuple.second.second;
                const VecLargeNum& target_RRsets = this->_node_to_RRsets[target_node];
                
                for (size_t RRset_id : target_RRsets) {
                    double rand_double = dsfmt_gv_genrand_open_close();
                    bool flag_add = false;
                    if (this->_casc_model == IC) flag_add = (rand_double < w);
                    else flag_add = (this->_RRsets[RRset_id].back() == target_node && rand_double < 1 - this->_weighted_in_degree[target_node]);
                    if (flag_add && !this->_vec_bool_RRset_cov_by_S[RRset_id]) {
                        this->_vec_bool_RRset_cov_by_S[RRset_id] = true;
                    } else continue;
                }
            }
            return;
        }

        double comp_cov_pmc(std::vector<UVWEdge>& selected_edges, string RRset_mode="val") {
            VecVecInt32& cur_RRsets = (RRset_mode == "select") ? this->_RRsets : this->_RRsets_val;
            VecVecLargeNum& cur_node_to_RRsets = (RRset_mode == "select")? this->_node_to_RRsets : this->_node_to_RRsets_val;
            VecBool& vec_bool_RRset_cov_by_S = (RRset_mode == "select")? this->_vec_bool_RRset_cov_by_S : this->_vec_bool_RRset_val_cov_by_S;
            ASSERT(cur_RRsets.size() > 0);

            ASSERT(vec_bool_RRset_cov_by_S.size() == cur_RRsets.size());
            // this->init_RRsets_cov_by_S(RRset_mode);
            double cov_num = 0.0;
            VecBool included_in_A (this->n, false);
            std::unordered_map<uint32_t, std::vector<double>> node_probs_map;
            for (auto& e_tuple: selected_edges) {
                uint32_t target_node = e_tuple.second.first;
                double w = e_tuple.second.second;
                node_probs_map[target_node].push_back(w);
                included_in_A[target_node] = true;
            }

            for (size_t i=0; i<cur_RRsets.size(); i++) {
                if (vec_bool_RRset_cov_by_S[i]) {
                    continue;
                }
                double prob_not_cov = 1.0;
                for (uint32_t node : cur_RRsets[i]) {
                    if (included_in_A[node]) {
                        for (double p : node_probs_map[node]) {
                            prob_not_cov *= (1.0 - p);
                        }
                    }
                }

                cov_num += (double)(1.0 - prob_not_cov);
            }
            
            printf("cov-num:%.3f\n", cov_num);
            return cov_num;
        } 

        void init_from_files(bool load_cand_edges=true) {
            ASSERT(this->_cur_working_folder != "");
            log_info("--- Initializing from previous files ---");

            string cand_edges_file_name = this->_cand_edges_filename;
            string seed_filename = this->_seed_filename;

            uint32_t u, v;
            double w;
            if (check_file_exist(seed_filename)) {
                log_info("--- File of seed nodes: " + seed_filename + " ---");
            } else {
                log_info("--- VITAL: Please select the seeds and store it in " + seed_filename + " ---");
                return;
            }
            // If we need to load all candidate edges, please uncomment the following loop
            if (load_cand_edges) {
                if (check_file_exist(cand_edges_file_name)) {
                    log_info("--- File of candidate edges: " + cand_edges_file_name + " ---");
                    ifstream cand_edges_file(cand_edges_file_name);
                    while (cand_edges_file >> u >> v >> w) {
                        this->_vec_cand_edges.push_back(make_pair(u, Edge(v, w)));
                    }
                } else {
                    log_info("--- VITAL: Please generate the candidate edges set and store it in " + cand_edges_file_name + " ---");
                    return;
                }
            }

            ifstream seed_file(seed_filename);            

            while (seed_file >> u) {
                this->_seed_set_to_augment.push_back(u);
                this->_vec_bool_seed[u] = true;
            }
            seed_file.close();
            
            return;
        }

        void organize_edges(string mode="heap") {
            assert(this->_vec_cand_edges.size()>0);
            if (mode == "heap")
            {
                this->_vec_node_prob_heap.resize(this->n);
                // Traverse cand_edges and organize probability
                for (size_t i = 0; i < this->_vec_cand_edges.size(); i++) {
                    uint32_t v = this->_vec_cand_edges[i].second.first;
                    double puv = this->_vec_cand_edges[i].second.second;
                    this->_vec_node_prob_heap[v].push(make_pair(i, puv));
                }
                log_info("Finish organizing prob into a heap");
            }
            else
            {
                this->_vec_node_prob_vec.resize(this->n);
                // Traverse cand_edges and organize probability
                for (size_t i = 0; i < this->_vec_cand_edges.size(); i++) {
                    uint32_t v = this->_vec_cand_edges[i].second.first;
                    double puv = this->_vec_cand_edges[i].second.second;
                    this->_vec_node_prob_vec[v].push_back(make_pair(i, puv));
                }
                log_info("Finish organizing prob into a vec");
            }
            
            return;
        }

        void fixed_ScaLIM(string log_filepath="") {
            ASSERT(this->_seed_set_to_augment.size() > 0);
            // preparation
            ResultInfo res;
            double sampling_time = 0.0, selection_time = 0.0;
            
            log_filepath = this->_cur_working_folder + "log_ScaLIM.txt";
            
            ofstream log_file(log_filepath, std::ios::app);
            
            log_info("--- Start Reading Candidate Edges ---");
            this->read_cand_edges();
            
            Timer timer("ScaLIM");
            
            log_info("--- Generating RR sets ---");
            this->build_RRsets(this->args.num_samples, true);
            this->update_RRsets_delta();
            log_info("RR sets delta", this->_RRsets_delta.size());
            sampling_time = timer.get_operation_time();
            log_info("Time for Sampling RR sets", sampling_time);

            log_info("--- Initializing coverage boolean arrays ---");
            
            this->init_seed_cov();
            log_file << "Original Inf: " << this->comp_inf_cov_by_S() << '\n';

            // Traverse cand_edges and organize probability
            // this->_vec_node_prob_heap.resize(this->n);
            // for (size_t i = 0; i < this->_vec_cand_edges.size(); i++) {
            //     uint32_t v = this->_vec_cand_edges[i].second.first;
            //     double puv = this->_vec_cand_edges[i].second.second;
            //     this->_vec_node_prob_heap[v].push(make_pair(i, puv));
            // }
            // log_info("Finish organizing prob");

            this->organize_edges("heap");
            
            // optimize
            double res_inf = this->edge_selection_greedy();
            selection_time = timer.get_operation_time();

            // write resulting edges
            string edges_save_path = this->_cur_working_folder+"selected_edges_ScaLIM.txt";
            log_info("--- Start writing edges ---");
            write_UVWEdges(this->_vec_selected_edges, edges_save_path);

            // Setting Results
            res.set_k_edges(this->args.k_edges);
            res.set_influence_original(res_inf * this->n / this->_cur_RRsets_num);
            res.set_sampling_time(sampling_time);
            res.set_selection_time(selection_time);
            res.set_RR_sets_size(this->_cur_RRsets_num);
            res.save_to_file(log_filepath);

            return;
        }

        void fixed_ScaLIM_minus(string log_filepath="") {
            ASSERT(this->_seed_set_to_augment.size() > 0);
            // preparation
            ResultInfo res;
            double sampling_time = 0.0, selection_time = 0.0;
            
            log_filepath = this->_cur_working_folder + "log_ScaLIM_minus.txt";
            
            ofstream log_file(log_filepath, std::ios::app);
            
            log_info("--- Start Reading Candidate Edges ---");
            this->read_cand_edges();
            
            Timer timer("ScaLIM_minus");
            
            log_info("--- Generating RR sets ---");
            this->build_RRsets(this->args.num_samples, true);
            this->update_RRsets_delta();
            sampling_time = timer.get_operation_time();
            log_info("Time for Sampling RR sets", sampling_time);

            log_info("--- Initializing coverage boolean arrays ---");
            
            this->init_seed_cov();
            log_file << "Original Inf: " << this->comp_inf_cov_by_S() << '\n';

            this->organize_edges("vec");
            
            // optimize PMC
            double res_inf = this->edge_selection_on_edge_RRsets();
            selection_time = timer.get_operation_time();

            // write resulting edges
            string edges_save_path = this->_cur_working_folder+"selected_edges_ScaLIM_minus.txt";
            write_UVWEdges(this->_vec_selected_edges, edges_save_path);

            // Setting Results
            res.set_k_edges(this->args.k_edges);
            res.set_influence_original(res_inf);
            res.set_sampling_time(sampling_time);
            res.set_selection_time(selection_time);
            res.set_RR_sets_size(this->_cur_RRsets_num);
            res.save_to_file(log_filepath);

            return;
        }

        void fixed_Greedy2(string log_filepath="") {
            ASSERT(this->_seed_set_to_augment.size() > 0);
            ASSERT(this->_vec_cand_edges.size() > 0);
            
            log_filepath = this->_cur_working_folder + "log_Greedy2.txt";
            ResultInfo res;
            double sampling_time = 0.0, sampling_2_time=0.0, selection_time = 0.0;
            ofstream log_file(log_filepath, std::ios::app);
            
            Timer timer("Greedy2");
            log_info("--- Sampling RR sets ---");
            this->build_RRsets(this->args.num_samples, true);
            sampling_time = timer.get_operation_time();
            log_info("Original Inf", this->comp_inf_cov_by_S());
            
            VecVecLargeNum RRsets_for_edges(this->_cur_RRsets_num);
            VecVecLargeNum edge_to_RRsets(this->_vec_cand_edges.size());

            // For every edge, initialize its RR sets. Re-organize RR sets for edges
            int cnt = 0, all_cnt = 0;
            for (size_t i = 0; i < this->_vec_cand_edges.size(); i++) {
                auto& cand_edge = this->_vec_cand_edges[i];
                uint32_t u = cand_edge.first, v = cand_edge.second.first;
                double puv = cand_edge.second.second;
                for (size_t RR_idx : this->_node_to_RRsets[v]) {
                    if (!this->_vec_bool_RRset_cov_by_S[RR_idx] && dsfmt_gv_genrand_open_close() < puv) {
                        RRsets_for_edges[RR_idx].push_back(i);
                        edge_to_RRsets[i].push_back(RR_idx);
                    }
                }
            }
            sampling_2_time += timer.get_operation_time();
            
            // Selection
            // VecLargeNum vec_selected_idx = max_cover_by_heap(edge_to_RRsets, RRsets_for_edges, this->args.k_edges);
            log_info("--- Start Edge Selection ---");
            VecLargeNum vec_selected_idx = max_cover_by_heap(edge_to_RRsets, RRsets_for_edges, this->args.k_edges);
            this->_vec_selected_edges.clear();
            for (auto i : vec_selected_idx) {
                this->_vec_selected_edges.push_back(this->_vec_cand_edges[i]);
            }
            selection_time = timer.get_operation_time();
            double inf_self = (comp_cov(edge_to_RRsets, RRsets_for_edges, vec_selected_idx) + this->_vec_bool_RRset_cov_by_S.count()) * this->n / this->_cur_RRsets_num;
            
            string edges_save_path = this->_cur_working_folder+"selected_edges_Greedy2.txt";
            write_UVWEdges(this->_vec_selected_edges, edges_save_path);
            
            // Setting Results
            res.set_k_edges(this->args.k_edges);
            res.set_influence_original(inf_self);
            res.set_sampling_time(sampling_time);
            res.set_sampling_2_time(sampling_2_time);
            res.set_selection_time(selection_time);
            res.set_RR_sets_size(this->_cur_RRsets_num);
            res.save_to_file(log_filepath);

            return;
        }

        void sample_edges(VecVecLargeNum& edge_to_RRsets, VecVecLargeNum& RRsets_for_edges, string RRset_mode="select") {
            assert(this->_vec_node_prob_vec.size()>0);
            VecVecInt32& cur_RRsets = (RRset_mode == "select") ? this->_RRsets : this->_RRsets_val;
            VecVecLargeNum& cur_node_to_RRsets = (RRset_mode == "select")? this->_node_to_RRsets : this->_node_to_RRsets_val;
            VecBool& vec_bool_cov_by_S = (RRset_mode == "select")? this->_vec_bool_RRset_cov_by_S : this->_vec_bool_RRset_val_cov_by_S;

            size_t RRsets_num = cur_RRsets.size();

            if (RRsets_for_edges.size() == RRsets_num) {
                return;
            }

            size_t edge_RRsets_num = RRsets_for_edges.size();

            for (size_t i = edge_RRsets_num; i < RRsets_num; i++) {
                RRsets_for_edges.push_back(VecLargeNum());
                if (vec_bool_cov_by_S[i]) {
                    continue;
                }

                VecuInt32& RRset = cur_RRsets[i];
                for (auto node : RRset) {
                    for (auto& idx_prob_pair : this->_vec_node_prob_vec[node]) {
                        size_t e_idx = idx_prob_pair.first;
                        double e_prob = idx_prob_pair.second;
                        if (dsfmt_gv_genrand_open_close() < e_prob) {
                            RRsets_for_edges[i].push_back(e_idx);
                            edge_to_RRsets[e_idx].push_back(i);
                        }
                    }
                }
            }
            return;
        }

        void update_RRsets_delta() {
            size_t prev_size = this->_RRsets_delta.size();
            size_t cur_idx = prev_size;
            // Re-index the RR sets not covered by S
            for (size_t i = prev_size; i < this->_RRsets.size(); i++)
            {
                if (this->_vec_bool_RRset_cov_by_S[i]) {
                    continue;
                }

                auto& RRset = this->_RRsets[i];
                VecuInt32 tmp_RRset;
                for(auto node: RRset) {
                    tmp_RRset.push_back(node);      // current RR set covers node
                    this->_node_to_RRsets_delta[node].push_back(cur_idx);
                }
                this->_RRsets_delta.push_back(tmp_RRset);
                cur_idx++;
            }
            ASSERT(cur_idx == this->_RRsets_delta.size());
            return;
        }

        void OPIM_ScaLIM() {
            ASSERT(this->_seed_set_to_augment.size() > 0);
            ASSERT(this->_vec_cand_edges.size() > 0);
            
            // Initialize parameters
            const double delta = this->args.delta;
            const double e = exp(1);
            const double approx = 1 - 1.0 / e;
            const double alpha = sqrt(log(6.0 / delta));
            const double beta = sqrt((1 - 1 / e) * (logcnk(this->n*this->args.k_seed, this->args.k_edges) + log(6.0 / delta)));
            const auto numRbase = size_t(2.0 * pow2((1 - 1 / e) * alpha + beta));
            const auto maxNumR = size_t(2.0 * this->n * pow2((1 - 1 / e) * alpha + beta) / this->args.k_edges / pow2(this->args.epsilon)) + 1;
            const auto numIter = (size_t)log2(maxNumR / numRbase) + 1;
            const double a1 = log(numIter * 3.0 / delta);
            const double a2 = log(numIter * 3.0 / delta);
            double min_u_bound = __DBL_MAX__;
            double last_bound = __DBL_MAX__;
            double sampling_time= 0.0, selection_time = 0.0;
            ResultInfo res;
            string log_filepath = this->_cur_working_folder + "log_ScaLIM.txt";
            Timer timer_ScaLIM("OPIM_ScaLIM");
            
            this->organize_edges("heap");
            timer_ScaLIM.log_operation_time("organizing");
            double last_approxC = 0.0;

            for (auto idx = 0; idx < numIter; idx++) {
                const auto numR = numRbase << idx;
                this->build_RRsets(numR, true, "select");
                this->build_RRsets(numR, true, "val");
                this->update_RRsets_delta();
                sampling_time += timer_ScaLIM.get_operation_time();

                VecLargeNum selected_idx_vec;
                double cov_num;
                // if args.pruned is true
                cov_num = this->edge_selection_greedy();
                
                selection_time += timer_ScaLIM.get_operation_time();

                double cov_val = this->comp_cov_pmc(this->_vec_selected_edges, "val");

                min_u_bound = std::min(min_u_bound, (double)cov_num);
                if (idx == 0) last_bound = cov_num;
                log_info("cov_prob", cov_num / this->_RRsets.size());
                const auto lowerSelect = pow2(sqrt(cov_val + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
                const auto upperOPT = pow2(sqrt(last_bound / approx + a2 / 2.0) + sqrt(a2 / 2.0));
                const auto approxOPIMC = lowerSelect / upperOPT;
                last_bound = cov_num;
                std::cout << " -->ScaLIM (" << idx + 1 << "/" << numIter << ") approx. (max-cover): " << approxOPIMC << ", #RR sets: " << this->_RRsets.size() << '\n';
                // double converge_err = 0.9;
                // if (approxOPIMC >= approx - this->args.epsilon && ((approxOPIMC - last_approxC < this->args.epsilon/2) || approxOPIMC > converge_err)){
                if (approxOPIMC >= approx - this->args.epsilon) {
                    // Setting Results
                    res.set_k_edges(this->args.k_edges);
                    res.set_approximation(approxOPIMC);
                    res.set_influence_original(this->comp_cov_pmc(this->_vec_selected_edges, "select") * this->n / this->_cur_RRsets_num);
                    res.set_influence(cov_val * this->n / this->_cur_RRsets_num);
                    res.set_sampling_time(sampling_time);
                    res.set_selection_time(selection_time);
                    res.set_RR_sets_size(this->_cur_RRsets_num*2);
                    res.save_to_file(log_filepath);
                    break;
                }
                last_approxC = approxOPIMC;
            }
            // this->count_RRsets_size();
            string edges_save_path = this->_cur_working_folder+"selected_edges_ScaLIM.txt";
            write_UVWEdges(this->_vec_selected_edges, edges_save_path);
            return;
        }

        void OPIM_ScaLIM_minus() {
            ASSERT(this->_seed_set_to_augment.size() > 0);
            ASSERT(this->_vec_cand_edges.size() > 0);
            
            // Initialize parameters
            const double delta = this->args.delta;
            const double e = exp(1);
            const double approx = 1 - 1.0 / e;
            const double alpha = sqrt(log(6.0 / delta));
            const double beta = sqrt((1 - 1 / e) * (logcnk(this->n*this->args.k_seed, this->args.k_edges) + log(6.0 / delta)));
            const auto numRbase = size_t(2.0 * pow2((1 - 1 / e) * alpha + beta));
            const auto maxNumR = size_t(2.0 * this->n * pow2((1 - 1 / e) * alpha + beta) / this->args.k_edges / pow2(this->args.epsilon)) + 1;
            const auto numIter = (size_t)log2(maxNumR / numRbase) + 1;
            const double a1 = log(numIter * 3.0 / delta);
            const double a2 = log(numIter * 3.0 / delta);
            double min_u_bound = __DBL_MAX__;
            double last_bound = __DBL_MAX__;
            double sampling_time= 0.0, selection_time = 0.0;
            ResultInfo res;
            string log_filepath = this->_cur_working_folder + "log_ScaLIM_minus.txt";
            Timer timer_ScaLIM("OPIM_ScaLIM_minus");
            
            this->organize_edges("vec");
            timer_ScaLIM.log_operation_time("organizing");
            double last_approxC = 0.0;

            for (auto idx = 0; idx < numIter; idx++) {
                const auto numR = numRbase << idx;
                this->build_RRsets(numR, true, "select");
                this->build_RRsets(numR, true, "val");
                this->update_RRsets_delta();
                sampling_time += timer_ScaLIM.get_operation_time();

                VecLargeNum selected_idx_vec;
                // selection
                log_info("---Start selection---");
                double cov_num = this->edge_selection_on_edge_RRsets();
                selection_time += timer_ScaLIM.get_operation_time();
                log_info("---Start validation---");
                double cov_val = this->comp_cov_pmc(this->_vec_selected_edges, "val");

                min_u_bound = std::min(min_u_bound, (double)cov_num);
                if (idx == 0) last_bound = cov_num;
                log_info("cov_prob", cov_num / this->_RRsets.size());
                const auto lowerSelect = pow2(sqrt(cov_val + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
                const auto upperOPT = pow2(sqrt(last_bound / approx + a2 / 2.0) + sqrt(a2 / 2.0));
                const auto approxOPIMC = lowerSelect / upperOPT;
                last_bound = cov_num;
                std::cout << " -->ScaLIM (" << idx + 1 << "/" << numIter << ") approx. (max-cover): " << approxOPIMC << ", #RR sets: " << this->_RRsets.size() << '\n';
                
                // if ((approxOPIMC >= approx - this->args.epsilon && approxOPIMC - last_approxC < this->args.epsilon/2) || approxOPIMC > 0.9) {
                if (approxOPIMC >= approx - this->args.epsilon) {
                    // Setting Results
                    res.set_k_edges(this->args.k_edges);
                    res.set_approximation(approxOPIMC);
                    res.set_influence_original(this->comp_cov_pmc(this->_vec_selected_edges, "select") * this->n / this->_cur_RRsets_num);
                    res.set_influence(cov_val * this->n / this->_cur_RRsets_num);
                    res.set_sampling_time(sampling_time);
                    res.set_selection_time(selection_time);
                    res.set_RR_sets_size(this->_cur_RRsets_num*2);
                    res.save_to_file(log_filepath);
                    break;
                }
                last_approxC = approxOPIMC;
            }
            // this->count_RRsets_size();
            string edges_save_path = this->_cur_working_folder+"selected_edges_ScaLIM_minus.txt";
            write_UVWEdges(this->_vec_selected_edges, edges_save_path);
            return;
        }

        void OPIM_Greedy2() {
            ASSERT(this->_seed_set_to_augment.size() > 0);
            ASSERT(this->_vec_cand_edges.size() > 0);
            this->organize_edges("vec");
            // Initialize parameters
            const double delta = this->args.delta;
            const double e = exp(1);
            const double approx = 1 - 1.0 / e;
            const double alpha = sqrt(log(6.0 / delta));
            const double beta = sqrt((1 - 1 / e) * (logcnk(this->n*this->args.k_seed, this->args.k_edges) + log(6.0 / delta)));
            const auto numRbase = size_t(2.0 * pow2((1 - 1 / e) * alpha + beta));
            const auto maxNumR = size_t(2.0 * this->n * pow2((1 - 1 / e) * alpha + beta) / this->args.k_edges / pow2(this->args.epsilon)) + 1;
            const auto numIter = (size_t)log2(maxNumR / numRbase) + 1;
            const double a1 = log(numIter * 3.0 / delta);
            const double a2 = log(numIter * 3.0 / delta);
            double min_u_bound = __DBL_MAX__;
            double last_bound = __DBL_MAX__;

            double sampling_time= 0.0, sampling_2_time = 0.0, selection_time = 0.0;
            ResultInfo res;
            string log_filepath = this->_cur_working_folder + "log_Greedy2.txt";
            
            Timer timer_Greedy2("OPIM_Greedy2");
            
            VecVecLargeNum edge_to_RRsets(this->_vec_cand_edges.size());
            VecVecLargeNum RRsets_for_edges;
            VecVecLargeNum edge_to_RRsets_val(this->_vec_cand_edges.size());
            VecVecLargeNum RRsets_for_edges_val;
            this->organize_edges("vec");
            timer_Greedy2.log_operation_time("organizing");

            double last_approxC = 0.0;
            
            for (auto idx = 0; idx < numIter; idx++) {
                const auto numR = numRbase << idx;
                this->build_RRsets(numR, true, "select");
                this->build_RRsets(numR, true, "val");
                sampling_time += timer_Greedy2.get_operation_time();
                this->sample_edges(edge_to_RRsets, RRsets_for_edges, "select");
                this->sample_edges(edge_to_RRsets_val, RRsets_for_edges_val, "val");
                sampling_2_time += timer_Greedy2.get_operation_time();

                VecLargeNum selected_idx_vec;
                double cov_num = this->edge_selection_by_max_cover(edge_to_RRsets, RRsets_for_edges, selected_idx_vec);
                selection_time += timer_Greedy2.get_operation_time();

                double cov_val = comp_cov(edge_to_RRsets_val, RRsets_for_edges_val, selected_idx_vec);

                min_u_bound = std::min(min_u_bound, (double)cov_num);
                if (idx == 0) last_bound = cov_num;
                log_info("cov_num", cov_num / (double)RRsets_for_edges.size());
                const auto lowerSelect = pow2(sqrt(cov_val + a1 * 2.0 / 9.0) - sqrt(a1 / 2.0)) - a1 / 18.0;
                const auto upperOPT = pow2(sqrt(last_bound / approx + a2 / 2.0) + sqrt(a2 / 2.0));
                const auto approxOPIMC = lowerSelect / upperOPT;
                last_bound = cov_num;
                std::cout << " -->Greedy2 (" << idx + 1 << "/" << numIter << ") approx. (max-cover): " << approxOPIMC << ", #RR sets: " << this->_RRsets.size() << '\n';
                // if ((approxOPIMC >= approx - this->args.epsilon && approxOPIMC - last_approxC < this->args.epsilon/2) || approxOPIMC > 0.9) {
                if (approxOPIMC >= approx - this->args.epsilon) {
                    // Setting Results
                    res.set_k_edges(this->args.k_edges);
                    res.set_approximation(approxOPIMC);
                    res.set_influence_original(comp_cov(edge_to_RRsets, RRsets_for_edges, selected_idx_vec) * this->n / this->_cur_RRsets_num);
                    res.set_influence(cov_val * this->n / this->_cur_RRsets_num);
                    res.set_sampling_time(sampling_time);
                    res.set_sampling_2_time(sampling_2_time);
                    res.set_selection_time(selection_time);
                    res.set_RR_sets_size(this->_cur_RRsets_num*2);
                    res.save_to_file(log_filepath);
                    break;
                }
                last_approxC = approxOPIMC;
            }
            string edges_save_path = this->_cur_working_folder+"selected_edges_Greedy2.txt";
            write_UVWEdges(this->_vec_selected_edges, edges_save_path);
            return;
        }
        // Entry to Greedy2
        void Greedy2() {
            if (this->args.num_samples == 0) {
                this->OPIM_Greedy2();
            } else
                this->fixed_Greedy2();
        }
        
        // Entry to ScaLIM
        void ScaLIM() {
            if (this->args.num_samples == 0) {
                this->OPIM_ScaLIM();
            } else
                this->fixed_ScaLIM();
        }
        //Entry to ScaLIM_minus
        void ScaLIM_minus() {
            if (this->args.num_samples == 0) {
                this->OPIM_ScaLIM_minus();
            } else
                this->fixed_ScaLIM_minus();
        }
        
        void create_param_dir() {
            if (!check_file_exist(this->folder + "params/")) {
                log_info("--- Creating params folder ---");
                make_dir(this->folder + "params/");
            }

            string cur_folder_name = this->folder + "params/" + to_string(this->args.k_edges) + "_" + to_string(this->args.num_cand_edges) + "_" + to_string(this->args.epsilon) + "_" + to_string(this->args.delta) + "_" + to_string(this->args.rand_seed) + "_" + this->_prob_mode;
            if (this->_seed_mode != "IM")
                cur_folder_name = this->folder + "params/" + to_string(this->args.k_edges) + "_" + to_string(this->args.num_cand_edges) + "_" + to_string(this->args.epsilon) + "_" + to_string(this->args.delta) + "_" + to_string(this->args.rand_seed) + "_" + this->_seed_mode + "seed_" + this->_prob_mode;
            if (this->args.beta > 1) {
                cur_folder_name = cur_folder_name + "_beta_" + to_string(this->args.beta);
            }
            if (this->args.num_samples != 0) {
                cur_folder_name = cur_folder_name + "_numsamples_" + to_string(this->args.num_samples);
            }
            bool flag = make_dir(cur_folder_name);
            if (flag)
            {
                log_info("--- Current experiments' data will be stored in " + cur_folder_name + " ---");
                this->_cur_working_folder = cur_folder_name + "/";
            }
            else log_info("--- Failed to create the dir ---");
            return;
        }
};
#endif