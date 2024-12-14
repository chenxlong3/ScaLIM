#ifndef EXPERIMENTS_H
#define EXPERIMENTS_H
#include "InfGraph.h"
#include "CommonFuncs.h"
#include "Timer.h"

Argument parse_args(int argn, char *argv[]) {
    Argument args;
    string folder_name = "";
    string probability_mode = "NONE";
    double epsilon = 0;
    string model = "";
    uint32_t k_seed=0, k_edges = 0, beta = 1;
    uint32_t rand_seed = 2023;
    double delta = 0.0;
    CascadeModel casc_model;
    bool fast_truncated = false;
    uint32_t num_cand_edges;
    string seed_mode = "IM";
    // For test
    model = "IC";

    for (int i = 0; i < argn; i++) {
        if (argv[i] == string("-dataset")) args.folder_name = string(argv[i + 1]) + "/";
        if (argv[i] == string("-graph_file")) args.graph_file = string(argv[i + 1]);
        if (argv[i] == string("-epsilon")) args.epsilon = atof(argv[i + 1]);
        if (argv[i] == string("-k_seed")) args.k_seed = atoi(argv[i + 1]);
        if (argv[i] == string("-k_edges")) args.k_edges = atoi(argv[i + 1]);
        if (argv[i] == string("-num_cand_edges")) args.num_cand_edges = atoi(argv[i + 1]);
        if (argv[i] == string("-delta")) args.delta = atof(argv[i+1]);
        if (argv[i] == string("-rand_seed")) args.rand_seed = atoi(argv[i + 1]);
        if (argv[i] == string("-beta")) args.beta = atoi(argv[i + 1]);
        if (argv[i] == string("-probability")) args.probability_mode = string(argv[i + 1]);
        if (argv[i] == string("-seed_mode")) args.seed_mode = string(argv[i + 1]);
        if (argv[i] == string("-method")) args.method = string(argv[i + 1]);
        if (argv[i] == string("-num_samples")) args.num_samples = size_t(atoi(argv[i + 1]));

        if (argv[i] == string("-fast")) {
            string str_fast = string(argv[i+1]);
            if (str_fast == string("True")) args.fast_truncated = true;
            else if (str_fast == string("False")) args.fast_truncated = false;
            else {
                ExitMessage("Illegal input of fast: should be True or False");
            }
        }
        if (argv[i] == string("-pruned")) {
            string str_fast = string(argv[i+1]);
            if (str_fast == string("True")) args.pruned = true;
            else if (str_fast == string("False")) args.pruned = false;
            else {
                ExitMessage("Illegal input of fast: should be True or False");
            }
        }
        if (argv[i] == string("-model")) {
            if (argv[i + 1] == string("LT")) {
                args.str_model = argv[i + 1];
                args.model = LT;
            } else if (argv[i + 1] == string("IC")) {
                args.str_model = argv[i + 1];
                args.model = IC;
            } else
                ExitMessage("model should be IC or LT");
        }
    }
    args.check_arguments_eligible();
    return args;
}

void format_graph(int argn, char *argv[]) {
    string folder_name = "";
    string filename = "";
    for (int i = 0; i < argn; i++) {
        if (argv[i] == string("-dataset")) folder_name = string(argv[i + 1]) + "/";
        if (argv[i] == string("-filename")) filename = string(argv[i + 1]);
    }

    if (folder_name == "" || filename == "")
        ExitMessage("argument dataset / filename missing");
    InfGraph g;
    g.format_graph(folder_name, filename);
    return;
}


void evaluate_G_inf_IMA(InfGraph& g, uint32_t num_mc, string edges_mode="IMA", uint32_t log_step = 5) {
    log_info("--- Start evaluating ---");
    VecDouble inf_spread;
    uint32_t upper = g.args.k_edges;
    VecuInt32 log_points;
    for (uint32_t i = 0; i <= upper; i += log_step) {
        log_points.push_back(i);
    }
    string selected_edges_filename = g._cur_working_folder + string("selected_edges_") + edges_mode + string(".txt");
    if (!check_file_exist(selected_edges_filename)) {
        log_info("VITAL: The edges_mode is not legal!!!");
        return;
    }

    std::vector<UVWEdge> vec_selected_edges = read_UVWEdges(selected_edges_filename);
    log_info("--- Using edges selected by " + edges_mode + " ---");
    
    ASSERT(g._seed_set_to_augment.size() > 0);
    ASSERT(log_points.back() <= vec_selected_edges.size());
    inf_spread.push_back(inf_evaluate(g._G, g._seed_set_to_augment, g.args.model, num_mc));
    if (vec_selected_edges.size() == 0) {
        log_info("--- No edge is selected ---");
        for (int i = 0; i < log_points.size() - 1; i++) {
            double spread = inf_spread[0];
            inf_spread.push_back(spread);
        }
        return;
    }

    for (int i = 0; i < log_points.size() - 1; i++) {
        log_info("--- Start adding edges --- k=" + to_string(log_points[i]));
        g.add_edges(vec_selected_edges, log_points[i], log_points[i + 1]);
        log_info("--- Start simulation ---");
        double spread = inf_evaluate(g._G, g._seed_set_to_augment, g._casc_model, num_mc);
        log_info(spread);
        inf_spread.push_back(spread);
    }

    string filename = g._cur_working_folder + string("k_inf_spread_") + edges_mode + string(".csv");
    ofstream outfile(filename);
    log_info("--- Saving Experiment Results ---");
    outfile << 'k' << ',' << "expected spread" << '\n';
    for (uint32_t i = 0; i < log_points.size(); i++)
    {
        outfile << log_points[i] << ',' << inf_spread[i] << '\n';
    }
    outfile.close();
    return;
}

void run_method(InfGraph& g) {
    g.create_param_dir();

    log_info("--- Setting seed set ---");
    if (check_file_exist(g._seed_filename)) {
        g.set_seed(g._seed_mode, true);
    } else {
        // If you choose influential seed set, please run an IM algorithm first to generate the seed set
        g.set_seed(g._seed_mode, false);
    }
    if (check_file_exist(g._cand_edges_filename) || check_file_exist(g._cand_edges_filename + ".vec")) {
        log_info("--- Reading the candidate edges ---");
        g.read_cand_edges();
    } else {
        log_info("--- Generating the candidate edges ---");
        if (g.args.num_cand_edges == 0) {
         g.generate_candidate_edges();
        } else
         g.generate_candidate_edges(g.args.num_cand_edges);
    }
    string method_name = g.args.method;
    string timer_name = method_name + "_" + g.folder;
    Timer timer = Timer(timer_name.c_str());
    // string param_folder = g._cur_working_folder;

    if (method_name == "OPIM")
    {
        // std::cout << g._inv_G[2][0].first << " " << g._inv_G[2][0].second << std::endl;
        g.clean_RRsets_InfGraph();
        g.Greedy2();
        g.clean_RRsets_InfGraph();
        g.ScaLIM();
    }
    // method is ScaLIM_minus
    if (method_name == "ScaLIM_minus")
    {
        g.ScaLIM_minus();
    }
    if (method_name == "ScaLIM") {
        g.ScaLIM();
    }
    if (method_name == "Greedy2")
    {
        g.Greedy2();
    }
}

void experiment_evaluate_from_file(InfGraph& g, string edges_mode="IMA", uint32_t log_step=5, size_t num_mc=0) {
    g.create_param_dir();
    g.init_from_files(false);
    evaluate_G_inf_IMA(g, 10000, edges_mode, log_step);
    return;
}

#endif