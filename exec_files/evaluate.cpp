#include "Experiments.h"
#include "CommonFuncs.h"

int main(int argc, char *argv[])
{
   Argument args = parse_args(argc, argv);
   dsfmt_gv_init_gen_rand(args.rand_seed);
   InfGraph g(args.folder_name, args.graph_file, args.probability_mode);
   g.set_args(args);    // setting arguments
   // InfGraph g = parseArg(argc, argv);
   string edges_mode = "";
   string eval_mode = "MC";
   uint32_t log_step = -1;
   size_t num_mc = 0;
   for (int i = 0; i < argc; i++) {
      if (argv[i] == string("-edges_mode")) edges_mode = string(argv[i + 1]);
      if (argv[i] == string("-log_step")) log_step = atoi(argv[i + 1]);
      if (argv[i] == string("-num_mc")) num_mc = atoi(argv[i + 1]);
   }
   if (edges_mode == "") {
      log_info("Please select edges mode");
      return 0;
   }
   if (log_step == -1) {
      ExitMessage("Please select log_step");
   }
   
   experiment_evaluate_from_file(g, edges_mode, log_step, num_mc);
   return 0;
}