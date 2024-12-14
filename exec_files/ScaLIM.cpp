#include <iostream>
#include "Experiments.h"
#include "CommonFuncs.h"

int main(int argc, char *argv[])
{
   // if needing to 
   string func = "";
   for (int i = 0; i < argc; i++) {
      if (argv[i] == string("-func")) func = string(argv[i + 1]);
   }
   if (func == "format") {
       format_graph(argc, argv);
       return 1;
   }

   Argument args = parse_args(argc, argv);
   dsfmt_gv_init_gen_rand(args.rand_seed);
   InfGraph g(args.folder_name, args.graph_file, args.probability_mode);
   g.set_args(args);    // setting arguments
   
   run_method(g);

   return 1;
}