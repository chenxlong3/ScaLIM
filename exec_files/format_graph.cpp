#include <iostream>
#include "Experiments.h"
#include "CommonFuncs.h"

int main(int argc, char *argv[])
{
   string folder_name = "";
   string filename = "";
   for (int i = 0; i < argc; i++) {
       if (argv[i] == string("-dataset")) folder_name = string(argv[i + 1]) + "/";
       if (argv[i] == string("-filename")) filename = string(argv[i + 1]);
   }

    if (folder_name == "" || filename == "")
        ExitMessage("argument dataset / filename missing");
   InfGraph g;
   g.format_graph(folder_name, filename);
   // experiment_efficacy(g);
   return 0;
}