#pragma once

#ifndef SET_PARAMETERS_H
#define SET_PARAMETERS_H

// to insert a new parameter there are four steps:
// 1- define it in the class
// 2- initialize it
// 3- set the flag
// 4- set the parameter

inline void general_program_statement(char * b)
{
  std::cout << "USAGE: " << b << " -f network.dat -uw(-w)" << std::endl << std::endl;
  std::cout << "-uw must be used if you want to use the unweighted null model; -w otherwise." << std::endl;
  std::cout << "network.dat is the list of edges. Please look at ReadMe.pdf for more details." << std::endl;

  std::cout << "\n\n\n";
  std::cout << "***************************************************************************************************************************************************" << std::endl;
  std::cout << "OPTIONS" << std::endl;
  std::cout << "\n  [-r R]:\t\t\tsets the number of runs for the first hierarchical level, bigger this value, more accurate the output (of course, it takes more). Default value is 10." << std::endl;
  std::cout << "\n  [-hr R]:\t\t\tsets the number of runs  for higher hierarchical levels. Default value is 50 (the method should be faster since the aggregated network is usually much smaller)." << std::endl;
  std::cout << "\n  [-seed m]:\t\t\tsets m equal to the seed for the random number generator. (instead of reading from time_seed.dat)" << std::endl;
  std::cout << "\n  [-hint filename]:\t\ttakes a partition from filename. The file is expected to have the nodes belonging to the same cluster on the same line." << std::endl;
  std::cout << "\n  [-load filename]:\t\ttakes modules from a tp file you already got in a previous run." << std::endl;
  std::cout << "\n  [-t T]:\t\t\tsets the threshold equal to T, default value is 0.1" << std::endl;
  std::cout << "\n  [-singlet]:\t\t\t finds singletons. If you use this flag, the program generally finds a number of nodes which are not assigned to any module.\n\t\t\t\t";
  std::cout << "the program will assign each node with at least one not homeless neighbor. This only applies to the lowest hierarchical level." << std::endl;
  std::cout << "\n  [-cp P]:\t\t\tsets a kind of resolution parameter equal to P. This parameter is used to decide if it is better to take some modules or their union.\n\t\t\t\tDefault value is 0.5. ";
  std::cout << "Bigger value leads to bigger clusters. P must be in the interval (0, 1)." << std::endl;
  std::cout << "\n  [-fast]:\t\t\tis equivalent to \"-r 1 -hr 1\" (the fastest possible execution)." << std::endl;
  std::cout << "\n  [-infomap runs]:\t\tcalls infomap and uses its output as a starting point. runs is the number of times you want to call infomap." << std::endl;
  std::cout << "\n  [-copra runs]:\t\tsame as above using copra." << std::endl;
  std::cout << "\n  [-louvain runs]:\t\tsame as above using louvain method." << std::endl;
  std::cout << "\n\nPlease look at ReadMe.pdf for a more detailed explanation." << std::endl;

  std::cout << "\n\n\n";
  std::cout << "***************************************************************************************************************************************************" << std::endl;
  std::cout << "OUTPUT FILES" << std::endl << std::endl;

  std::cout << "The program will create a directory called \"[network.dat]_oslo_files\". If the directory is not empty it will cleared, so be careful if you want to save some previous output files.\n" << std::endl;
  std::cout << "All the files will be written in this directory. " << std::endl;
  std::cout << "The first level partition will be written in a file called \"tp\", the next ";
  std::cout << "hierchical network will be recorded as \"net1\", " << std::endl;
  std::cout << "the second level partition will be called \"tp1\" and so on." << std::endl;
  std::cout << "For convenience, the first level partition will be also written in a file called \"tp\" located in the same folder where the program is." << std::endl;
  std::cout << "***************************************************************************************************************************************************" << std::endl;

  std::cout << std::endl << std::endl;
  std::cout << "PLEASE LOOK AT ReadMe.pdf for more details. Thanks!" << std::endl << std::endl << std::endl;
}

inline void error_statement(char * b)
{
  std::cerr << "\n\n************************************************************" << std::endl;
  std::cerr << "ERROR while reading parameters from command line... Please read program instructions or type: \n" << b << std::endl;
  std::cerr << "************************************************************" << std::endl;
}

class Parameters
{
public:

  Parameters();
  ~Parameters() = default;
  void print();
  bool _set_(int argc, char * argv[]);

  //*******************************************************
  std::string file1;
  std::string file2;
  std::string file_load;

  int seed_random;

  double threshold;

  int clean_up_runs;
  int inflate_runs;
  int inflate_stopper;
  double equivalence_parameter;
  int CUT_Off;

  int maxborder_nodes;
  double maxbg_ordinary;
  int iterative_stopper;
  int minimality_stopper;
  double hierarchy_convergence;

  int Or;
  int hier_gather_runs;

  double coverage_inclusion_module_collection;
  double coverage_percentage_fusion_left;
  double check_inter_p;
  double coverage_percentage_fusion_or_submodules;

  bool print_flag_subgraph;

  bool value;
  bool value_load;
  bool fast;
  bool weighted;
  bool homeless_anyway;

  int max_iteration_convergence;

  std::deque<std::string> to_run;
  std::deque<std::string> to_run_part;

  int infomap_runs;
  int copra_runs;
  int louvain_runs;

private:
  std::map<std::string, int> command_flags;
  bool set_flag_and_number(double & threshold, int & argct, int argc, char * argv[], double min_v, double max_v,
    const std::string& warning);
  bool set_flag_and_number(int &, int & argct, int argc, char * argv[], int min_v, int max_v,
    const std::string& warning);
  bool set_flag_and_number_external_program(const std::string& program_name, int & argct, int & number_to_set, int argc, char * argv[]);
};

inline void Parameters::print()
{
  std::cout << "**************************************" << std::endl;
  std::cout << "Threshold:\t\t\t" << threshold << std::endl;
  std::cout << "Network file:\t\t\t" << file1 << std::endl;

  if (weighted)
    std::cout << "Weighted: yes" << std::endl;
  else
    std::cout << "Weighted: no" << std::endl;

  if (fast)
    std::cout << "-fast option selected" << std::endl;

  if (value)
    std::cout << "Hint from file:\t\t\t" << file2 << std::endl;
  if (value_load)
    std::cout << "tp-file:\t\t\t" << file_load << std::endl;

  std::cout << "First Level Runs:\t\t\t" << Or << std::endl;
  std::cout << "Higher Level Runs:\t\t\t" << hier_gather_runs << std::endl;
  std::cout << "-cp:\t\t\t" << coverage_percentage_fusion_or_submodules << std::endl;

  if (seed_random != -1)
    std::cout << "Random number generator seed:\t\t\t" << seed_random << std::endl;

  if (!homeless_anyway)
    std::cout << "-singlet option selected" << std::endl;

  for (unsigned i = 0; i < to_run.size(); i++)
    std::cout << "String to run: [" << to_run[i] << "]\t\t\t\t\t\tModule file: [" << to_run_part[i] << "]" << std::endl;

  std::cout << "**************************************" << std::endl << std::endl;
}

inline bool Parameters::set_flag_and_number(double & number_to_set, int & argct, int argc, char * argv[], double min_v, double max_v,
  const std::string& warning)
{
  argct++;
  if (argct == argc) {
    std::cout << "you didn't set any number for the " << warning << std::endl;
    error_statement(argv[0]);
    return false;
  }

  std::string tt = argv[argct];
  double ttt;
  if (!cast_string_to_double(tt, ttt))
  {
    std::cout << "you didn't set any number for the " << warning << std::endl;
    error_statement(argv[0]);
    return false;
  }

  number_to_set = ttt;

  if (number_to_set<min_v || number_to_set>max_v)
  {
    std::cout << "the " << warning << " must be between " << min_v << " and " << max_v << std::endl;
    error_statement(argv[0]);
    return false;
  }

  return true;
}

inline bool Parameters::set_flag_and_number(int & number_to_set, int & argct, int argc, char * argv[], int min_v, int max_v,
  const std::string& warning)
{
  argct++;
  if (argct == argc) {
    std::cout << "you didn't set any number for the " << warning << std::endl;
    error_statement(argv[0]);
    return false;
  }

  std::string tt = argv[argct];
  double ttt;
  if (!cast_string_to_double(tt, ttt))
  {
    std::cout << "you didn't set any number for the " << warning << std::endl;
    error_statement(argv[0]);
    return false;
  }

  number_to_set = cast_int(ttt);

  if (number_to_set<min_v || number_to_set>max_v)
  {
    std::cout << "the " << warning << " must be between " << min_v << " and " << max_v << std::endl;
    error_statement(argv[0]);
    return false;
  }
  return true;
}

inline Parameters::Parameters()
{
  //**************************************************************************

  seed_random = -1;

  threshold = 0.1;											// this is the P-value for the significance of the module

  clean_up_runs = 25;										// the number of runs in the clean up procedure
  inflate_runs = 3;											// the number of runs in the clean up of the inflate procedure
  inflate_stopper = 5;										// the number of runs in the inflate procedure
  equivalence_parameter = 0.33;								// this parameters tells when nodes are considered equivalent in the clean up procedure
  CUT_Off = 200;											// this is used in the inflate function

  maxborder_nodes = 100;									// this is to speed up the code in looking for "reasonably good" neighbors
  maxbg_ordinary = 0.1;										// same as above
  iterative_stopper = 10;									// this is to prevent the iterative procedure to last too long. this can happen in case of strong backbones (just an idea, not sure)
  minimality_stopper = 10;									// this is to prevent too many minimality checks
  hierarchy_convergence = 0.05;								// this parameter is used to stop the hierarchical process when not enough modules are found

  Or = 10;													// this is the number of global runs in the gather function		(first level)
  hier_gather_runs = 50;									// this is the number of global runs in the gather function		(higher level)

  coverage_inclusion_module_collection = 0.49999;			// this is used to see if two modules are higly similar in processing the clusters (big_module)
  coverage_percentage_fusion_left = 0.8;					// this is used to see when fusing clusters how much is left
  check_inter_p = 0.05;										// this parameter is a check parameter for the fusion of quite similar clusters
  coverage_percentage_fusion_or_submodules = 0.5;			// this is the resolution parameter to decide between split clusters or unions, if you increase this value the program tends to find bigger clusters

  print_flag_subgraph = true;										// this flag is used to print things when necessary

  /* these are some flags to read input files */
  value = false;
  value_load = false;
  fast = false;
  weighted = false;
  homeless_anyway = true;

  //********************* collect_groups
  max_iteration_convergence = 10;							// parameter for the convergence of the collect_groups function

  infomap_runs = 0;
  copra_runs = 0;
  louvain_runs = 0;

  command_flags.emplace("-w", 1);
  command_flags.emplace("-uw", 2);
  command_flags.emplace("-singlet", 3);
  command_flags.emplace("-f", 4);
  command_flags.emplace("-hint", 5);
  command_flags.emplace("-load", 6);
  command_flags.emplace("-t", 7);
  command_flags.emplace("-r", 8);
  command_flags.emplace("-hr", 9);
  command_flags.emplace("-seed", 10);
  command_flags.emplace("-cp", 11);
  command_flags.emplace("-fast", 12);
  command_flags.emplace("-infomap", 13);
  command_flags.emplace("-copra", 14);
  command_flags.emplace("-louvain", 15);
}

inline bool Parameters::set_flag_and_number_external_program(const std::string& program_name, int & argct, int & number_to_set, int argc, char * argv[])
{
  argct++;
  if (argct == argc)
  {
    std::cout << "you didn't set the number of " << program_name << std::endl;

    error_statement(argv[0]);
    return false;
  }

  std::string tt = argv[argct];
  double ttt;
  if (!cast_string_to_double(tt, ttt))
  {
    std::cout << "you didn't set the number of " << program_name << std::endl;

    error_statement(argv[0]);
    return false;
  }

  number_to_set = cast_int(ttt);

  if (number_to_set < 0)
  {
    std::cout << " the number of " << program_name << " must be positive" << std::endl;

    error_statement(argv[0]);
    return false;
  }

  return true;
}

inline bool Parameters::_set_(int argc, char * argv[])
{
  int argct = 0;

  if (argc <= 1)			/* if no arguments, return error_statement about program print_usage.*/
  {
    error_statement(argv[0]);
    return false;
  }

  bool f_set = false;
  bool set_weighted = false;

  while (++argct < argc)			// input file name
  {
    std::cout << "setting " << argv[argct] << std::endl;
    std::string temp = argv[argct];
    const auto itf = command_flags.find(temp);

    if (itf == command_flags.end())
    {
      error_statement(argv[0]);
      return false;
    }

    int vp = itf->second;

    switch (vp)
    {
    case 1:
      weighted = true;
      set_weighted = true;
      break;
    case 2:
      weighted = false;
      set_weighted = true;
      break;
    case 3:
      homeless_anyway = false;
      break;
    case 4:
      argct++;
      if (argct == argc)
      {
        error_statement(argv[0]);
        return false;
      }
      file1 = argv[argct];
      f_set = true;
      break;
    case 5:
      argct++;
      if (argct == argc)
      {
        error_statement(argv[0]);
        return false;
      }
      file2 = argv[argct];
      value = true;
      break;
    case 6:
      argct++;
      if (argct == argc)
      {
        error_statement(argv[0]);
        return false;
      }
      file_load = argv[argct];
      value_load = true;
      break;
    case 7:
      if (!set_flag_and_number(threshold, argct, argc, argv, 0., 1., "threshold"))
      {
        return false;
      }
      break;
    case 8:
      if (!set_flag_and_number(Or, argct, argc, argv, 0, R2_IM2, "runs"))
      {
        return false;
      }
      break;
    case 9:
      if (!set_flag_and_number(hier_gather_runs, argct, argc, argv, 0, R2_IM2, "higher-level runs"))
      {
        return false;
      }
      break;
    case 10:
      if (!set_flag_and_number(seed_random, argct, argc, argv, 1, R2_IM2, "seed of the random number generator"))
      {
        return false;
      }
      break;
    case 11:
      if (!set_flag_and_number(coverage_percentage_fusion_or_submodules, argct, argc, argv, 0., 1., "resolution parameter"))
      {
        return false;
      }
      break;
    case 12:
      fast = true;
      break;
    case 13:
      if (!set_flag_and_number_external_program("runs for infomap", argct, infomap_runs, argc, argv))
      {
        return false;
      }
      break;
    case 14:
      if (!set_flag_and_number_external_program("runs for copra", argct, copra_runs, argc, argv))
      {
        return false;
      }
      break;
    case 15:
      if (!set_flag_and_number_external_program("runs for louvain method", argct, louvain_runs, argc, argv))
      {
        return false;
      }
      break;
    default:
      error_statement(argv[0]);
      return false;
    }
  }
  /*******************************************************************/

  if (!f_set)
  {
    std::cerr << "\n\n************************************************************" << std::endl;
    std::cout << "ERROR: you didn't set the file with the network.  Please read program instructions or type: \n" << argv[0] << std::endl;
    std::cerr << "************************************************************" << std::endl;

    return false;
  }

  if (!set_weighted)
  {
    std::cerr << "\n\n************************************************************" << std::endl;
    std::cout << "ERROR: you didn't set the option -w (weighted network) or -uw (unweighted network).  Please read program instructions or type: \n" << argv[0] << std::endl;
    std::cerr << "************************************************************" << std::endl;

    return false;
  }

  if (seed_random == -1)
  {
    srand_file();
  }
  else
  {
    srand5(seed_random);
  }

  if (fast)
  {
    Or = 1;
    hier_gather_runs = 1;
  }

  for (int i = 0; i < infomap_runs; i++)
  {
    char number_r[1000];

    //cout<<"************** "<<string(argv[0])<< std::endl;
    std::string pros(argv[0]);

    if (pros == "./oslom_undir")
      sprintf(number_r, "./infomap_undir_script NETx %d %d", irand(10000000), 1);
    else
      sprintf(number_r, "./infomap_dir_script NETx %d %d", irand(10000000), 1);

    std::string sr(number_r);

    //cout<<"here "<< std::endl;
    to_run.push_back(sr);
    to_run_part.emplace_back("infomap.part");
  }

  for (int i = 0; i < copra_runs; i++)
  {
    char number_r[1000];
    sprintf(number_r, "java -cp copra.jar COPRA NETx -v 5 -w");
    std::string sr(number_r);

    to_run.push_back(sr);
    to_run_part.emplace_back("clusters-NETx");
  }

  for (int i = 0; i < louvain_runs; i++)
  {
    char number_r[1000];
    sprintf(number_r, "./louvain_script -f NETx");
    std::string sr(number_r);

    to_run.push_back(sr);
    to_run_part.emplace_back("louvain.part");
  }
  return true;
}

#endif // SET_PARAMETERS_H
