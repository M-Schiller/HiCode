int get_partition_from_file_list_pajek(std::string s, std::deque<std::deque<int>> & ten, std::deque<int> & oldlabels)
{
  ten.clear();

  std::ifstream inb(s);
  std::string sst;
  getline(inb, sst);

  int a = 1;
  int bb;

  std::map<int, int> id_value;
  std::map <int, int> value_com;
  while (inb >> bb) {
    value_com.insert(std::make_pair(bb, value_com.size()));
    id_value.insert(std::make_pair(a, bb));
    a++;
  }

  int com = 0;
  for (auto it = value_com.begin(); it != value_com.end(); ++it)
    it->second = com++;

  std::deque <int> first;
  for (int i = 0; i < value_com.size(); i++)
    ten.push_back(first);

  for (auto itm = id_value.begin(); itm != id_value.end(); ++itm)
    ten[value_com[itm->second]].push_back(oldlabels[itm->first - 1]);

  return 0;
}

int pajek_format(std::string filename, bool directed)
{
  // filename is the edge list, directed=true if directed

  std::map<int, int> newlabels;
  std::deque<int> old_labels;

  {
    int label = 0;

    std::ifstream inb(filename);

    std::string ins;

    while (getline(inb, ins)) if (!ins.empty() && ins[0] != '#') {
      std::deque<double> ds;
      cast_string_to_doubles(ins, ds);

      if (ds.size() < 2) {
        std::cout << "From file " << filename << ": string not readable:" << ins << std::endl;
      }

      else {
        int innum1 = cast_int(ds[0]);
        int innum2 = cast_int(ds[1]);

        auto itf = newlabels.find(innum1);
        if (itf == newlabels.end()) {
          newlabels.insert(std::make_pair(innum1, label++));
          old_labels.push_back(innum1);
        }

        itf = newlabels.find(innum2);
        if (itf == newlabels.end()) {
          newlabels.insert(std::make_pair(innum2, label++));
          old_labels.push_back(innum2);
        }
      }
    }
  }

  int dim = newlabels.size();
  std::ofstream out("net.paj");
  out << "*Vertices " << dim << std::endl;

  for (int i = 0; i < dim; i++) {
    out << i + 1 << " \"" << old_labels[i] << "\"" << std::endl;
  }

  if (!directed)
    out << "*Edges " << std::endl;
  else
    out << "*Arcs " << std::endl;

  {
    std::ifstream inb(filename);
    std::string ins;
    while (getline(inb, ins)) if (!ins.empty() && ins[0] != '#') {
      std::deque<std::string> ds;
      separate_strings(ins, ds);

      int innum1 = cast_int(cast_string_to_double(ds[0]));
      int innum2 = cast_int(cast_string_to_double(ds[1]));

      double w = 1;
      if (ds.size() > 2) {
        if (cast_string_to_double(ds[2], w) == false)
          w = 1;
      }

      int new1 = newlabels[innum1];
      int new2 = newlabels[innum2];

      //if(directed || new1<=new2)
      out << new1 + 1 << " " << new2 + 1 << " " << w << std::endl;
    }
  }

  std::ofstream lout("labels.dat");
  prints(old_labels, lout);

  //systemCall("unixdos network2.dat network.net");
  //systemCall("rm network2.dat");

  return 0;
}
