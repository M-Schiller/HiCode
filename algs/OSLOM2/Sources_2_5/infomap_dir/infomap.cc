#include "infomap.h"

unsigned stou(char *s) {
  return strtoul(s, (char **)NULL, 10);
}

void printTree(std::string s, std::multimap<double, treeNode, std::greater<>>::iterator it_tM, std::ofstream *outfile, bool flip);
void repeated_partition(MTRand *R, Node ***node, GreedyBase *greedy, bool silent, int Ntrials);
void partition(MTRand *R, Node ***node, GreedyBase *greedy, bool silent);

// Call: trade <seed> <Ntries>
int main(int argc, char *argv[])
{
  if (argc != 4)
  {
    std::cout << "Call: ./infomap <seed> <network.net> <# attempts>" << std::endl;
    exit(-1);
  }

  int Ntrials = atoi(argv[3]);  // Set number of partition attempts
  std::string infile = std::string(argv[2]);
  std::string networkName(infile.begin(), infile.begin() + infile.find('.'));
  std::string line;
  std::string buf;

  MTRand *R = new MTRand(stou(argv[1]));
  std::ostringstream oss;
  std::ofstream outfile;

  /* Read network in Pajek format with nodes ordered 1, 2, 3, ..., N,            */
  /* each directed link occurring only once, and link weights > 0.               */
  /* For more information, see http://vlado.fmf.uni-lj.si/pub/networks/pajek/.   */
  /* Node weights are optional and sets the relative proportion to which         */
  /* each node receives teleporting random walkers. Default value is 1.          */
  /* Example network with three nodes and four directed and weighted links:      */
  /* *Vertices 3                                                                 */
  /* 1 "Name of first node" 1.0                                                  */
  /* 2 "Name of second node" 2.0                                                 */
  /* 3 "Name of third node" 1.0                                                  */
  /* *Arcs 4                                                                     */
  /* 1 2 1.0                                                                     */
  /* 1 3 1.7                                                                     */
  /* 2 3 2.0                                                                     */
  /* 3 2 1.2                                                                     */

  std::cout << "Reading network " << argv[2] << "..." << std::flush;
  std::ifstream net(argv[2]);
  int Nnode = 0;
  std::istringstream ss;
  while (Nnode == 0)
  {
    if (!getline(net, line))
    {
      std::cout << "the network file is not in Pajek format...exiting" << std::endl;
      exit(-1);
    }
    else {
      ss.clear();
      ss.str(line);
      ss >> buf;
      if (buf == "*Vertices" || buf == "*vertices" || buf == "*VERTICES")
      {
        ss >> buf;
        Nnode = std::stoi(buf);
      }
      else {
        std::cout << "the network file is not in Pajek format...exiting" << std::endl;
        exit(-1);
      }
    }
  }

  std::vector<std::string> nodeNames(Nnode);
  std::vector<double> nodeWeights = std::vector<double>(Nnode, 1.0);
  double totNodeWeights = 0.0;

  // Read node names, assuming order 1, 2, 3, ...
  for (int i = 0; i < Nnode; i++)
  {
    getline(net, line);
    int nameStart = line.find_first_of('\"');
    int nameEnd = line.find_last_of('\"');
    if (nameStart < nameEnd)
    {
      nodeNames[i] = std::string(line.begin() + nameStart + 1, line.begin() + nameEnd);
      line = std::string(line.begin() + nameEnd + 1, line.end());
      ss.clear();
      ss.str(line);
    }
    else
    {
      ss.clear();
      ss.str(line);
      ss >> buf;
      ss >> nodeNames[i];
    }

    buf = "1";
    ss >> buf;
    double nodeWeight = std::stof(buf);
    if (nodeWeight <= 0.0)
    {
      nodeWeight = 1.0;
    }
    nodeWeights[i] = nodeWeight;
    totNodeWeights += nodeWeight;
  }

  // Read the number of links in the network
  getline(net, line);
  ss.clear();
  ss.str(line);
  ss >> buf;

  if (buf != "*Edges"
    && buf != "*edges"
    && buf != "*Arcs"
    && buf != "*arcs")
  {
    std::cout << std::endl << "Number of nodes not matching, exiting" << std::endl;
    exit(-1);
  }

  int Nlinks = 0;
  int NdoubleLinks = 0;
  std::map<int, std::map<int, double>> Links;

  // Read links in format "from to weight", for example "1 3 0.7"
  while (getline(net, line))
  {
    ss.clear();
    ss.str(line);
    ss >> buf;
    int linkEnd1 = std::stof(buf);
    ss >> buf;
    int linkEnd2 = std::stoi(buf);
    buf.clear();
    ss >> buf;
    double linkWeight;
    if (buf.empty()) // If no information
    {
      linkWeight = 1.0;
    }
    else
    {
      linkWeight = std::stof(buf);
    }
    linkEnd1--; // Nodes start at 1, but C++ arrays at 0.
    linkEnd2--;

    // Aggregate link weights if they are definied more than once
    auto fromLink_it = Links.find(linkEnd1);
    if (fromLink_it == Links.end())
    { // new link
      std::map<int, double> toLink;
      toLink.emplace(linkEnd2, linkWeight);
      Links.emplace(linkEnd1, toLink);
      Nlinks++;
    }
    else
    {
      auto toLink_it = fromLink_it->second.find(linkEnd2);
      if (toLink_it == fromLink_it->second.end()) // new link
      {
        fromLink_it->second.emplace(linkEnd2, linkWeight);
        Nlinks++;
      }
      else
      {
        toLink_it->second += linkWeight;
        NdoubleLinks++;
      }
    }
  }
  net.close();

  std::cout << "done! (found " << Nnode << " nodes and " << Nlinks << " links";
  if (NdoubleLinks > 0)
    std::cout << ", aggregated " << NdoubleLinks << " link(s) defined more than once";

  /////////// Partition network /////////////////////
  Node **node = new Node*[Nnode];
  for (int i = 0; i < Nnode; i++) {
    node[i] = new Node(i, nodeWeights[i] / totNodeWeights);
  }

  int NselfLinks = 0;
  for (auto fromLink_it = Links.begin(); fromLink_it != Links.end(); ++fromLink_it)
  {
    for (auto toLink_it = fromLink_it->second.begin(); toLink_it != fromLink_it->second.end(); ++toLink_it)
    {
      int from = fromLink_it->first;
      int to = toLink_it->first;
      double weight = toLink_it->second;
      if (weight > 0.0)
      {
        if (from == to)
        {
          NselfLinks++;
        }
        else
        {
          node[from]->outLinks.emplace_back(to, weight);
          node[to]->inLinks.emplace_back(from, weight);
        }
      }
    }
  }

  std::cout << ", ignoring " << NselfLinks << " self link(s)." << std::endl;
  //cout << ", including " <<  NselfLinks << " self link(s))." << endl;

  //Swap vector to free memory
  std::map<int, std::map<int, double>>().swap(Links);

  GreedyBase * greedy = new Greedy(R, Nnode, node, Nnode);
  greedy->initiate();

  std::vector<double> size(Nnode);
  for (int i = 0; i < Nnode; i++)
    size[i] = node[i]->size;

  // vector<int> forceclu(Nnode,0);
  // for(int i=3;i<Nnode;i++)
  // 	forceclu[i] = 1;
  // greedy->determMove(forceclu);
  // greedy->level(&node,true);

  std::cout << "Now partition the network:" << std::endl;
  repeated_partition(R, &node, greedy, false, Ntrials);
  int Nmod = greedy->Nnode;
  std::cout << "Done! Code length " << greedy->codeLength / log(2.0) << " in " << Nmod << " modules." << std::endl;

  // Order links by size
  std::vector<double> exit(Nmod, 0.0);
  std::multimap<double, std::pair<int, int>, std::greater<>> sortedLinks;
  for (int i = 0; i < Nmod; i++)
  {
    int NoutLinks = node[i]->outLinks.size();
    for (int j = 0; j < NoutLinks; j++)
    {
      double linkFlow = node[i]->outLinks[j].second / greedy->beta;
      sortedLinks.emplace(linkFlow, std::make_pair(i + 1, node[i]->outLinks[j].first + 1));
      exit[i] += linkFlow;
    }
  }

  // Order modules by size
  std::multimap<double, treeNode, std::greater<>> treeMap;
  for (int i = 0; i < greedy->Nnode; i++)
  {
    treeNode tmp_tN;
    auto it_tM = treeMap.emplace(node[i]->size, tmp_tN);
    it_tM->second.exit = exit[i];
    for (int& member : node[i]->members)
    {
      it_tM->second.members.emplace(size[member], std::make_pair(member, nodeNames[member]));
    }
  }

  //Print partition in format "module:rank size name"
  oss.str("");
  oss << networkName << ".tree";
  outfile.open(oss.str().c_str());
  outfile << "# Code length " << greedy->codeLength / log(2.0) << " in " << Nmod << " modules." << std::endl;
  int k = 1;
  for (auto it = treeMap.begin(); it != treeMap.end(); ++it)
  {
    std::string s;
    s.append(to_string(k));
    s.append(":");
    printTree(s, it, &outfile, false);
    k++;
  }
  outfile.close();

  // Create cluster vector
  std::vector<int> clusterVec = std::vector<int>(Nnode);
  int clusterNr = 0;
  for (auto mod = treeMap.begin(); mod != treeMap.end(); ++mod)
  {
    for (auto mem = mod->second.members.begin(); mem != mod->second.members.end(); ++mem)
    {
      clusterVec[mem->second.first] = clusterNr;
    }
    clusterNr++;
  }
  // Print partition in Pajek's .clu format
  oss.str("");
  oss << networkName << ".clu";
  outfile.open(oss.str().c_str());
  outfile << "*Vertices " << Nnode << "\x0D\x0A";
  for (int i = 0; i < Nnode; i++)
    outfile << clusterVec[i] + 1 << "\x0D\x0A";
  outfile.close();

  // Print map in Pajek's .net format (links sorted in descending order)
  oss.str("");
  oss << networkName << "_map.net";
  outfile.open(oss.str().c_str());
  outfile << "*Vertices " << Nmod << "\x0D\x0A";
  for (int i = 0; i < Nmod; i++)
    outfile << i + 1 << " \"" << i + 1 << "\"" << "\x0D\x0A";
  outfile << "*Arcs " << sortedLinks.size() << "\x0D\x0A";
  for (auto it = sortedLinks.begin(); it != sortedLinks.end(); ++it)
    outfile << "  " << it->second.first << " " << it->second.second << " " << it->first << "\x0D\x0A";
  outfile.close();

  // Print size of modules in Pajek's .vec format
  oss.str("");
  oss << networkName << "_map.vec";
  outfile.open(oss.str().c_str());
  outfile << "*Vertices " << Nmod << "\x0D\x0A";
  for (int i = 0; i < Nmod; i++)
    outfile << node[i]->size << "\x0D\x0A";
  outfile.close();

  // Print map in .map format for the Map Generator at www.mapequation.org
  oss.str("");
  oss << networkName << ".map";
  outfile.open(oss.str().c_str());
  outfile << "# modules: " << Nmod << std::endl;
  outfile << "# modulelinks: " << sortedLinks.size() << std::endl;
  outfile << "# nodes: " << Nnode << std::endl;
  outfile << "# links: " << Nlinks << std::endl;
  outfile << "# codelength: " << greedy->codeLength / log(2.0) << std::endl;
  outfile << "*Directed" << std::endl;
  outfile << "*Modules " << Nmod << std::endl;
  k = 0;
  for (std::multimap<double, treeNode, std::greater<double>>::iterator it = treeMap.begin(); it != treeMap.end(); ++it)
  {
    outfile << k + 1 << " \"" << it->second.members.begin()->second.second << "\" " << it->first << " " << it->second.exit << std::endl;
    k++;
  }
  outfile << "*Nodes " << Nnode << std::endl;
  k = 1;
  for (std::multimap<double, treeNode, std::greater<double>>::iterator it = treeMap.begin(); it != treeMap.end(); ++it)
  {
    std::string s;
    s.append(to_string(k));
    s.append(":");
    printTree(s, it, &outfile, true);
    k++;
  }
  outfile << "*Links " << sortedLinks.size() << std::endl;
  for (std::multimap<double, std::pair<int, int>, std::greater<double>>::iterator it = sortedLinks.begin(); it != sortedLinks.end(); ++it)
    outfile << it->second.first << " " << it->second.second << " " << 1.0*it->first << std::endl;
  outfile.close();

  //   // print size of vertices (imported as a vector in Pajek)
  //   strcpy(netname,"");
  //   sprintf(netname,"%s-%dmodules.vec",outfile1,greedy->Nmodules);
  //   ofw1 = fopen(netname,"w");
  //   fprintf(ofw1,"*Vertices %d\015\012",greedy->Nmodules);
  //   for(int i=0;i<greedy->Nmodules;i++)
  //       fprintf(ofw1,"%e \015\012",module[i]->prob - module[i]->exit);
  //   fclose(ofw1);

  for (int i = 0; i < greedy->Nnode; i++)
  {
    delete node[i];
  }
  delete[] node;

  delete greedy;
  delete R;
}

void partition(MTRand *R, Node ***node, GreedyBase *greedy, bool silent)
{
  int Nnode = greedy->Nnode;
  Node **cpy_node = new Node*[Nnode];
  for (int i = 0; i < Nnode; i++)
  {
    cpy_node[i] = new Node();
    cpyNode(cpy_node[i], (*node)[i]);
  }

  int iteration = 0;
  double outer_oldCodeLength;
  do {
    outer_oldCodeLength = greedy->codeLength;

    if (iteration > 0
      && iteration % 2 == 0
      && greedy->Nnode > 1)  // Partition the partition
    {
      if (!silent)
        std::cout << "Iteration " << iteration + 1 << ", moving ";

      Node **rpt_node = new Node*[Nnode];
      for (int i = 0; i < Nnode; i++) {
        rpt_node[i] = new Node();
        cpyNode(rpt_node[i], cpy_node[i]);
      }
      std::vector<int> subMoveTo(Nnode);
      std::vector<int> moveTo(Nnode);
      int subModIndex = 0;

      for (int i = 0; i < greedy->Nnode; i++) {
        int sub_Nnode = (*node)[i]->members.size();

        if (sub_Nnode > 1) {
          Node **sub_node = new Node*[sub_Nnode];
          std::set<int> sub_mem;
          for (int j = 0; j < sub_Nnode; j++)
            sub_mem.insert((*node)[i]->members[j]);
          auto it_mem = sub_mem.begin();
          std::vector<int> sub_renumber = std::vector<int>(Nnode);
          std::vector<int> sub_rev_renumber = std::vector<int>(sub_Nnode);
          for (int j = 0; j < sub_Nnode; j++)
          {
            int orig_nr = *it_mem;
            int orig_NoutLinks = cpy_node[orig_nr]->outLinks.size();
            int orig_NinLinks = cpy_node[orig_nr]->inLinks.size();
            sub_renumber[orig_nr] = j;
            sub_rev_renumber[j] = orig_nr;
            sub_node[j] = new Node(j, cpy_node[orig_nr]->teleportWeight);
            sub_node[j]->selfLink = cpy_node[orig_nr]->selfLink; // Take care of self-link
            for (int k = 0; k < orig_NoutLinks; k++)
            {
              int orig_link = cpy_node[orig_nr]->outLinks[k].first;
              int orig_link_newnr = sub_renumber[orig_link];
              double orig_weight = cpy_node[orig_nr]->outLinks[k].second;
              if (orig_link < orig_nr)
              {
                if (sub_mem.find(orig_link) != sub_mem.end())
                {
                  sub_node[j]->outLinks.emplace_back(orig_link_newnr, orig_weight);
                  sub_node[orig_link_newnr]->inLinks.emplace_back(j, orig_weight);
                }
              }
            }
            for (int k = 0; k < orig_NinLinks; k++)
            {
              int orig_link = cpy_node[orig_nr]->inLinks[k].first;
              int orig_link_newnr = sub_renumber[orig_link];
              double orig_weight = cpy_node[orig_nr]->inLinks[k].second;
              if (orig_link < orig_nr)
              {
                if (sub_mem.find(orig_link) != sub_mem.end())
                {
                  sub_node[j]->inLinks.emplace_back(orig_link_newnr, orig_weight);
                  sub_node[orig_link_newnr]->outLinks.emplace_back(j, orig_weight);
                }
              }
            }
            ++it_mem;
          }

          GreedyBase * sub_greedy = new Greedy(R, sub_Nnode, sub_node, sub_Nnode);
          sub_greedy->initiate();
          partition(R, &sub_node, sub_greedy, true);
          for (int j = 0; j < sub_greedy->Nnode; j++)
          {
            int Nmembers = sub_node[j]->members.size();
            for (int k = 0; k < Nmembers; k++) {
              subMoveTo[sub_rev_renumber[sub_node[j]->members[k]]] = subModIndex;
            }
            moveTo[subModIndex] = i;
            subModIndex++;
            delete sub_node[j];
          }

          delete[] sub_node;
          delete sub_greedy;
        }
        else {
          subMoveTo[(*node)[i]->members[0]] = subModIndex;
          moveTo[subModIndex] = i;

          subModIndex++;
        }
      }

      for (int i = 0; i < greedy->Nnode; i++)
        delete (*node)[i];
      delete[] * node;

      greedy->Nnode = Nnode;
      greedy->Nmod = Nnode;
      greedy->Ndanglings = 0;
      greedy->node = rpt_node;
      greedy->calibrate();
      greedy->determineMove(subMoveTo);
      greedy->level(node, false);
      greedy->determineMove(moveTo);
      *node = rpt_node;

      outer_oldCodeLength = greedy->codeLength;

      if (!silent)
        std::cout << greedy->Nnode << " modules, looping ";
    }
    else if (iteration > 0) {
      if (!silent)
        std::cout << "Iteration " << iteration + 1 << ", moving " << Nnode << " nodes, looping ";

      Node **rpt_node = new Node*[Nnode];
      for (int i = 0; i < Nnode; i++) {
        rpt_node[i] = new Node();
        cpyNode(rpt_node[i], cpy_node[i]);
      }

      std::vector<int>moveTo(Nnode);
      for (int i = 0; i < greedy->Nnode; i++)
      {
        int Nmembers = (*node)[i]->members.size();
        for (int j = 0; j < Nmembers; j++)
        {
          moveTo[(*node)[i]->members[j]] = i;
        }
      }
      for (int i = 0; i < greedy->Nnode; i++)
        delete (*node)[i];
      delete[] * node;

      greedy->Nnode = Nnode;
      greedy->Nmod = Nnode;
      greedy->Ndanglings = 0;
      greedy->node = rpt_node;
      greedy->calibrate();
      greedy->determineMove(moveTo);
      *node = rpt_node;
    }
    else
    {
      if (!silent)
        std::cout << "Iteration " << iteration + 1 << ", moving " << Nnode << " nodes, looping ";
    }

    double oldCodeLength;
    do
    {
      oldCodeLength = greedy->codeLength;
      bool moved = true;
      int Nloops = 0;
      int count = 0;
      while (moved) {
        moved = false;
        double inner_oldCodeLength = greedy->codeLength;
        greedy->move(moved);
        Nloops++;
        count++;
        if (greedy->codeLength - inner_oldCodeLength < 1.0e-10)
          moved = false;

        if (count == 10) {
          greedy->tune();
          count = 0;
        }
      }

      greedy->level(node, true);

      if (!silent)
        std::cout << Nloops << " ";
    } while (oldCodeLength - greedy->codeLength > 1.0e-10);

    iteration++;
    if (!silent)
      std::cout << "times between mergings to code length " << greedy->codeLength / log(2.0) << " in " << greedy->Nmod << " modules." << std::endl;
  } while (outer_oldCodeLength - greedy->codeLength > 1.0e-10);

  for (int i = 0; i < Nnode; i++)
    delete cpy_node[i];
  delete[] cpy_node;
}

void repeated_partition(MTRand *R, Node ***node, GreedyBase *greedy, bool silent, int Ntrials)
{
  double shortestCodeLength = 1000.0;
  int Nnode = greedy->Nnode;
  std::vector<int> cluster(Nnode);

  for (int trial = 0; trial < Ntrials; trial++)
  {
    if (!silent)
      std::cout << "Attempt " << trial + 1 << "/" << Ntrials << std::endl;

    Node **cpy_node = new Node*[Nnode];
    for (int i = 0; i < Nnode; i++) {
      cpy_node[i] = new Node();
      cpyNode(cpy_node[i], (*node)[i]);
    }

    greedy->Nnode = Nnode;
    greedy->Nmod = Nnode;
    greedy->Ndanglings = 0;
    greedy->node = cpy_node;
    greedy->calibrate();

    partition(R, &cpy_node, greedy, silent);

    if (greedy->codeLength < shortestCodeLength)
    {
      shortestCodeLength = greedy->codeLength;

      // Store best partition
      for (int i = 0; i < greedy->Nnode; i++) {
        for (auto mem = cpy_node[i]->members.begin(); mem != cpy_node[i]->members.end(); ++mem)
        {
          cluster[(*mem)] = i;
        }
      }
    }

    for (int i = 0; i < greedy->Nnode; i++)
    {
      delete cpy_node[i];
    }
    delete[] cpy_node;
  }

  // Commit best partition
  greedy->Nnode = Nnode;
  greedy->Nmod = Nnode;
  greedy->Ndanglings = 0;
  greedy->node = *node;
  greedy->calibrate();
  greedy->determineMove(cluster);
  greedy->level(node, true);
}

void printTree(std::string s, std::multimap<double, treeNode, std::greater<>>::iterator it_tM, std::ofstream *outfile, bool flip)
{
  if (!it_tM->second.nextLevel.empty())
  {
    int i = 1;
    for (auto it = it_tM->second.nextLevel.begin(); it != it_tM->second.nextLevel.end(); ++it)
    {
      std::string cpy_s(s + to_string(i) + ":");
      printTree(cpy_s, it, outfile, flip);
      i++;
    }
  }
  else
  {
    int i = 1;
    for (auto mem = it_tM->second.members.begin(); mem != it_tM->second.members.end(); ++mem)
    {
      if (flip)
      {
        std::string cpy_s(s + to_string(i) + " \"" + mem->second.second + "\" " + to_string(mem->first));
        *outfile << cpy_s << std::endl;
      }
      else
      {
        std::string cpy_s(s + to_string(i) + " " + to_string(mem->first) + " \"" + mem->second.second + "\"");
        *outfile << cpy_s << std::endl;
      }
      i++;
    }
  }
}
