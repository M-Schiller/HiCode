#include "infomap.h"

unsigned stou(char *s) {
  return std::strtoul(s, (char **)NULL, 10);
}

void partition(MTRand *R, Node ***node, GreedyBase *greedy, bool silent);
void repeated_partition(MTRand *R, Node ***node, GreedyBase *greedy, bool silent, int Ntrials);
void printTree(const std::string& s, std::multimap<double, treeNode, std::greater<>>::iterator it_tM, std::ofstream *outfile, bool flip);

// Call: trade <seed> <Ntries>
int main(int argc, char *argv[])
{
  if (argc != 4)
  {
    std::cout << "Call: ./infomap <seed> <network.net> <# attempts>" << std::endl;
    exit(-1);
  }

  const int Ntrials = std::atoi(argv[3]);  // Set number of partition attempts
  std::string infile = std::string(argv[2]);
  const std::string networkName(infile.begin(), infile.begin() + infile.find(".net"));
  std::string line;
  std::string buf;

  MTRand *R = new MTRand(stou(argv[1]));

  /* Read network in Pajek format with nodes ordered 1, 2, 3, ..., N,            */
  /* each undirected link occurring only once, and link weights > 0.             */
  /* (if a link is defined more than once, weights are aggregated)               */
  /* For more information, see http://vlado.fmf.uni-lj.si/pub/networks/pajek/.   */
  /* Example network with three nodes and                                        */
  /* three undirected weighted links:                                            */
  /* *Vertices 3                                                                 */
  /* 1 "Name of first node"                                                      */
  /* 2 "Name of second node"                                                     */
  /* 3 "Name of third node"                                                      */
  /* *Edges 3                                                                    */
  /* 1 2 1.0                                                                     */
  /* 1 3 3.3                                                                     */
  /* 2 3 2.2                                                                     */

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
    else
    {
      ss.clear();
      ss.str(line);
      ss >> buf;
      if (buf == "*Vertices"
        || buf == "*vertices"
        || buf == "*VERTICES")
      {
        ss >> buf;
        Nnode = std::stoi(buf);
      }
      else
      {
        std::cout << "the network file is not in Pajek format...exiting" << std::endl;
        exit(-1);
      }
    }
  }

  auto *nodeNames = new std::string[Nnode];

  // Read node names, assuming order 1, 2, 3, ...
  for (int i = 0; i < Nnode; i++)
  {
    getline(net, line);
    const int nameStart = line.find_first_of('\"');
    const int nameEnd = line.find_last_of('\"');
    if (nameStart < nameEnd)
    {
      nodeNames[i] = std::string(line.begin() + nameStart + 1, line.begin() + nameEnd);
    }
    else
    {
      ss.clear();
      ss.str(line);
      ss >> buf;
      ss >> nodeNames[i];
    }
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
  // Read links in format "from to weight", for example "1 3 2" (all integers) and each undirected link only ones (weight is optional).
  while (getline(net, line))
  {
    ss.clear();
    ss.str(line);
    ss >> buf;
    int linkEnd1 = std::stoi(buf);
    ss >> buf;
    int linkEnd2 = std::stoi(buf);
    buf.clear();
    ss >> buf;
    double linkWeight;
    if (buf.empty()) // If no information
      linkWeight = 1.0;
    else
      linkWeight = std::stof(buf);

    linkEnd1--; // Nodes start at 1, but C++ arrays at 0.
    linkEnd2--;

    if (linkEnd2 < linkEnd1)
    {
      int tmp = linkEnd1;
      linkEnd1 = linkEnd2;
      linkEnd2 = tmp;
    }

    // Aggregate link weights if they are definied more than once
    auto fromLink_it = Links.find(linkEnd1);
    if (fromLink_it == Links.end()) // new link
    {
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
  double totalDegree = 0.0;
  std::vector<double> degree(Nnode);
  auto node = new Node*[Nnode];

  for (int i = 0; i < Nnode; i++)
  {
    node[i] = new Node(i);
    degree[i] = 0.0;
  }

  int NselfLinks = 0;
  for (auto& link : Links)
  {
    for (auto toLink_it = link.second.begin(); toLink_it != link.second.end(); ++toLink_it)
    {
      int from = link.first;
      int to = toLink_it->first;
      double weight = toLink_it->second;
      if (weight > 0.0)
      {
        if (from == to)
        {
          NselfLinks++;
        }
        else {
          node[from]->links.emplace_back(to, weight);
          node[to]->links.emplace_back(from, weight);
          totalDegree += 2 * weight;
          degree[from] += weight;
          degree[to] += weight;
        }
      }
    }
  }
  if (NselfLinks > 0)
  {
    std::cout << ", ignoring " << NselfLinks << " self link(s)." << std::endl;
  }
  else
  {
    std::cout << ")" << std::endl;
  }

  //Swap maps to free memory
  for (auto& Link : Links)
  {
    std::map<int, double>().swap(Link.second);
  }
  std::map<int, std::map<int, double>>().swap(Links);

  GreedyBase * greedy = new Greedy(R, Nnode, totalDegree, node);
  greedy->initiate();

  const double uncompressedCodeLength = -(greedy->nodeDegree_log_nodeDegree);

  std::cout << "Now partition the network:" << std::endl;
  repeated_partition(R, &node, greedy, false, Ntrials);
  const int Nmod = greedy->Nnode;
  std::cout << "Done! Code length " << greedy->codeLength << " in " << Nmod << " modules." << std::endl;
  std::cout << "Compressed by " << 100.0*(1.0 - greedy->codeLength / uncompressedCodeLength) << " percent." << std::endl;

  // Order modules by size
  std::multimap<double, treeNode, std::greater<>> treeMap;
  for (int i = 0; i < greedy->Nnode; i++)
  {
    const int Nmembers = node[i]->members.size();
    treeNode tmp_tN;
    auto it_tM = treeMap.emplace(node[i]->degree / totalDegree, tmp_tN);
    for (int j = 0; j < Nmembers; j++)
    {
      it_tM->second.members.emplace(degree[node[i]->members[j]] / totalDegree, std::make_pair(node[i]->members[j], nodeNames[node[i]->members[j]]));
    }
  }

  // Order links by size
  std::multimap<double, std::pair<int, int>, std::greater<>> sortedLinks;
  for (int i = 0; i < Nmod; i++)
  {
    for (unsigned j = 0; j < node[i]->links.size(); j++)
    {
      if (i <= node[i]->links[j].first)
      {
        sortedLinks.emplace(node[i]->links[j].second / totalDegree, std::make_pair(i + 1, node[i]->links[j].first + 1));
      }
    }
  }

  //Print partition in format "module:rank size name"
  std::ostringstream oss;
  oss << networkName << ".tree";
  std::ofstream outfile;
  outfile.open(oss.str().c_str());
  outfile << "# Code length " << greedy->codeLength << " in " << Nmod << " modules." << std::endl;
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

  // Print partitions in Pajek's .clu format
  std::vector<int> clusterVec(Nnode);
  int clusterNr = 0;
  for (auto& mod : treeMap)
  {
    for (auto& member : mod.second.members)
    {
      clusterVec[member.second.first] = clusterNr;
    }
    clusterNr++;
  }
  oss.str("");
  oss << networkName << ".clu";
  outfile.open(oss.str().c_str());
  outfile << "*Vertices " << Nnode << "\x0D\x0A";

  for (int i = 0; i < Nnode; i++)
  {
    outfile << clusterVec[i] + 1 << "\x0D\x0A";
  }
  outfile.close();

  // Print map in Pajek's .net format (links sorted in descending order)
  oss.str("");
  oss << networkName << "_map.net";
  outfile.open(oss.str().c_str());
  outfile << "*Vertices " << Nmod << "\x0D\x0A";
  for (int i = 0; i < Nmod; i++)
  {
    outfile << i + 1 << " \"" << i + 1 << "\"" << "\x0D\x0A";
  }
  outfile << "*Edges " << sortedLinks.size() << "\x0D\x0A";
  for (auto& sortedLink : sortedLinks)
  {
    outfile
      << "  " << sortedLink.second.first
      << " " << sortedLink.second.second
      << " " << 1.0 * sortedLink.first / totalDegree
      << "\x0D\x0A";
  }
  outfile.close();

  // Print size of modules in Pajek's .vec format
  oss.str("");
  oss << networkName << "_map.vec";
  outfile.open(oss.str().c_str());
  outfile << "*Vertices " << Nmod << "\x0D\x0A";
  for (int i = 0; i < Nmod; i++)
  {
    outfile << 1.0*node[i]->degree / totalDegree << "\x0D\x0A";
  }
  outfile.close();

  // Print map in .map format for the Map Generator at www.mapequation.org
  oss.str("");
  oss << networkName << ".map";
  outfile.open(oss.str().c_str());
  outfile << "# modules: " << Nmod << std::endl;
  outfile << "# modulelinks: " << sortedLinks.size() << std::endl;
  outfile << "# nodes: " << Nnode << std::endl;
  outfile << "# links: " << Nlinks << std::endl;
  outfile << "# codelength: " << greedy->codeLength << std::endl;
  outfile << "*Undirected" << std::endl;
  outfile << "*Modules " << Nmod << std::endl;
  k = 0;

  for (auto& it : treeMap)
  {
    outfile << k + 1 << " \"" << it.second.members.begin()->second.second << "\" " << it.first << " " << node[k]->exit / totalDegree << std::endl;
    k++;
  }
  outfile << "*Nodes " << Nnode << std::endl;
  k = 1;
  for (auto it = treeMap.begin(); it != treeMap.end(); ++it)
  {
    std::string s;
    s.append(to_string(k));
    s.append(":");
    printTree(s, it, &outfile, true);
    k++;
  }
  outfile << "*Links " << sortedLinks.size() << std::endl;
  for (auto& sortedLink : sortedLinks)
  {
    outfile << sortedLink.second.first << " " << sortedLink.second.second << " " << 1.0 * sortedLink.first << std::endl;
  }
  outfile.close();

  delete[] nodeNames;
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
  do
  {
    outer_oldCodeLength = greedy->codeLength;

    if (iteration > 0
      && (iteration % 2 == 0)
      && greedy->Nnode > 1)  // Partition the partition
    {
      if (!silent)
      {
        std::cout << "Iteration " << iteration + 1 << ", moving " << std::flush;
      }

      auto rpt_node = new Node*[Nnode];
      for (int i = 0; i < Nnode; i++)
      {
        rpt_node[i] = new Node();
        cpyNode(rpt_node[i], cpy_node[i]);
      }
      std::vector<int> subMoveTo(Nnode);
      std::vector<int> moveTo(Nnode);
      int subModIndex = 0;

      for (int i = 0; i < greedy->Nnode; i++)
      {
        const int sub_Nnode = (*node)[i]->members.size();

        if (sub_Nnode > 1)
        {
          Node **sub_node = new Node*[sub_Nnode];
          std::set<int> sub_mem;
          for (int j = 0; j < sub_Nnode; j++)
            sub_mem.insert((*node)[i]->members[j]);
          auto it_mem = sub_mem.begin();
          int *sub_renumber = new int[Nnode];
          int *sub_rev_renumber = new int[sub_Nnode];
          double totalDegree = 0.0;
          for (int j = 0; j < sub_Nnode; j++)
          {
            //    fprintf(stderr,"%d %d\n",j,(*it_mem));
            int orig_nr = (*it_mem);
            int orig_Nlinks = cpy_node[orig_nr]->links.size(); // ERROR HERE
            sub_renumber[orig_nr] = j;
            sub_rev_renumber[j] = orig_nr;
            sub_node[j] = new Node(j);
            for (int k = 0; k < orig_Nlinks; k++)
            {
              int orig_link = cpy_node[orig_nr]->links[k].first;
              int orig_link_newnr = sub_renumber[orig_link];
              double orig_weight = cpy_node[orig_nr]->links[k].second;
              if (orig_link < orig_nr)
              {
                if (sub_mem.find(orig_link) != sub_mem.end())
                {
                  sub_node[j]->links.emplace_back(orig_link_newnr, orig_weight);
                  sub_node[orig_link_newnr]->links.emplace_back(j, orig_weight);
                  totalDegree += 2.0*orig_weight;
                }
              }
            }
            ++it_mem;
          }

          GreedyBase * sub_greedy = new Greedy(R, sub_Nnode, totalDegree, sub_node);
          sub_greedy->initiate();
          partition(R, &sub_node, sub_greedy, true);
          for (int j = 0; j < sub_greedy->Nnode; j++)
          {
            int Nmembers = sub_node[j]->members.size();
            for (int k = 0; k < Nmembers; k++)
            {
              subMoveTo[sub_rev_renumber[sub_node[j]->members[k]]] = subModIndex;
            }
            moveTo[subModIndex] = i;
            subModIndex++;
            delete sub_node[j];
          }

          delete[] sub_node;
          delete sub_greedy;
          delete[] sub_renumber;
          delete[] sub_rev_renumber;
        }
        else
        {
          subMoveTo[(*node)[i]->members[0]] = subModIndex;
          moveTo[subModIndex] = i;

          subModIndex++;
        }
      }

      for (int i = 0; i < greedy->Nnode; i++)
      {
        delete (*node)[i];
      }
      delete[](*node);

      greedy->Nnode = Nnode;
      greedy->Nmod = Nnode;
      greedy->node = rpt_node;
      greedy->initiate();
      greedy->determineMove(subMoveTo);
      greedy->level(node, false);
      greedy->determineMove(moveTo);
      (*node) = rpt_node;

      outer_oldCodeLength = greedy->codeLength;

      if (!silent)
        std::cout << greedy->Nnode << " modules, looping " << std::flush;
    }
    else if (iteration > 0)
    {
      if (!silent)
        std::cout << "Iteration " << iteration + 1 << ", moving " << Nnode << " nodes, looping " << std::flush;

      const auto rpt_node = new Node*[Nnode];
      for (int i = 0; i < Nnode; i++)
      {
        rpt_node[i] = new Node();
        cpyNode(rpt_node[i], cpy_node[i]);
      }

      std::vector<int> moveTo(Nnode);
      for (int i = 0; i < greedy->Nnode; i++)
      {
        const int Nmembers = (*node)[i]->members.size();
        for (int j = 0; j < Nmembers; j++)
        {
          moveTo[(*node)[i]->members[j]] = i;
        }
      }

      for (int i = 0; i < greedy->Nnode; i++)
      {
        delete (*node)[i];
      }
      delete[](*node);

      greedy->Nnode = Nnode;
      greedy->Nmod = Nnode;
      greedy->node = rpt_node;
      greedy->initiate();
      greedy->determineMove(moveTo);

      (*node) = rpt_node;
    }
    else
    {
      if (!silent)
        std::cout << "Iteration " << iteration + 1 << ", moving " << Nnode << " nodes, looping " << std::flush;
    }

    double oldCodeLength;
    do
    {
      oldCodeLength = greedy->codeLength;
      bool moved = true;
      int Nloops = 0;
      int count = 0;
      while (moved)
      {
        moved = false;
        const double inner_oldCodeLength = greedy->codeLength;
        greedy->move(moved);
        Nloops++;
        count++;
        if (inner_oldCodeLength - greedy->codeLength < 1.0e-10)
        {
          moved = false;
        }

        if (count == 10)
        {
          greedy->tune();
          count = 0;
        }
        // 	if(!silent){
        // 	  cerr << Nloops;
        // 	  int loopsize = to_string(Nloops).length();
        // 	  for(int i=0;i<loopsize;i++)
        // 	    cerr << "\b";
        // 	}
      }
      greedy->level(node, true);

      if (!silent)
      {
        std::cout << Nloops << " " << std::flush;
      }
    } while (oldCodeLength - greedy->codeLength > 1.0e-10);

    iteration++;
    if (!silent)
    {
      std::cout << "times between mergings to code length " << greedy->codeLength << " in " << greedy->Nmod << " modules." << std::endl;
    }
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
    greedy->node = cpy_node;
    greedy->initiate();

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
  greedy->node = (*node);
  greedy->initiate();
  greedy->determineMove(cluster);
  greedy->level(node, true);
}

void printTree(const std::string& s, std::multimap<double, treeNode, std::greater<>>::iterator it_tM, std::ofstream *outfile, bool flip)
{
  if (!it_tM->second.nextLevel.empty())
  {
    int i = 1;
    for (auto it = it_tM->second.nextLevel.begin(); it != it_tM->second.nextLevel.end(); ++it)
    {
      const std::string cpy_s(s + to_string(i) + ":");
      printTree(cpy_s, it, outfile, flip);
      i++;
    }
  }
  else
  {
    int i = 1;
    for (auto& member : it_tM->second.members)
    {
      if (flip)
      {
        const std::string cpy_s(s + to_string(i) + " \"" + member.second.second + "\" " + to_string(member.first));
        (*outfile) << cpy_s << std::endl;
      }
      else
      {
        const std::string cpy_s(s + to_string(i) + " " + to_string(member.first) + " \"" + member.second.second + "\"");
        (*outfile) << cpy_s << std::endl;
      }
      i++;
    }
  }
}
