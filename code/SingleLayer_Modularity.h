#pragma once

#ifndef SINGLELAYER_MODULARITY_H
#define SINGLELAYER_MODULARITY_H

#include <string>
#include <iostream>

class SingleLayer_Modularity : public SingleLayer_Method {
public:
  bool checkRequire() {
    if (!checkRequiredFile(config["Modularity_Dir"] + "convert")) return false;
    if (!checkRequiredFile(config["Modularity_Dir"] + "community")) return false;
    if (!checkRequiredFile(config["Modularity_Dir"] + "hierarchy")) return false;
    return true;
  }
  void generateCommunity(std::string graphFile, std::string communityFile, int truth_NC)
  {
#ifdef DEBUG
    std::cout << "into Modularity SingleLayer" << std::endl;
#endif
    std::string singleEdgeGraphFile = config["TMP_DIR"] + "singleEdgeGraph";
    std::string modularityGraphFile = config["TMP_DIR"] + "ModularityGraph";
    std::string modularityWeightsFile = modularityGraphFile + ".weights";
    std::string modularityCommunityFile = config["TMP_DIR"] + "ModularityCommunity";
    Graph g;
    g.load(graphFile);
    if (config["WeightedGraph"] != "TRUE")
      g.sampleByWeights();
    g.saveSingleEdge(singleEdgeGraphFile);
    systemCall((config["Modularity_Dir"] + "convert -i " + singleEdgeGraphFile
      + " -o " + modularityGraphFile + " -w " + modularityWeightsFile).c_str());
    systemCall((config["Modularity_Dir"] + "community " + modularityGraphFile
      + " -l -1 -w " + modularityWeightsFile + " > " + modularityCommunityFile + " ").c_str());
    //systemCall(("cp "+modularityCommunityFile+" "+communityFile+".Modularity").c_str());
    Community comms;
    comms.loadModularity(modularityCommunityFile);
#ifdef PAUSEDEBUG
    getchar();
#endif
    comms.save(communityFile);
#ifdef PAUSEDEBUG
    getchar();
#endif
    puts("Modularity Finished");
  }
};
#endif // SINGLELAYER_MODULARITY_H
