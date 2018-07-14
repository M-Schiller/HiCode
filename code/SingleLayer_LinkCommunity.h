#pragma once

#ifndef SINGLELAYER_LINKCOMMUNITY_H
#define SINGLELAYER_LINKCOMMUNITY_H

#include <string>
#include <iostream>

class SingleLayer_LinkCommunity : public SingleLayer_Method
{
public:
  bool checkRequire() override
  {
    if (!checkRequiredFile(config["LinkCommunity_Dir"] + "calcAndWrite_Jaccards")) return false;
    if (!checkRequiredFile(config["LinkCommunity_Dir"] + "callCluster")) return false;
    return true;
  }

  void generateCommunity(std::string graphFile, std::string communityFile, int truth_NC)
  {
    std::string LCGraphFile = config["TMP_DIR"] + "LC.graph";
    std::cout << graphFile << std::endl;
    std::string JaccardsFile = config["TMP_DIR"] + "LC.jaccs";
    //cout<<graphFile<< std::endl;
    Graph g;
    g.load(graphFile);
    if (config["WeightedGraph"] != "TRUE")
      g.sampleByWeights();
    g.savePairs(LCGraphFile);
    const std::string LCCommunityFile = config["TMP_DIR"] + "LC.cluster";

    systemCall(config["LinkCommunity_Dir"] + "calcAndWrite_Jaccards " + LCGraphFile + " " + JaccardsFile);
    systemCall(config["LinkCommunity_Dir"] + "callCluster " + LCGraphFile + " " + JaccardsFile + " " + LCCommunityFile + " threshhold");

    //systemCall(config["LinkCommunity_Dir"]+"calcJaccards "+LCGraphFile+" "+JaccardsFile);
    //systemCall(config["LinkCommunity_Dir"]+"clusterJaccards "+LCGraphFile+" "+JaccardsFile+" "+LCCommunityFile+" "+config["TMP_DIR"]+"LC.stats"+" "+config["LinkCommunity_Thres"]);

    Community comms;
    comms.loadLinkCommunity(LCCommunityFile);
    comms.save(communityFile);

    puts("LC Finished");
    systemCall("rm -f " + config["TMP_DIR"] + "*");
  }
};
#endif // SINGLELAYER_LINKCOMMUNITY_H
