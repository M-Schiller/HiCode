#pragma once

#ifndef SINGLELAYER_OSLOM_H
#define SINGLELAYER_OSLOM_H

#include <string>
#include <cstdlib>

class SingleLayer_OSLOM : public SingleLayer_Method
{
public:
  bool checkRequire() override
  {
    return checkRequiredFile(config["OSLOM_Dir"] + "oslom_undir");
  }

  void generateCommunity(std::string graphFile, std::string communityFile, int truth_NC)
  {
    std::string oslomGraphFile = config["TMP_DIR"] + "oslomGraph";
    Graph g;
    g.load(graphFile);
    if (config["WeightedGraph"] != "TRUE")
      g.sampleByWeights();
    g.savePairs(oslomGraphFile);

    systemCall(config["OSLOM_Dir"] + "oslom_undir -f " + oslomGraphFile + " -uw");
    std::string OSLOMCommunityFile = oslomGraphFile + "_oslo_files/tp";
    std::string communityFilename = communityFile.substr(0, communityFile.size() - 4);
    std::string suffix = communityFile.substr(communityFile.size() - 4);
    Community comms;
    comms.loadOSLOM(OSLOMCommunityFile);
    comms.save(communityFilename + int2str(0) + suffix);
    int abs_diff = abs(comms.NC - truth_NC);
    int min_diff = abs_diff;
    int sim_comm = 0;
    for (int i = 1;; i++)
    {
      if (!checkFileExist(OSLOMCommunityFile + int2str(i))) break;
      Community comms;
      comms.loadOSLOM(OSLOMCommunityFile + int2str(i));
      comms.save(communityFilename + int2str(i) + suffix);
      abs_diff = abs(comms.NC - truth_NC);
      if (abs_diff < min_diff) {
        min_diff = abs_diff;
        sim_comm = i;
      }
    }
    std::cout << "choosing Hirachy : " << sim_comm << std::endl;
    systemCall("cp " + communityFilename + int2str(sim_comm) + suffix + " " + communityFile);
    puts("OSLOM Finished");
  }
};
#endif // SINGLELAYER_OSLOM_H
