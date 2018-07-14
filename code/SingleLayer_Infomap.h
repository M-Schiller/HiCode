#pragma once

#ifndef SINGLELAYER_INFOMAP_H
#define SINGLELAYER_INFOMAP_H

class SingleLayer_Infomap : public SingleLayer_Method
{
public:
  bool checkRequire() override
  {
    return checkRequiredFile(config["Infomap_Dir"] + "Infomap");
  }

  void generateCommunity(std::string graphFile, std::string communityFile, int truth_NC)
  {
    std::string infomapGraphFile = config["TMP_DIR"] + "InfomapGraph";
    std::string infomapCommunityDir = config["TMP_DIR"];
    std::string infomapCommunityFile = config["TMP_DIR"] + "InfomapGraph.clu";

    Graph g;
    g.load(graphFile);
    g.saveSingleEdge(infomapGraphFile);

    std::string tags = " -i link-list --clu -s " + int2str(rand());
    systemCall((config["Infomap_Dir"] + "Infomap "
      + infomapGraphFile + " " + infomapCommunityDir + tags).c_str());

    Community comms{};
    comms.loadInfomap(infomapCommunityFile);
    comms.save(communityFile);
  }
};
#endif // SINGLELAYER_INFOMAP_H
