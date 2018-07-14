#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <string>
#include <utility>
#include <vector>
#include <numeric>
#include <queue>
#include <set>
#include <map>
#include <sstream>
#include <cctype>
#include <fstream>
//#include <ctime>
#include <sys/time.h>
#include <unistd.h>
//#include <assert.h>

std::map<std::string, std::string> config;

#define overlap //using overlap headfiles
//#define DEBUG

#include "basic.h"
#include "Graph.h"
#include "Community.h"
#include "SingleLayer_Method.h"
#include "SingleLayer_Modularity.h"
#include "SingleLayer_Infomap.h"
#include "SingleLayer_OSLOM.h"
#include "SingleLayer_LinkCommunity.h"
#include "Framework.h"
#include "Framework_Remove.h"
#include "Framework_ReducePP.h"

#ifdef overlap
#include "Metrics_Overlapping.h"
#else
#include "Metrics.h"
#endif // overlap
#define lowbit(x) ((x)&(-(x)))
#define sqr(x) ((x)*(x))

#define LINESEP "###########################\n"

SingleLayer_Method *SLM, *SLM_Modularity;
Framework *FW1;
Framework *FW2;

void loadConfig(char* filename)
{
  std::ifstream in_file(filename);
  std::string line;

  while (std::getline(in_file, line))
  {
    if (line.find_first_not_of(' ') == line.npos)
    {
      continue;
    }
    std::string key, value;
    std::istringstream tmp(line);
    tmp >> key >> value;

    config[key] = value;
    std::cout << "K: " << key << " V: " << value << std::endl;
  }
  in_file.close();

  exit(0);
}

SingleLayer_Method* getSingleLayerMethod(const std::string& name)
{
  if (name == "Modularity")
    return new SingleLayer_Modularity();
  if (name == "Infomap")
    return new SingleLayer_Infomap();
  if (name == "OSLOM")
    return new SingleLayer_OSLOM();
  if (name == "LinkCommunity")
    return new SingleLayer_LinkCommunity();
  return new SingleLayer_Method();
}

bool checkRequire()
{
  if (!checkRequiredFile(config["DATA_DIR"] + "graph"))
  {
    return false;
  }
  //cout<<config["SingleLayer_Method"]<<":here"<< std::endl;
  SLM_Modularity = new SingleLayer_Modularity();
  SLM_Modularity->setConfig(config);

  if (config["SingleLayer_Method"] == "Modularity")
    SLM = new SingleLayer_Modularity();
  else if (config["SingleLayer_Method"] == "Infomap")
    SLM = new SingleLayer_Infomap();
  else if (config["SingleLayer_Method"] == "LinkCommunity")
    SLM = new SingleLayer_LinkCommunity();
  else if (config["SingleLayer_Method"] == "OSLOM")
    SLM = new SingleLayer_OSLOM();
  else
  {
    puts("UNKNOWN SingleLayer_Method");
    return false;
  }

  SLM->setConfig(config);
  if (!SLM->checkRequire())
  {
    return false;
  }

  if (config["Framework"] == "Remove")
    FW1 = FW2 = new Framework_Remove();
  else if (config["Framework"] == "Reduce++")
    FW1 = FW2 = new Framework_ReducePP();
  else
  {
    puts("UNKNOWN Framework");
  }
  FW1->setConfig(config);
  FW2->setConfig(config);
  if (!FW1->checkRequire())
  {
    return false;
  }
  if (!FW2->checkRequire())
  {
    return false;
  }

  return true;
}

void generateLayerGraph(
  const std::string& graphFile,
  std::vector<std::string> otherCommFiles,
  const std::string& outputFile,
  bool isFirst)
{
  Graph cur;
  std::cout << graphFile << std::endl;
  if (!checkFileExist(graphFile))
    std::cout << graphFile << " not exist" << std::endl;
  cur.load(std::move(graphFile));
  std::vector<Community> comms;
  comms.clear();
  for (int i = 0; i < otherCommFiles.size(); i++)
  {
    Community comm;
    std::cout << i << " : " << otherCommFiles[i] << std::endl;
    comm.load(otherCommFiles[i]);
    if (config.find("Framework_CommunitySizeThres") != config.end())
      comm.removeSmallComms(std::stoi(config["Framework_CommunitySizeThres"]));
    comms.push_back(comm);
  }
  Graph layerGraph;
  if (isFirst)
  {
    std::cout << "isFirst" << std::endl;
    layerGraph = FW2->calcLayerGraph(cur, comms);
  }
  else
  {
    std::vector<Community> c1;
    c1.clear();
    c1.push_back(comms[0]);
    for (int i = 1; i < comms.size(); i++) comms[i - 1] = comms[i];
    comms.pop_back();
    std::cout << "notFirst 1" << std::endl;
    layerGraph = FW1->calcLayerGraph(cur, c1);
    std::cout << "notFirst 2" << std::endl;
    std::cout << "comms.size() = " << comms.size();
    layerGraph = FW1->calcLayerGraph(layerGraph, comms);
    std::cout << "notFirst 3" << std::endl;
  }
  layerGraph.save(std::move(outputFile));
}

void generateLayerGraphAll(
  const std::string& graphFile,
  std::vector<std::string> otherCommFiles,
  const std::string& outputFile,
  bool isFirst)
{
  Graph cur;
  cur.load(std::move(graphFile));
  std::vector<Community> comms;
  comms.clear();
  for (int i = 0; i < otherCommFiles.size(); i++)
  {
    Community comm;
    comm.load(otherCommFiles[i]);
    if (config.find("Framework_CommunitySizeThres") != config.end())
      comm.removeSmallComms(std::stoi(config["Framework_CommunitySizeThres"]));
    comms.push_back(comm);
  }
  Graph layerGraph;
  if (isFirst)
    layerGraph = FW2->calcLayerGraphAll(cur, comms);
  else
  {
    //vector<Community> c1;c1.clear();c1.push_back(comms[0]);
    //for (int i=1;i<comms.size();i++) comms[i-1]=comms[i];
    //comms.pop_back();
    //layerGraph=FW1->calcLayerGraphAll(cur,c1);
    layerGraph = FW2->calcLayerGraphAll(cur, comms);
  }
  layerGraph.save(std::move(outputFile));
}

void generateCommunity(const std::string& graphFile, const std::string& communityFile, int layerid, std::vector<int> Truth_NC)
{
  //Truth_NC[0] = -1, Truth_NC[i] = truth[i-1].NC, 1<=i<truth.size() //if there is k truths, Truth_NC.size() = k+1 and truth.size() = k
  //layerid [1, nlayer]
      /*if (layerid==1)
        SLM_Modularity->generateCommunity(graphFile,commnityFile, Truth_NC);
      else */
  std::cout << graphFile << " " << communityFile << std::endl;
  int truth_NC;
  if (layerid < Truth_NC.size())
  {
    truth_NC = Truth_NC[layerid];
  }
  else
  {
    /*int randid = ceil(getRand(0, Truth_NC.size()-1));
    std::cout<<"random choosing similar truth : "<<randid<< std::endl;
    truth_NC = Truth_NC[randid];*/
    std::vector<int>::iterator result = max_element(Truth_NC.begin(), Truth_NC.end());
    int truth_NC_idx = distance(Truth_NC.begin(), result);
    truth_NC = Truth_NC[truth_NC_idx];
    //getchar();
  }
  std::cout << "truth Number of C = " << truth_NC << std::endl;
  SLM->generateCommunity(std::move(communityFile), truth_NC);
}

std::vector<std::pair<std::string, std::string>> showCommunityStat(const std::string& filename)
{
  std::vector<std::pair<std::string, std::string>> stat;

  Community comms;
  comms.load(filename);
  stat.emplace_back("#Comm", int2str(comms.NC));
  stat.emplace_back("Avg Comm Size", double2str(comms.N / (double)comms.NC));
  stat.emplace_back("Comm Entropy", double2str(calcEntropy(comms)));
  std::vector<std::vector<int >> cid = comms.getPartitionIDOver();
  double counter = 0, OM_counter = 0;
  for (int cid_i = 1; cid_i < cid.size(); cid_i++)
  {
    if (cid[cid_i].size() != 0)
    {
      counter++;
      OM_counter += cid[cid_i].size();
    }
  }
  stat.emplace_back("Coverage", double2str(counter / (double)comms.N));
  stat.emplace_back("OM", double2str(OM_counter / (double)comms.N));
  return stat;
}

std::vector<std::pair<std::string, std::string>> showGraphStat(std::string filename)
{
  std::vector<std::pair<std::string, std::string>> stat;
  stat.clear();
  Graph g;
  g.load(filename);
  double W = 0;
  for (int i = 0; i < g.m_edges.size(); i++) W += g.m_edges[i].w;
  stat.emplace_back("Graph Density", double2str(W / (double)(g.N*g.N)));
  return stat;
}

void showConfig(FILE* fout, FILE* fout_tex)
{
  std::vector<std::pair<std::string, std::string>> data;
  data.insert(data.end(), config.begin(), config.end());

  Tex_Table(VSS2VVS(data), "Experiment Setting", fout_tex);
}

std::vector<std::string> truth, frameworks;
std::vector<std::pair<std::string, std::vector<std::vector<std::string>>>> iterationDetail;
std::map<std::string, std::vector<std::string>> iterativeStat;
std::vector<std::string> iterativeStatNames;
std::vector<bool> iterativeStatNames_graphic;

//for save max origal and layer
double maxmodOriG = 0;
double maxmodLayerG = 0;
double tmpmodOriG = 0;
double tmpmodLayerG = 0;
int maxOriIter = 0;
int maxLayIter = 0;

int nLayer;

bool isLayer(const std::string& layerName)
{
  return layerName.substr(0, 5) == "Layer";
}

void autoShrinkLayer()//自动确定层数的策略，已经废弃
{
  int cnt = 1;
  for (int i = 1; i <= nLayer; i++)
  {
    double ma = 0;
    for (int z = 1; z < cnt; z++)//cnt=1???????
    {
      ma = std::max(ma, calcNMI(config["RESULT_DIR"] + "Layer" + int2str(i) + ".gen", config["RESULT_DIR"] + "Layer" + int2str(z) + ".gen"));
    }
    if (ma < str2double(config["ShrinkLayerThres"]))
    {
      systemCall("cp " + config["RESULT_DIR"] + "Layer" + int2str(i) + ".gen " + config["RESULT_DIR"] + "Layer" + int2str(cnt) + ".gen");
      cnt++;
    }
  }
  nLayer = cnt - 1;
}

std::vector<std::vector<double>> preVScur;//for draw the NMI between the current and previous

std::pair<std::string, std::vector<std::vector<std::string>>> analyze(std::string name)
{
  std::cout << "analyze begin" << std::endl;
  std::vector<std::pair<std::string, std::string>> comms;
  comms.clear();
  for (int i = 0; i < truth.size(); i++)
  {
    comms.emplace_back(truth[i], config["DATA_DIR"] + truth[i]);
  }
  std::cout << "comes 1" << std::endl;
  for (int i = 1; i <= nLayer; i++)
  {
    comms.emplace_back("Layer" + int2str(i), config["RESULT_DIR"] + "Layer" + int2str(i) + ".gen");
  }
  std::cout << "comes 2" << std::endl;
  std::vector<std::vector<std::string>> result;
  result.clear();

  std::vector<double> _modOriG;
  _modOriG.clear();
  std::vector<double> _modLayerG;
  _modLayerG.clear();
  std::vector<double> _NMICL;
  _NMICL.clear();
  std::vector<double> _NMIPL;
  _NMIPL.clear();
  double score;

  for (int i = 0; i < comms.size(); i++)
  {
    std::cout << "comes for start " << int2str(i) << std::endl;
    std::vector<std::string> row;
    row.clear();
    row.push_back(comms[i].first);
    std::vector<std::pair<std::string, std::string>> res = showCommunityStat(comms[i].second);
    //vs Previous Layer
    if (checkFileExist(comms[i].second + ".last"))
    {
      res.emplace_back("NMI(Last Iteration)", double2str(score = calcNMI(comms[i].second, comms[i].second + ".last")));
      _NMIPL.push_back(score);
    }

    else
    {
      res.emplace_back("NMI(Last Iteration)", double2str(log(0)));
    }
    //getchar();
    //Comm VS Graph
    std::string originalGraphFile = config["DATA_DIR"] + "graph";
    std::string layerGraphFile = config["RESULT_DIR"] + comms[i].first + ".graph";
    //cout<<"main3 analyze here"<< std::endl;
    //cout<<"comm file:"<<comms[i].second<< std::endl;
    res.emplace_back("Modularity(Original)", double2str(score = calcModularity(originalGraphFile, comms[i].second)));
    //cout<<"main2 analyze here"<< std::endl;
    if (isLayer(comms[i].first))
    {
      _modOriG.push_back(score);
    }
    //cout<<"main3 analyze here"<< std::endl;
    res.emplace_back("Modularity(Layer Graph)", double2str(score = calcModularity(layerGraphFile, comms[i].second)));

    if (isLayer(comms[i].first))
      _modLayerG.push_back(score);
    //Comm VS Comm

    for (int j = 0; j < comms.size(); j++)
    {
      res.emplace_back("NMI(" + comms[j].first + ")", double2str(score = calcNMI(comms[i].second, comms[j].second)));
      //cout<<" Overlapping NMI = "<<score<< std::endl;
      if (j > i
        && isLayer(comms[i].first)
        && isLayer(comms[j].first))
        _NMICL.push_back(score);
    }

    if (i == 0)
    {
      std::vector<std::string> title;
      title.clear();
      title.push_back("");
      for (int z = 0; z < res.size(); z++)
        title.push_back(res[z].first);
      result.push_back(title);
    }
    for (int z = 0; z < res.size(); z++)
      row.push_back(res[z].second);
    result.push_back(row);
    std::cout << "comes for end " << int2str(i) << std::endl;
  }

  //for save max
  std::cout << "push_back start" << std::endl;
  tmpmodOriG = getAverage(_modOriG);
  tmpmodLayerG = getAverage(_modLayerG);
  preVScur.push_back(_NMIPL);
  std::string tmpstr;
  iterativeStat["Name"].push_back(name);
  for (int index = 0; index < nLayer; index++)
  {
    tmpstr = "Mod_Layer" + int2str(index + 1) + "_in_" + "OriG";
    iterativeStat[tmpstr].push_back(double2str(_modOriG[index]));
  }
  for (int index = 0; index < nLayer; index++)
  {
    tmpstr = "Mod_Layer" + int2str(index + 1) + "_in_" + "LayerG";
    iterativeStat[tmpstr].push_back(double2str(_modLayerG[index]));
  }
  for (int index = 0; index < nLayer; index++)
  {
    tmpstr = "NMI_Layer" + int2str(index + 1) + "_vs_" + "inPreIter";
    if (_NMIPL.empty())
      iterativeStat[tmpstr].push_back(double2str(getAverage(_NMIPL)));
    else
      iterativeStat[tmpstr].push_back(double2str(_NMIPL[index]));
  }
  iterativeStat["AvgMod_Layer_in_OriG"].push_back(double2str(getAverage(_modOriG)));
  iterativeStat["AvgMod_Layer_in_LayerG"].push_back(double2str(getAverage(_modLayerG)));
  iterativeStat["AvgNMI_Iter_vs_PrevIter"].push_back(double2str(getAverage(_NMIPL)));
  iterativeStat["AvgNMI_CrossLayer"].push_back(double2str(getAverage(_NMICL)));
  std::cout << "push_back end" << std::endl;
  return std::make_pair(name, transVVS(result));
}

void generateGraphic()
{
  FILE* fout = fopen("graphicData", "w");
  std::vector<std::vector<std::string>> data;
  data.clear();
  for (int i = 0; i < iterativeStatNames.size(); i++)
  {
    if (iterativeStatNames_graphic[i])
    {
      data.push_back(iterativeStat[iterativeStatNames[i]]);
    }
  }
  std::cout << "VVS start" << std::endl;
  data = transVVS(data);
  std::cout << "VVS stop" << std::endl;
  for (int i = 0; i < data.size(); i++)
  {
    for (int j = 0; j < data[i].size(); j++)
    {
      if (j) fprintf(fout, ",");
      fprintf(fout, "%s", data[i][j].c_str());
    }
    fprintf(fout, "\n");
  }
  fclose(fout);
  systemCall("Rscript tools/graphic/genGraphics.R");//systemCall("Rscript tools/graphic/genGraphics.R");
}

void showAnalyze()
{
  FILE* fout = fopen((config["RESULT_DIR"] + "analyze").c_str(), "w");
  FILE* fout_tex = fopen((config["RESULT_DIR"] + "analyze.tex").c_str(), "w");
  fprintf(fout_tex, "\\documentclass[4pt]{article}\n");
  fprintf(fout_tex, "\\usepackage{graphicx}");
  fprintf(fout_tex, "\\usepackage[papersize={30cm,50cm}]{geometry}");
  fprintf(fout_tex, "\\begin{document}\n");

  showConfig(fout, fout_tex);

  std::vector<std::vector<std::string>> expSummary;
  expSummary.clear();
  std::cout << iterativeStatNames.size() << std::endl;
  for (int i = 0; i < iterativeStatNames.size(); i++)
  {
    std::cout << iterativeStatNames[i] << std::endl;
    expSummary.push_back(iterativeStat[iterativeStatNames[i]]);
  }

  //Tex_Table(transVVS(expSummary),"Experiment Summary",fout_tex);

  generateGraphic();
  int includeG_i = 0;
  for (int i = 0; i < iterativeStatNames.size(); i++)
  {
    if (iterativeStatNames_graphic[i])
    {
      Tex_includeGraphics(iterativeStatNames[i] + ".pdf",
        iterativeStatNames[i], 5, fout_tex);
      fprintf(fout, "%s\n", iterativeStatNames[i].c_str());
      for (int z = 0; z < iterativeStat[iterativeStatNames[i]].size(); z++)
      {
        fprintf(fout, "%s\t", iterativeStat[iterativeStatNames[i]][z].c_str());
      }
      fprintf(fout, "\n");
      includeG_i++;
      if (includeG_i % 3 == 0)
      {
        fprintf(fout_tex, "\\clearpage\n");
      }
    }
  }

  int layers = std::stoi(config["Number_Of_Layers"]);
  for (int i = 1; i <= layers; ++i)
  {
    std::string str = "layer" + int2str(i) + "_vs_pre";
    Tex_includeGraphics(str + ".pdf",
      str.c_str(), 5, fout_tex);
  }
  for (int i = 1; i <= layers; ++i)
  {
    for (int t = 1; t <= truth.size(); t++)
    {
      std::string str2 = "truth" + int2str(t) + "_vs_layer" + int2str(i);
      //string str2 = truth[t-1]+"_vs_layer"+int2str(i);
      Tex_includeGraphics(str2 + ".pdf", str2.c_str(), 5, fout_tex);
    }
    fprintf(fout_tex, "\\clearpage\n");
  }

  //Tex_includeGraphics("deptconnectgen_vs_layer1.pdf", "deptconnect.gen_vs_layer1.pdf",5,fout_tex);

  int cnt = 0;
  for (int i = 0; i < iterationDetail.size(); i++)
  {
    if (i % 1 == 0)
    {
      cnt++;
      if (cnt % 5 == 0) fprintf(fout_tex, "\\clearpage\n");
      Tex_Table(iterationDetail[i].second, iterationDetail[i].first, fout_tex);
    }
  }
  Tex_Table(iterationDetail[maxOriIter].second, iterationDetail[maxOriIter].first, fout_tex);//maxOriIter
  Tex_Table(iterationDetail[maxLayIter].second, iterationDetail[maxLayIter].first, fout_tex);//maxLayIter
  fprintf(fout_tex, "\\end{document}\n");

  fclose(fout);
  fclose(fout_tex);

  systemCall("cat " + config["RESULT_DIR"] + "analyze");
  systemCall("xelatex " + config["RESULT_DIR"] + "analyze.tex");
  systemCall("cp analyze.pdf " + config["RESULT_DIR"] + "analyze.pdf");
}

void PrintCurAndPre()
{
  FILE *fout = fopen("Cur_vs_Pre", "w");
  int layers = std::stoi(config["Number_Of_Layers"]);
  for (int i = 0; i < layers; ++i)
  {
    std::string str = "layer" + int2str(i + 1) + "_vs_pre";
    if (i != layers - 1)
      str += ",";
    else
      str += "\n";
    fprintf(fout, str.c_str());
  }
  //cout<<"total: "<<preVScur.size()<< std::endl;
  for (int i = 0; i < preVScur.size(); i++)
  {
    for (int j = 0; j < preVScur[i].size(); j++)
    {
      if (j != preVScur[i].size() - 1) fprintf(fout, "%f,", preVScur[i][j]);
      else fprintf(fout, "%f", preVScur[i][j]);
      //cout<<" llss ";
    }
    fprintf(fout, "\n");
  }
  fclose(fout);
  //cout<<"here"<< std::endl;
  systemCall("Rscript tools/graphic/genCurVsPre.R");//systemCall("Rscript tools/graphic/genCurVsPre.R");
}

void PrintTruthAndLayer()
{
  FILE* fout = fopen("Truth_vs_layer", "w");
  int cnt = 0;
  std::vector<std::vector<std::string>> data;
  data.clear();
  int layers = std::stoi(config["Number_Of_Layers"]);
  for (int k = 0; k < layers; k++)
  {
    for (int t = 1; t <= truth.size(); t++)
    {
      std::vector<std::string>  vstmp;
      std::string str = "truth" + int2str(t) + "_vs_layer" + int2str(k + 1);
      //cout<<truth[t-1]<< std::endl;
      vstmp.push_back(str);
      //int index1=1;
      //double in1 = 0;
      //vector<vector<string>> tmp = iterationDetail[iterationDetail.size()-1].second;
      //找层数与groundtruth最像的；
      //for(int t=1; t<=truth.size(); ++t)
      //    if(str2double(tmp[9+k+truth.size()][t]) > in1)
      //    {
      //        in1 = str2double(tmp[9+k+truth.size()][t]);
      //        index1 = t;
      //    }
      //puts("hereim1");

      for (int i = 0; i < iterationDetail.size(); i++)//what the hell it is??????????
      {
        std::vector<std::vector<std::string>> table = iterationDetail[i].second;
        //puts("hereim2");
        vstmp.push_back(table[9 + k + truth.size()][t].c_str());//if you want to add table elements, do not forget to change the number 9.
        //puts("hereim3");
        //if(i!=iterationDetail.size()-1)
        // fprintf(fout,",");
      }
      //fprintf(fout,"\n");
      data.push_back(vstmp);
      //puts("hereim4");
    }
  }
  data = transVVS(data);
  for (int i = 0; i < data.size(); i++)
  {
    for (int j = 0; j < data[i].size(); j++)
    {
      if (j) fprintf(fout, ",");
      fprintf(fout, "%s", data[i][j].c_str());
    }
    fprintf(fout, "\n");
  }
  //puts("hereim5");
  fclose(fout);
  systemCall("Rscript tools/graphic/genTruth.R");//systemCall("Rscript tools/graphic/genCurVsPre.R");
}

int main(int argc, char** argv)
{
  if (argc != 2)
  {
    puts("Argument Error");
    return 0;
  }
  loadConfig(argv[1]);
  config["TMP_DIR"] = config["RESULT_DIR"] + "tmp/";

  systemCall("mkdir -p " + config["RESULT_DIR"] + "Iteration/");//new add
  systemCall(("rm -f -r " + config["RESULT_DIR"] + "*").c_str());
  systemCall("mkdir -p " + config["TMP_DIR"]);

  truth = split(config["Ground_Truth"]);
  frameworks = split(config["Frameworks"]);

  systemCall("mkdir -p " + config["RESULT_DIR"] + "BasicAlgs");
  systemCall("mkdir -p " + config["RESULT_DIR"] + "Time");
  //for save max
  systemCall("mkdir -p " + config["RESULT_DIR"] + "maxOriginal");
  systemCall("mkdir -p " + config["RESULT_DIR"] + "maxLayer");

  //read NC in truths
  std::vector<int> Truth_NC;
  Truth_NC.push_back(-1);//Truth_NC[0]=-1, Truth_NC[i]=truth[i-1].NC, 1<=i<=truth.size()
#ifdef DEBUG
  std::cout << "I'm here" << endl;
#endif
  for (int truth_i = 0; truth_i < truth.size(); truth_i++)
  {
    Community truth_com;
    truth_com.load(config["DATA_DIR"] + truth[truth_i]);
    std::cout << truth[truth_i] << " NC: " << truth_com.NC << std::endl;
    Truth_NC.push_back(truth_com.NC);
  }
  //time_t start, end;
  struct timeval start, end;

  double span = 0.0;

  //systemCall("cp "+config["RESULT_DIR"]+"Basictime.txt "+config["RESULT_DIR"]+"Time/");
  //for (;;);

  nLayer = std::stoi(config["Number_Of_Layers"]);

  systemCall("touch " + config["RESULT_DIR"] + "Time/" + config["SingleLayer_Method"] + "_Frameworktime.txt");
  std::ofstream outfile2;
  outfile2.open((config["RESULT_DIR"] + "Time/" + config["SingleLayer_Method"] + "_Frameworktime.txt").c_str(), std::ios::out);

  for (int nf = 0; nf < frameworks.size(); nf++)
  {
    std::cout << "running time of algorithm " + frameworks[nf] << std::endl;
    outfile2 << "running time of " + frameworks[nf] + " in seconds:" << std::endl;
    //start = time(NULL);
    gettimeofday(&start, NULL);
    config["Framework"] = frameworks[nf];
    if (!checkRequire()) continue;
    iterationDetail.clear();
    iterativeStat.clear();
    iterativeStatNames.clear();
    iterativeStatNames_graphic.clear();
    iterativeStatNames.push_back("Name");
    iterativeStatNames_graphic.push_back(false);
    iterativeStatNames.push_back("AvgMod_Layer_in_OriG");
    iterativeStatNames_graphic.push_back(true);
    iterativeStatNames.push_back("AvgMod_Layer_in_LayerG");
    iterativeStatNames_graphic.push_back(true);
    iterativeStatNames.push_back("AvgNMI_CrossLayer");
    iterativeStatNames_graphic.push_back(true);
    iterativeStatNames.push_back("AvgNMI_Iter_vs_PrevIter");
    iterativeStatNames_graphic.push_back(true);
    std::string tmpstr;
    for (int i = 1; i <= nLayer; i++)
    {
      tmpstr = "NMI_Layer" + int2str(i) + "_vs_" + "inPreIter";
      iterativeStatNames.push_back(tmpstr);
      iterativeStatNames_graphic.push_back(true);
      tmpstr = "Mod_Layer" + int2str(i) + "_in_" + "LayerG";
      iterativeStatNames.push_back(tmpstr);
      iterativeStatNames_graphic.push_back(true);
      tmpstr = "Mod_Layer" + int2str(i) + "_in_" + "OriG";
      iterativeStatNames.push_back(tmpstr);
      iterativeStatNames_graphic.push_back(true);
    }
    for (int i = 0; i < iterativeStatNames.size(); i++)
    {
      iterativeStat[iterativeStatNames[i]] = std::vector<std::string>(0);
      iterativeStat[iterativeStatNames[i]].push_back(iterativeStatNames[i]);
    }
    systemCall(("cp " + config["DATA_DIR"] + "graph " + config["RESULT_DIR"] + "Layer1.graph").c_str());
#ifdef DEBUG
    std::cout << "here i come" << endl;
#endif

    //Identification stage
    std::cout << "nLayer = " << nLayer << std::endl;
    for (int i = 1; i <= nLayer; i++)
    {
      std::cout << "i = " << i << " nLayer = " << nLayer << std::endl;
      if (i != 1)
      {
        std::vector< std::string> otherCommFiles;
        for (int j = 1; j < i; j++)
        {
          if (i != j)
          {
            otherCommFiles.push_back(config["RESULT_DIR"] + "Layer" + int2str(j) + ".gen");
          }
        }
        generateLayerGraph(config["DATA_DIR"] + "graph", otherCommFiles, config["RESULT_DIR"] + "Layer" + int2str(i) + ".graph", i == 1);
      }

#ifdef DEBUG
      std::cout << "here i come 2" << endl;
#endif
      generateCommunity(config["RESULT_DIR"] + "Layer" + int2str(i) + ".graph"
        , config["RESULT_DIR"] + "Layer" + int2str(i) + ".gen", i, Truth_NC);
#ifdef DEBUG
      getchar();
#endif
      if (config.find("Save_CommunitySizeThres") != config.end())
      {
        Community comm;
        comm.load(config["RESULT_DIR"] + "Layer" + int2str(i) + ".gen");
        comm.removeSmallComms(std::stoi(config["Save_CommunitySizeThres"]));
        comm.save(config["RESULT_DIR"] + "Layer" + int2str(i) + ".gen");
      }
      systemCall("mv LC_output.txt " + config["RESULT_DIR"] + "DSvalue" + int2str(i));
    }
    systemCall("mkdir -p " + config["RESULT_DIR"] + "Iteration/0/");//new add
    systemCall("cp " + config["RESULT_DIR"] + "Layer*.gen " + config["RESULT_DIR"] + "Iteration/0/");
    systemCall("cp " + config["RESULT_DIR"] + "Layer*.graph " + config["RESULT_DIR"] + "Iteration/0/");
    systemCall("cat " + config["RESULT_DIR"] + "Iteration/0/Layer?.gen > " + config["RESULT_DIR"] + "Iteration/0/AllLayer.gen");
    systemCall("cp " + config["RESULT_DIR"] + "DSvalue* " + config["RESULT_DIR"] + "Iteration/0/");
#ifdef DEBUG
    std::cout << "come 2" << endl;
#endif
    iterationDetail.push_back(analyze("Initialize Iteration"));
#ifdef DEBUG
    std::cout << "come 3" << endl;
#endif
    //for save maxoriginal and layer
    if (tmpmodOriG > maxmodOriG)
    {
      maxmodOriG = tmpmodOriG;
      for (int nk = 1; nk <= nLayer; nk++)
      {
        systemCall("cp " + config["RESULT_DIR"] + "*.gen " + config["RESULT_DIR"] + "maxOriginal");
        systemCall("cp " + config["RESULT_DIR"] + "*.graph " + config["RESULT_DIR"] + "maxOriginal");
      }
    }
    if (tmpmodLayerG > maxmodLayerG)
    {
      maxmodLayerG = tmpmodLayerG;
      for (int nk = 1; nk <= nLayer; nk++)
      {
        systemCall("cp " + config["RESULT_DIR"] + "*.gen " + config["RESULT_DIR"] + "maxLayer");
        systemCall("cp " + config["RESULT_DIR"] + "*.graph " + config["RESULT_DIR"] + "maxLayer");
      }
    }
#ifdef DEBUG
    std::cout << "come1" << endl;
#endif

    //Refinement stage
    for (int r = 1; r <= std::stoi(config["Number_Of_Iteration"]); r++)
    {
      fprintf(stderr, "####### Running Iteration %d with %d Layers #########\n", r, nLayer);
      for (int i = 1; i <= nLayer; i++)
        systemCall(("cp " + config["RESULT_DIR"] + "Layer" + int2str(i) + ".gen  "
          + config["RESULT_DIR"] + "Layer" + int2str(i) + ".gen.last").c_str());
      std::vector<int> order;
      order.clear();
      for (int i = 1; i <= nLayer; i++) order.push_back(i);
      //random_shuffle(order.begin(),order.end());

      for (int z = 0; z < order.size(); z++)
      {
        int i = order[z];
        std::vector<std::string> otherCommFiles;
        for (int j = 0; j < order.size(); j++)
        {
          if (i != order[j])
          {
            otherCommFiles.push_back(config["RESULT_DIR"] + "Layer" + int2str(order[j]) + ".gen");
          }
        }
        std::cout << "gen LG start" << std::endl;
        generateLayerGraph(config["DATA_DIR"] + "graph", otherCommFiles, config["RESULT_DIR"] + "Layer" + int2str(i) + ".graph", i == 1);
        std::cout << "gen LG end" << std::endl;
        generateCommunity(config["RESULT_DIR"] + "Layer" + int2str(i) + ".graph"
          , config["RESULT_DIR"] + "Layer" + int2str(i) + ".gen",
          i, Truth_NC);
        std::cout << "gen C end" << std::endl;
        //getchar();
        if (config.find("Save_CommunitySizeThres") != config.end())
        {
          Community comm;
          comm.load(config["RESULT_DIR"] + "Layer" + int2str(i) + ".gen");
          comm.removeSmallComms(std::stoi(config["Save_CommunitySizeThres"]));
          comm.save(config["RESULT_DIR"] + "Layer" + int2str(i) + ".gen");
        }
        systemCall("mv LC_output.txt " + config["RESULT_DIR"] + "DSvalue" + int2str(i));
        std::cout << "come end" << std::endl;
      }
#ifdef DEBUG
      std::cout << "come 5" << endl;
#endif
      //autoShrinkLayer();//自动确定层数的策略，已经废弃。
      iterationDetail.push_back(analyze("Iteration " + int2str(r)));
#ifdef DEBUG
      std::cout << "come 6" << endl;
#endif
      //for save max original and layer
      if (tmpmodOriG > maxmodOriG)
      {
        maxOriIter = r;
        maxmodOriG = tmpmodOriG;
        for (int nk = 1; nk <= nLayer; nk++)
        {
          systemCall("cp " + config["RESULT_DIR"] + "*.gen " + config["RESULT_DIR"] + "maxOriginal");
          systemCall("cp " + config["RESULT_DIR"] + "*.graph " + config["RESULT_DIR"] + "maxOriginal");
        }
      }
      if (tmpmodLayerG > maxmodLayerG)
      {
        maxLayIter = r;
        maxmodLayerG = tmpmodLayerG;
        for (int nk = 1; nk <= nLayer; nk++)
        {
          systemCall("cp " + config["RESULT_DIR"] + "*.gen " + config["RESULT_DIR"] + "maxLayer");
          systemCall("cp " + config["RESULT_DIR"] + "*.graph " + config["RESULT_DIR"] + "maxLayer");
        }
      }
      systemCall("mkdir -p " + config["RESULT_DIR"] + "Iteration/" + int2str(r) + "/");//new add
      systemCall("cp " + config["RESULT_DIR"] + "Layer*.gen " + config["RESULT_DIR"] + "Iteration/" + int2str(r) + "/");
      systemCall("cp " + config["RESULT_DIR"] + "Layer*.graph " + config["RESULT_DIR"] + "Iteration/" + int2str(r) + "/");
      systemCall("cat " + config["RESULT_DIR"] + "Iteration/" + int2str(r) + "/Layer?.gen > " + config["RESULT_DIR"] + "Iteration/" + int2str(r) + "/AllLayer.gen");
      systemCall("cp " + config["RESULT_DIR"] + "DSvalue* " + config["RESULT_DIR"] + "Iteration/" + int2str(r) + "/");
      //getchar();
      }
    PrintTruthAndLayer();
#ifdef DEBUG
    std::cout << "here2" << endl;
#endif
    PrintCurAndPre();
#ifdef DEBUG
    std::cout << "here3" << endl;
#endif
    showAnalyze();
#ifdef DEBUG
    std::cout << "here4" << endl;
#endif

    systemCall("mkdir " + config["RESULT_DIR"] + "Framework_" + config["Framework"]);
    systemCall("mv " + config["RESULT_DIR"] + "analyze " + config["RESULT_DIR"] + "Framework_" + config["Framework"] + "/");
    systemCall("mv " + config["RESULT_DIR"] + "analyze.pdf " + config["RESULT_DIR"] + "Framework_" + config["Framework"] + "/");
    systemCall("mv " + config["RESULT_DIR"] + "*.gen " + config["RESULT_DIR"] + "Framework_" + config["Framework"] + "/");
    systemCall("mv layer?_vs_pre.pdf " + config["RESULT_DIR"] + "Framework_" + config["Framework"] + "/");
    systemCall("mv Avg*.pdf " + config["RESULT_DIR"] + "Framework_" + config["Framework"] + "/");
    systemCall("mv Mod_Layer?_in_*.pdf " + config["RESULT_DIR"] + "Framework_" + config["Framework"] + "/");
    systemCall("mv NMI_Layer?_vs_inPreIter.pdf " + config["RESULT_DIR"] + "Framework_" + config["Framework"] + "/");
    systemCall("mv truth?_vs_layer?.pdf " + config["RESULT_DIR"] + "Framework_" + config["Framework"] + "/");
    systemCall("mv Truth_vs_layer " + config["RESULT_DIR"] + "Framework_" + config["Framework"] + "/");
    systemCall("mv graphicData " + config["RESULT_DIR"] + "Framework_" + config["Framework"] + "/");
    systemCall("mv Cur_vs_Pre  " + config["RESULT_DIR"] + "Framework_" + config["Framework"] + "/");
    systemCall("cat " + config["RESULT_DIR"] + "maxLayer/Layer?.gen > " + config["RESULT_DIR"] + "maxLayer/AllLayer.gen");
    systemCall("cat " + config["RESULT_DIR"] + "maxOriginal/Layer?.gen > " + config["RESULT_DIR"] + "maxOriginal/AllLayer.gen");
    //copy graph to Framwork directory
    systemCall("mv " + config["RESULT_DIR"] + "*.graph " + config["RESULT_DIR"] + "Framework_" + config["Framework"] + "/");
    systemCall("rm " + config["RESULT_DIR"] + "*.gen.last");

    //systemCall("rm layer?_vs_pre.pdf Mod_Layer?_in_*.pdf NMI_Layer?_vs_inPreIter.pdf truth?_vs_layer?.pdf");
    //systemCall("rm -f "+config["RESULT_DIR"]+"*");

    //end = time(NULL);
    //outfile2<<difftime(end, start)<< std::endl;
    gettimeofday(&end, NULL);
    span = end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec) / 1000000.0;
    outfile2 << span << std::endl;
    }
  outfile2.close();
  //systemCall("cp "+config["RESULT_DIR"]+"Frameworktime.txt "+config["RESULT_DIR"]+"Time/");

  //systemCall("rm -f -r "+config["TMP_DIR"]);
  systemCall("mkdir " + config["RESULT_DIR"] + "Data");
  systemCall("cp " + config["DATA_DIR"] + "graph " + config["RESULT_DIR"] + "Data/");
  for (int i = 0; i < truth.size(); i++)
  {
    systemCall("cp " + config["DATA_DIR"] + truth[i] + " " + config["RESULT_DIR"] + "Data/");
  }

  puts("Finished");

  return 0;
    }
