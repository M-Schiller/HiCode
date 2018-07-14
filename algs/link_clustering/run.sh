#! /usr/bin/env bash

g++ -o ./algs/LinkCommunity/dataTrans ./algs/LinkCommunity/dataTrans.cpp -std=c++1z

#./algs/LinkCommunity/dataTrans ./data/grad/graph

#copy ./algs/LC/graphNew ./data/grad/graphNew
g++ -o ./algs/LinkCommunity/calcAndWrite_Jaccards ./algs/LinkCommunity/calcAndWrite_Jaccards.cpp -std=c++1z

#./algs/LinkCommunity/calcAndWrite_Jaccards Graph.pairs Graph.jaccs

g++ -o ./algs/LinkCommunity/callCluster ./algs/LinkCommunity/callCluster.cpp -std=c++1z


#sh ./algs/LinkCommunity/loop_thresholds.sh
