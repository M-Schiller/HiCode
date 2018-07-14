#pragma once

#ifndef METRIC_OVERLAPPING_NMI_H
#define METRIC_OVERLAPPING_NMI_H

inline int intersect_size(std::vector<int> Xi, std::vector<int> Yj)
{
  int Xi_size = Xi.size();
  int Yj_size = Yj.size();
  if (Xi_size == 0 || Yj_size == 0)
    return 0;
  std::vector<int> v(std::max(Xi_size, Yj_size));
  std::vector<int>::iterator it;

  int Xi_array[Xi_size];
  int Yj_array[Yj_size];
  int i, j;
  for (i = 0; i < Xi_size; i++)
    Xi_array[i] = Xi[i];
  for (j = 0; j < Yj_size; j++)
    Yj_array[i] = Yj[i];
  int* pXi = Xi_array;
  int* pYj = Yj_array;
  std::sort(pXi, pXi + Xi_size);
  std::sort(pYj, pYj + Yj_size);
  it = std::set_intersection(pXi, pXi + Xi_size, pYj, pYj + Yj_size, v.begin());
  v.resize(it - v.begin());
  return v.size();
}

inline double Entropy_h(double p)
/*
    require: p>=0
    return:
         { -plog(p),    if p>0
    h(p)=|
         {  0      ,    otherwise
*/
{
  return p > 0 ? (p*log(p)) : 0;
}

inline double Entropy_H(const std::vector<int>& Xi, int X_size)
/*
    calculate H(Xi)
    require: X_size>|Xi|>0
    return:  H(Xi) = h[P(Xi∈X)] + h[P(Xi∉X)]
*/

{
  double PXi_1 = (double)Xi.size() / (double)X_size; //P(Xi∈X)
  std::cout << "Xi.size = " << Xi.size() << std::endl;
  double PXi_0 = 1 - PXi_1;        //P(Xi∉X)
  std::cout << "PXi_0 = " << PXi_0 << std::endl;
  std::cout << "PXi_1 = " << PXi_1 << std::endl;
  return Entropy_h(PXi_1) + Entropy_h(PXi_0);
}

inline double each_right_Conditional_Entropy(
  const std::vector<int>& Xi,
  const std::vector<int>& Yj,
  int X_size,
  int Y_size)
  /*
      Xi 就是 X_left_comm_vector， 同样的Yj 就是 Y_right_comm_vector
      calculate H(Xi|Yj)
      require:  X_size > |Xi| > 0
                Y_size > |Yj| > 0
      return:
                           {   H(Xi,Yj) - H(Yj), if constraint satisfied
                H(Xi|Yj) = |
                           {   -9999(-INF)     , otherwise
  */
{
  //计算H(Xi,Yj) = h[P(1,1)] + h[P(1,0)] + h[P(0,1)] + h[P(0,0)]
  int intersectionXi_Yj = intersect_size(Xi, Yj);
  int unionXi_Yj = Xi.size() - intersectionXi_Yj + Yj.size();
  double P_1_1 = intersect_size(Xi, Yj) / (double)X_size;//calculate P(1,1) = |Xi ∩ Yj|/n
  double P_1_0 = (Xi.size() - intersectionXi_Yj) / (double)X_size;// P(1,0) = (|Xi|-|Xi ∩ Yj|)/n
  double P_0_1 = (Yj.size() - intersectionXi_Yj) / (double)X_size;// P(0,1) = (|Yi|-|Xi ∩ Yj|)/n
  double P_0_0 = 1 - (unionXi_Yj) / (double)X_size;             // P(0,0) = |Xi ∪ Yj|/n
  double h11 = Entropy_h(P_1_1);//h[P(1,1)]
  double h10 = Entropy_h(P_1_0);
  double h01 = Entropy_h(P_0_1);
  double h00 = Entropy_h(P_0_0);

  if ((h11 + h10) > (h01 + h00)) //(satisfaction check)
      //判断是否满足限制(constraint)h[P(1,1)] + h[P(1,1)] > h[P(0,1)] + h[P(1,0)]
    return h11 + h10 + h01 + h00 - Entropy_H(Yj, Y_size);
  else
    //表示不满足该限制条件，设置为-9999就可以在each_left_Conditional_Entropy中判断是不是每个j都不满足该限制条件
    return -9999;

  //因为只有都不满足，each_left_Conditional_Entropy中的MinH_XiYj才会是-9999
}

inline double each_left_Conditional_Entropy(
  const std::vector<int>& X_left_comm_vector,
  Community right_commsY,
  int X_size)
  /*
      Y = right_commsY
      calculate H(Xi|Y)
      require: X_size > |Xi| > 0
           and right_commsY  ≠ ∅
      return:
                         {    min_{Yj∈Y}H(Xi|Yj), if min_{Yi∈Y}H(Xi|Yj) ≠ -INF(-9999)
               H(Xi|Y) = |
                         {    H(Xi)             , otherwise
  */
{
  int k = right_commsY.NC;
  std::vector<int> Y_right_comm_vector = right_commsY.id[1];
  double H_XiYi;
  double MinH_XiYj = each_right_Conditional_Entropy(X_left_comm_vector,
    Y_right_comm_vector,
    X_size,
    right_commsY.N);
  for (int j = 2; j <= k; j++)
  {
    Y_right_comm_vector = right_commsY.id[j];
    double H_XiYj = each_right_Conditional_Entropy(X_left_comm_vector,
      Y_right_comm_vector,
      X_size,
      right_commsY.N);

    if (MinH_XiYj > H_XiYj)
      MinH_XiYj = H_XiYj;
  }
  if (MinH_XiYj == -9999)//表示对每个j，都不满足限制h[P(1,1)]+h[P(0,0)]>h[P(0,1)]+h[P(1,0)]
    MinH_XiYj = Entropy_H(X_left_comm_vector, X_size);//在这种情况下，赋值H(Xi)

  return MinH_XiYj;
}

inline double Conditional_Entropy(Community left_comms, const Community& right_comms)
/*
    calculate Conditional Entropy H(X|Y)
    require:    left_comms ≠ ∅ and right_comms ≠ ∅
    return:               1   _k_   H(Xi|Y)
                H(X|Y) = ---  \  ` ---------
                          k   /__,   H(Xi)
                              i=1
*/
{
  int k = left_comms.NC;
  int X_size = left_comms.N;
  double sum = 0.0;
  for (int i = 1; i <= k; i++)
  {
    std::vector<int> X_left_comm_vector = left_comms.id[i];
    sum += each_left_Conditional_Entropy(X_left_comm_vector, right_comms, X_size);
  }
  return sum / (double)k;
}

inline double calcNMI_real(Community commsX, Community commsY)
/*
    calculate Overlapping NMI
    require:?????????? 怎样的两个社团可以拿来比？需要做什么预处理。
            目前需要
            不管原来社团集合中是否允许存在不在任何社团中的点，传入的两个社团必须做到每个点至少在一个社团中，
            并且每个社团的size就是处理后的size。
    return:
            Overlapping NMI = 1 - 0.5*[H(X|Y) + H(Y|X)]
*/
{
  if (commsX.N <= 0
    || commsY.N <= 0)
  {
    std::cout << "Info: one of the two community sets is empty" << std::endl;
    return 0;
  }
  else
  {
    double H_X_in_Y = Conditional_Entropy(commsX, commsY);
    double H_Y_in_X = Conditional_Entropy(commsY, commsX);
    return (1 - 0.5*(H_X_in_Y + H_Y_in_X));// 1 - 0.5*[H(X|Y) + H(Y|X)]}
  }
}

inline double calcNMI(const std::string& commsX_file, const std::string& commsY_file)
//对用户提供的计算overlapping NMI的接口
{
  std::cout << commsX_file << ' ' << commsY_file << std::endl;
  Community X;
  X.load(commsX_file);
  Community Y;
  Y.load(commsY_file);
  if (config.find("Metric_CommunitySizeThres") != config.end())
  {
    X.removeSmallComms(std::stoi(config["Metric_CommunitySizeThres"]));
    Y.removeSmallComms(std::stoi(config["Metric_CommunitySizeThres"]));
  }
  double score = calcNMI_real(X, Y);//调用真正计算Overlapping NMI的接口,返回OverlappingNMI值
  std::cout << commsX_file << " vs " << commsY_file << " Overlapping NMI " << score << std::endl;
  return score;//返回OverlappingNMI值给用户
}

inline double calcEntropy(Community comm)//commnities set's Entropy
{
  double E = 0;
  for (int i = 1; i <= comm.NC; i++)
  {
    double p = comm.id[i].size() / (double)(comm.N);
    if (p > 0) E -= p * log(p);
  }
  return E;
}
#endif // METRIC_OVERLAPPING_NMI_H
