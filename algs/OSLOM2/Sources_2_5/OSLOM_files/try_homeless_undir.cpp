#include <map>
#include <set>

int oslom_net_global::try_to_assign_homeless(module_collection & Mcoll, bool anyway)
{
  Mcoll.put_gaps();

  //if(paras.print_cbs)
    //cout<<"checking homeless nodes "<< std::endl;

  std::deque<int> homel;
  Mcoll.homeless(homel);

  int before_procedure = homel.size();

  if (homel.empty())
    return before_procedure;

  /*cout<<"homel"<< std::endl;
  print_id(homel, std::cout);*/

  std::set<int> called;						// modules connected to homeless nodes
  std::map<int, std::set<int>> homel_module;		// maps the homeless node with the modules it's connected to

  for (int i : homel)
  {
    std::set<int> thish;
    for (int j = 0; j < vertices[homel[i]]->links->size(); j++)
    {
      int & neigh = vertices[homel[i]]->links->l[j];

      for (auto itk : Mcoll.memberships[neigh])
      {
        called.insert(itk);
        thish.insert(itk);
      }
    }

    if (!thish.empty())
      homel_module[i] = thish;
  }

  std::map<int, int> module_kin;
  std::map<int, int> module_ktot;

  for (auto its : called)
  {
    module_kin[its] = cast_int(kin_m(Mcoll.modules[its]));
    module_ktot[its] = cast_int(ktot_m(Mcoll.modules[its]));
  }

  std::map<int, std::deque<int>> to_check;			// module - homeless nodes added to that
  for (auto& itm : homel_module)
  {
    double cmin = 1.1;
    int belongs_to = -1;

    //cout<<"homeless node: "<<id_of(itm->first)<< std::endl;

    for (auto its = itm.second.begin(); its != itm.second.end(); ++its)
    {
      int kin_node = cast_int(vertices[itm.first]->kplus_m(Mcoll.modules[*its]));

      /*cout<<"module: "<<*its<<" kin: "<<module_kin[*its]<<"  ktot: "<<module_ktot[*its]<<" kin h "<<kin_node<< std::endl;
      print_ri(Mcoll.modules[*its]);*/

      int kout_g = module_ktot[*its] - module_kin[*its];
      int tm = oneM - module_ktot[*its];
      //double rh= compute_r_hyper(kin_node, kout_g, tm, vertices[itm->first]->stub_number);
      double kinw = vertices[itm.first]->kplus_w(Mcoll.modules[*its]);
      //double weight_part= log_together(kinw, kin_node);

      double rh = compute_global_fitness_randomized_short(kin_node, kout_g, tm, vertices[itm.first]->stub_number, kinw);

      //double cs=  1 - pow(1 - rh, dim - Mcoll.modules[*its].size());
      //cout<<"rh: "<<rh<<" ..."<< std::endl;
      if (rh < cmin)
      {
        cmin = rh;
        belongs_to = *its;
      }
    }

    if (belongs_to != -1)
    {
      if (to_check.find(belongs_to) == to_check.end())
      {
        std::deque<int> void_d;
        to_check[belongs_to] = void_d;
      }
      to_check[belongs_to].push_back(itm.first);
    }
    //if(paras.print_cbs)
      //cout<<"homeless node: "<<id_of(itm->first)<<" belongs_to "<<belongs_to<<" cmin... "<<cmin<< std::endl;
      //cherr();
  }

  //if(paras.print_cbs)
    //cout<<"homeless node: "<<homel.size()<<" try_to_assign: "<<homel_module.size()<<" modules to check: "<<to_check.size()<< std::endl;

  // **** try the groups with the homeless //******************
  bool something = false;

  for (auto& itm : to_check)
  {
    std::deque<int> union_deque = Mcoll.modules[itm.first];

    for (unsigned i = 0; i < itm.second.size(); i++)
      union_deque.push_back(itm.second[i]);

    if (anyway)
    {
      something = true;
      Mcoll.insert(union_deque, ran4() + paras.threshold);
    }
    else
    {
      std::deque<int> grbe;
      double bs = CUP_check(union_deque, grbe);

      //cout<<"union_deque after "<<itm->first<<" size: "<<grbe.size()<< std::endl;
      if (grbe.size() > 1)
      {
        something = true;
        Mcoll.insert(grbe, bs);
      }
    }
  }

  if (something)
  {
    Mcoll.compute_inclusions();
  }

  return before_procedure;
}
