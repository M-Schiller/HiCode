#pragma once

#ifndef SINGLELAYER_METHOD_H
#define SINGLELAYER_METHOD_H

class SingleLayer_Method
{
public:
  virtual ~SingleLayer_Method() = default;
  std::map<std::string, std::string> config{};
  void setConfig(std::map<std::string, std::string> config)
  {
    this->config = config;
  }
  virtual bool checkRequire()
  {
    puts("ERR");
    return true;
  }
  virtual void generateCommunity(const std::string& communityFile, int truth_NC)
  {
    puts("ERR");
  }
};
#endif // SINGLELAYER_METHOD_H
