#pragma once

#ifndef HISTOGRAMS_INCLUDED
#define HISTOGRAMS_INCLUDED

#include <deque>
#include <algorithm>

int intlog_binning(std::deque<int> c, int number_of_bins, std::deque<double> & Xs, std::deque<double> & Ys, std::deque<double> & var)
{
  // this function is to make a log_histogram along the x axis and to compute y-errors
  //prints(c);

  Xs.clear();
  Ys.clear();
  var.clear();

  std::deque<int> d;
  for (int i : c)
  {
    if (i > 0)
    {
      d.push_back(i);
    }
  }
  c.clear();
  c = d;
  std::sort(c.begin(), c.end());

  int max = c[c.size() - 1];
  int min = c[0];

  std::deque <double> bins;
  double step = log(min);
  if (max == min)
  {
    max++;
  }

  const double bin_width = (log(max) - log(min)) / number_of_bins;

  while (step <= log(max) + 4 * bin_width)
  {
    //cout<<"step: "<<(exp(step))<< std::endl;

    bins.push_back(exp(step));
    step += bin_width;
  }

  //cout<<"bins"<< std::endl;
  //prints(bins);

  std::deque<int> hist;
  std::deque<double> hist2;

  int index = 0;

  for (int i : c)
  {
    while (i - bins[index] > -1e-6)
    {
      index++;
      hist.push_back(0);
      hist2.push_back(0);
    }
    //cout<<index<<" "<<c[i]<<" "<<bins[index]<<" "<< std::endl;
    hist[index - 1]++;
    hist2[index - 1] += double(i);
  }

  std::deque<int> integers;
  index = 0;
  for (int i = min; i < bins[bins.size() - 1] - 1; i++)
  {
    while (i - bins[index] > -1e-6)
    {
      index++;
      integers.push_back(0);
    }

    //cout<<i<<" "<<index<<" "<<integers[index-1]<<" ********"<< std::endl;
    integers[index - 1]++;
  }

  for (int i = 0; i < hist.size(); i++)
  {
    double h1 = bins[i];
    double h2 = bins[i + 1];

    if (hist[i] > 0)
    {
      Xs.push_back(hist2[i] / hist[i]);
      double y = double(hist[i]) / (c.size()*integers[i]);
      Ys.push_back(y);

      //cout<<h1<<" "<<h2<<" "<<y<<" "<<(hist2[i]/hist[i])<<" ***"<< std::endl;

      var.push_back(double(hist[i]) / (c.size()*c.size()*integers[i] * integers[i]));
      //cout<<"-> "<<hist[i]<<" "<<double(hist[i])/(c.size()*c.size()*integers[i]*integers[i])<< std::endl;
    }
  }
  return 0;
}

template <typename T>
int xybinning(
  std::deque<T> &c,
  std::deque<T> &d,
  int number_of_bins,
  std::deque<double> & xs,
  std::deque<double> & ys,
  std::deque<double> & var,
  std::deque<int> & nums)
{
  // so, this function takes two datasets (c and d) and gathers the data in bin, takes xs and ys as the average in each bin, var is the variance of the y average
  // the difference with the same stuff called not_norm_histogram is that the other one averages x with y weights.

  xs.clear();
  ys.clear();
  var.clear();
  nums.clear();

  double min = double(c[0]);
  double max = double(c[0]);

  for (int i = 0; i < c.size(); i++)
  {
    if (min > double(c[i]))
    {
      min = double(c[i]);
    }
    if (max<double(c[i]))
    {
      max = double(c[i]);
    }
  }

  min -= 1e-6;
  max += 1e-6;

  if (max == min)
  {
    max += 1e-3;
  }

  std::deque<std::deque<double>> hist_x;		// x values in the bin
  std::deque<std::deque<double>>  hist_y;		// y values in the bin

  double step = min;
  double bin = (max - min) / number_of_bins;		// bin width

  std::deque<double> f;
  while (step <= max + 2 * bin)
  {
    hist_x.push_back(f);
    hist_y.push_back(f);
    step += bin;
  }

  for (int i = 0; i < c.size(); i++)
  {
    double data = double(c[i]);

    if (data >= min && data <= max)
    {
      int index = int((data - min) / bin);
      //cout<<data<<" "<<exp(data)<<" "<<index<< std::endl;

      hist_x[index].push_back(double(c[i]));
      hist_y[index].push_back(double(d[i]));
    }
  }

  for (int i = 0; i < hist_x.size() - 1; i++)
  {
    double x = average_func(hist_x[i]);
    double y = average_func(hist_y[i]);

    //cout<<x<<" "<<exp(x)<<" "<<y<< std::endl;

    if (!hist_y[i].empty())
    {
      xs.push_back(x);
      ys.push_back(y);
      var.push_back(variance_func(hist_y[i]) / double(hist_y[i].size()));
      nums.push_back(hist_y[i].size());
      //cout<<x<<" "<<exp(x)<<" "<<y<<" "<<(hist_y[i].size())<<" "<<variance_func(hist_y[i])<< std::endl;
    }
  }

  for (int i = 0; i < var.size(); i++)
  {
    if (var[i] < 1e-8)
    {
      var[i] = 1e-8;
    }
  }

  return 0;
}

template <typename T>
int xybinning(std::deque<T> &c, std::deque<T> &d, int number_of_bins, std::deque<double> & xs, std::deque<double> & ys, std::deque<double> & var)
{
  std::deque<int> nums;
  return xybinning(c, d, number_of_bins, xs, ys, var, nums);
}

double compute_quantiles(double q, std::deque<double> & y, std::deque<double> & qs)
{
  int qv = cast_int((1 - q) / 2 * y.size());
  qv = std::max(qv, 0);

  if (qv >= y.size())
  {
    qv = y.size() - 1;
  }
  qs.push_back(y[qv]);

  qv = cast_int((1 + q) / 2 * y.size());
  qv = std::max(qv, 0);

  if (qv >= y.size())
  {
    qv = y.size() - 1;
  }
  qs.push_back(y[qv]);
  return 0;
}

template <typename T>
int xybinning_quantiles(
  std::deque<T> &c,
  std::deque<T> &d,
  int number_of_bins,
  std::deque<double> & xs,
  std::deque<double> & ys,
  std::deque<double> & var,
  std::deque<int> & nums,
  std::deque<std::deque<double>> & Mq,
  double qa,
  double qb)
{
  // so, this function takes two datasets (c and d) and gathers the data in bin, takes xs and ys as the average in each bin, var is the variance of the y average
  // the difference with the same stuff called not_norm_histogram is that the other one averages x with y weights.

  xs.clear();
  ys.clear();
  var.clear();
  nums.clear();
  Mq.clear();

  double min = double(c[0]);
  double max = double(c[0]);

  for (int i = 0; i < c.size(); i++)
  {
    if (min > double(c[i]))
    {
      min = double(c[i]);
    }
    if (max<double(c[i]))
    {
      max = double(c[i]);
    }
  }

  min -= 1e-6;
  max += 1e-6;

  if (max == min)
  {
    max += 1e-3;
  }

  std::deque<std::deque<double>> hist_x;		// x values in the bin
  std::deque<std::deque<double>>  hist_y;		// y values in the bin

  double step = min;
  const double bin_width = (max - min) / number_of_bins;

  std::deque<double> f;
  while (step <= max + 2 * bin_width) {
    hist_x.push_back(f);
    hist_y.push_back(f);
    step += bin_width;
  }

  for (int i = 0; i < c.size(); i++)
  {
    double data = double(c[i]);

    if (data >= min && data <= max)
    {
      int index = int((data - min) / bin_width);
      //cout<<data<<" "<<exp(data)<<" "<<index<< std::endl;

      hist_x[index].push_back(double(c[i]));
      hist_y[index].push_back(double(d[i]));
    }
  }

  for (int i = 0; i < hist_x.size() - 1; i++)
  {
    double x = average_func(hist_x[i]);
    double y = average_func(hist_y[i]);

    //cout<<x<<" "<<exp(x)<<" "<<y<< std::endl;

    if (!hist_y[i].empty())
    {
      xs.push_back(x);
      ys.push_back(y);
      var.push_back(variance_func(hist_y[i]) / double(hist_y[i].size()));
      nums.push_back(hist_y[i].size());
      sort(hist_y[i].begin(), hist_y[i].end());

      std::deque<double> qs;
      compute_quantiles(qa, hist_y[i], qs);
      compute_quantiles(qb, hist_y[i], qs);

      Mq.push_back(qs);
      //cout<<x<<" "<<exp(x)<<" "<<y<<" "<<(hist_y[i].size())<<" "<<variance_func(hist_y[i])<< std::endl;
    }
  }

  for (int i = 0; i < var.size(); i++)
  {
    if (var[i] < 1e-8)
    {
      var[i] = 1e-8;
    }
  }
  return 0;
}

template <typename T>
int log_histogram(std::deque<T> &c, std::ostream & out, int number_of_bins)		// c is the set od data, min is the lower bound, max is the upper one
{
  std::deque <T> d;
  for (int i = 0; i < c.size(); i++)
  {
    if (c[i] > 0)
    {
      d.push_back(c[i]);
    }
  }

  c.clear();
  c = d;

  double min = double(c[0]);
  double max = double(c[0]);

  for (int i = 0; i < c.size(); i++)
  {
    if (min > double(c[i]))
    {
      min = double(c[i]);
    }

    if (max<double(c[i]))
    {
      max = double(c[i]);
    }
  }

  std::deque <int> hist;
  std::deque <double> hist2;
  std::deque <double> bins;
  double step = log(min);
  if (max == min)
    max++;

  double bin_width = (log(max) - log(min)) / number_of_bins;		// bin width

  while (step <= log(max) + 2 * bin_width) {
    bins.push_back(exp(step));
    hist.push_back(0);
    hist2.push_back(0);
    step += bin_width;
  }

  for (int i = 0; i < c.size(); i++) {
    int index = bins.size() - 1;
    for (int j = 0; j < bins.size() - 1; j++) if ((fabs(double(c[i]) - bins[j]) < 1e-7) || (double(c[i]) > bins[j] && double(c[i]) < bins[j + 1])) {
      // this could be done in a more efficient way

      index = j;
      break;
    }

    //cout<<hist[index]<<" "<<index<< std::endl;

    hist[index]++;
    hist2[index] += double(c[i]);
  }

  for (int i = 0; i < hist.size() - 1; i++)
  {
    double h1 = bins[i];
    double h2 = bins[i + 1];

    double x = hist2[i] / hist[i];
    double y = double(hist[i]) / (c.size()*(h2 - h1));

    if (fabs(y) > 1e-10)
    {
      out << x << "\t" << y << std::endl;
    }
  }
  return 0;
}

int log_histogram(std::deque<double> &c, std::deque<double> &c2, std::ostream & out, int number_of_bins)		// c is the set od data, min is the lower bound, max is the upper one
{
  // c must be sorted and must be c[i]>0

  double min = double(c[0]);
  auto max = double(c[0]);

  for (int i = 0; i < c.size(); i++) {
    if (min > double(c[i]))
    {
      min = double(c[i]);
    }

    if (max<double(c[i]))
    {
      max = double(c[i]);
    }
  }

  std::deque <int> hist0;
  std::deque <double> hist;
  std::deque <double> hist2;
  std::deque <double> bins;
  double step = log(min);
  if (max == min)
    max++;

  const double bin_width = (log(max) - log(min)) / number_of_bins;

  while (step <= log(max) + 2 * bin_width)
  {
    bins.push_back(exp(step));
    hist0.push_back(0);
    hist.push_back(0);
    hist2.push_back(0);
    step += bin_width;
  }

  for (int i = 0; i < c.size(); i++)
  {
    int index = bins.size() - 1;
    for (int j = 0; j < bins.size() - 1; j++)
    {
      if (abs(double(c[i]) - bins[j]) < 1e-7
        || (double(c[i]) > bins[j] && double(c[i]) < bins[j + 1]))
      {
        // this could be done in a more efficient way

        index = j;
        break;
      }
    }

    //cout<<hist[index]<<" "<<index<< std::endl;
    hist0[index]++;
    hist[index] += c2[i];
    hist2[index] += double(c[i])*c2[i];
  }

  for (int i = 0; i < hist.size() - 1; i++)
  {
    double h1 = bins[i];
    double h2 = bins[i + 1];

    double x = hist2[i] / hist[i];
    double y = double(hist[i]) / (hist0[i]);

    if (abs(y) > 1e-10)
    {
      out << x << "\t" << y << std::endl;
    }
  }
  return 0;
}

template <typename T>
int histogram(std::vector <T> &c, std::ostream & out, int number_of_bins, double b1, double b2)
{
  // this should be OK
  // c is the set of data, b1 is the lower bound, b2 is the upper one (if they are equal, default limits are used)

  double min = double(c[0]);
  double max = double(c[0]);

  for (int i = 0; i < c.size(); i++)
  {
    if (min > double(c[i]))
    {
      min = double(c[i]);
    }

    if (max<double(c[i]))
    {
      max = double(c[i]);
    }
  }

  min -= 1e-6;
  max += 1e-6;

  if (b1 != b2)
  {
    min = b1;
    max = b2;
  }

  if (max == min)
  {
    max += 1e-3;
  }

  std::deque <int> hist;
  std::deque <double> hist2;

  double step = min;
  double bin_width = (max - min) / number_of_bins;

  while (step <= max + 2 * bin_width)
  {
    hist.push_back(0);
    hist2.push_back(0);
    step += bin_width;
  }

  for (int i = 0; i < c.size(); i++)
  {
    double data = double(c[i]);

    if (data > min && data <= max)
    {
      int index = int((data - min) / bin_width);

      hist[index]++;
      hist2[index] += double(c[i]);
    }
  }

  for (int i = 0; i < hist.size() - 1; i++)
  {
    double x = hist2[i] / hist[i];
    double y = double(hist[i]) / (c.size()*bin_width);

    if (fabs(y) > 1e-10)
      out << x << "\t" << y << std::endl;
  }
  return 0;
}

template <typename T>
int not_norm_histogram_correlated(std::deque<T> &c, std::deque<T> &d, std::ostream & out, int number_of_bins, double b1, double b2)
{
  // c is the x axis, d the y, b1 is the lower bound, b2 is the upper one (if they are equal, default limits are used)
  double min = double(c[0]);
  double max = double(c[0]);

  for (int i = 0; i < c.size(); i++)
  {
    if (min > double(c[i]))
    {
      min = double(c[i]);
    }

    if (max<double(c[i]))
    {
      max = double(c[i]);
    }
  }

  min -= 1e-6;
  max += 1e-6;

  if (b1 != b2)
  {
    min = b1;
    max = b2;
  }

  if (max == min)
  {
    max += 1e-3;
  }

  std::deque <int> hist;			// frequency in the bin
  std::deque <double> hist_x;		// x sum in the bin
  std::deque <double> hist_y;		// y sum in the bin

  double step = min;
  double bin = (max - min) / number_of_bins;		// bin width

  while (step <= max + 2 * bin) {
    hist.push_back(0);
    hist_x.push_back(0);
    hist_y.push_back(0);
    step += bin;
  }

  for (int i = 0; i < c.size(); i++)
  {
    double data = double(c[i]);

    if (data > min && data <= max)
    {
      int index = int((data - min) / bin);

      hist[index]++;
      hist_x[index] += double(c[i]);
      hist_y[index] += double(d[i]);
    }
  }

  for (int i = 0; i < hist.size() - 1; i++)
  {
    double x = hist_x[i] / hist[i];
    double y = hist_y[i] / hist[i];;

    if (fabs(y) > 1e-10)
    {
      out << x << "\t" << y << std::endl;
    }
  }
  return 0;
}

template <typename T>
int histogram(std::deque <T> &c, std::ostream & out, int number_of_bins, double b1, double b2)
{
  // this should be OK
  // c is the set of data, b1 is the lower bound, b2 is the upper one (if they are equal, default limits are used)

  double min = double(c[0]);
  double max = double(c[0]);

  for (int i = 0; i < c.size(); i++)
  {
    if (min > double(c[i]))
    {
      min = double(c[i]);
    }

    if (max<double(c[i]))
    {
      max = double(c[i]);
    }
  }

  min -= 1e-6;
  max += 1e-6;

  if (b1 != b2)
  {
    min = b1;
    max = b2;
  }

  if (max == min)
  {
    max += 1e-3;
  }

  std::deque <int> hist;
  std::deque <double> hist2;

  double step = min;
  double bin_width = (max - min) / number_of_bins;

  while (step <= max + 2 * bin_width)
  {
    hist.push_back(0);
    hist2.push_back(0);
    step += bin_width;
  }

  for (int i = 0; i < c.size(); i++)
  {
    double data = double(c[i]);

    if (data > min && data <= max)
    {
      int index = int((data - min) / bin_width);
      hist[index]++;
      hist2[index] += double(c[i]);
    }
  }

  for (int i = 0; i < hist.size() - 1; i++)
  {
    double x = hist2[i] / hist[i];
    double y = double(hist[i]) / (c.size() * bin_width);

    if (fabs(y) > 1e-10)
    {
      out << x << "\t" << y << std::endl;
    }
  }
  return 0;
}

//ofstream pout("hist.dat");
//not_norm_histogram(c, pout, 10, 0, 0);

template <typename T>
int not_norm_histogram(std::vector<T> &c, std::ostream & out, int number_of_bins, double b1, double b2)
{
  // this should be OK
  // c is the set of data, b1 is the lower bound, b2 is the upper one (if they are equal, default limits are used)

  double min = double(c[0]);
  double max = double(c[0]);

  for (int i = 0; i < c.size(); i++)
  {
    if (min > double(c[i]))
    {
      min = double(c[i]);
    }

    if (max<double(c[i]))
    {
      max = double(c[i]);
    }
  }

  min -= 1e-6;
  max += 1e-6;

  if (b1 != b2)
  {
    min = b1;
    max = b2;
  }

  if (max == min)
  {
    max += 1e-3;
  }

  std::deque <int> hist;
  std::deque <double> hist2;

  double step = min;
  const double bin_width = (max - min) / number_of_bins;

  while (step <= max + 2 * bin_width) {
    hist.push_back(0);
    hist2.push_back(0);
    step += bin_width;
  }

  for (int i = 0; i < c.size(); i++)
  {
    double data = double(c[i]);

    if (data > min && data <= max)
    {
      int index = int((data - min) / bin_width);

      hist[index]++;
      hist2[index] += double(c[i]);
    }
  }

  for (int i = 0; i < hist.size() - 1; i++)
  {
    double x = hist2[i] / hist[i];
    double y = double(hist[i]) / (c.size());

    if (fabs(y) > 1e-10)
    {
      out << x << "\t" << y << std::endl;
    }
  }
  return 0;
}

template <typename T>
int not_norm_histogram(std::deque<T> &c, std::ostream & out, int number_of_bins, double b1, double b2)
{
  // this should be OK
  // c is the set of data, b1 is the lower bound, b2 is the upper one (if they are equal, default limits are used)

  double min = double(c[0]);
  double max = double(c[0]);

  for (int i = 0; i < c.size(); i++)
  {
    if (min > double(c[i]))
    {
      min = double(c[i]);
    }
    if (max<double(c[i]))
    {
      max = double(c[i]);
    }
  }

  min -= 1e-6;
  max += 1e-6;

  if (b1 != b2)
  {
    min = b1;
    max = b2;
  }

  if (max == min)
  {
    max += 1e-3;
  }

  std::deque<int> hist;
  std::deque<double> hist2;

  double step = min;
  const double bin_width = (max - min) / number_of_bins;

  while (step <= max + 2 * bin_width) {
    hist.push_back(0);
    hist2.push_back(0);
    step += bin_width;
  }

  for (int i = 0; i < c.size(); i++)
  {
    double data = double(c[i]);

    if (data > min && data <= max)
    {
      int index = int((data - min) / bin_width);

      hist[index]++;
      hist2[index] += double(c[i]);
    }
  }

  for (int i = 0; i < hist.size() - 1; i++)
  {
    double x = hist2[i] / hist[i];
    double y = double(hist[i]) / (c.size());

    if (fabs(y) > 1e-10)
      out << x << "\t" << y << "\t" << sqrt(hist[i]) / c.size() << std::endl;
  }
  return 0;
}

int histogram(std::deque<double> &c, std::deque<double> &c2, std::ostream & out, int number_of_bins, double b1, double b2)
{
  // this should be OK
  // c is the set of data, b1 is the lower bound, b2 is the upper one (if they are equal, default limits are used)
  double min = c[0];
  double max = c[0];

  for (double i : c)
  {
    if (min > double(i))
    {
      min = double(i);
    }
    if (max < double(i))
    {
      max = double(i);
    }
  }

  min -= 1e-6;
  max += 1e-6;

  if (b1 != b2)
  {
    min = b1;
    max = b2;
  }

  if (max == min)
  {
    max += 1e-3;
  }

  std::deque <int> hist0;
  std::deque <double> hist;
  std::deque <double> hist2;

  double step = min;
  const double bin_width = (max - min) / number_of_bins;

  while (step <= max + 2 * bin_width)
  {
    hist0.push_back(0);
    hist.push_back(0);
    hist2.push_back(0);
    step += bin_width;
  }

  for (int i = 0; i < c.size(); i++)
  {
    double data = c[i];

    if (data > min && data <= max) {
      auto index = int((data - min) / bin_width);

      hist0[index]++;
      hist[index] += c2[i];
      hist2[index] += double(c[i] * c2[i]);
    }
  }

  for (int i = 0; i < hist.size() - 1; i++)
  {
    double x = hist2[i] / hist[i];
    auto y = double(hist[i] / hist0[i] / bin_width);

    if (fabs(y) > 1e-10)
    {
      out << x << "\t" << y << std::endl;
    }
  }
  return 0;
}

int not_norm_histogram(std::deque<double> &c, std::deque<double> &c2, std::ostream & out, int number_of_bins, double b1, double b2)
{
  // this should be OK
  // c is the set of data, b1 is the lower bound, b2 is the upper one (if they are equal, default limits are used)

  double min = c[0];
  double max = c[0];

  for (double i : c)
  {
    if (min > i)
    {
      min = c[i];
    }
    if (max < c[i])
    {
      max = c[i];
    }
  }

  min -= 1e-6;
  max += 1e-6;

  if (b1 != b2)
  {
    min = b1;
    max = b2;
  }

  if (max == min)
  {
    max += 1e-3;
  }

  std::deque <int> hist0;
  std::deque <double> hist;
  std::deque <double> hist2;

  double step = min;
  double bin_width = (max - min) / number_of_bins;

  while (step <= max + 2 * bin_width)
  {
    hist0.push_back(0);
    hist.push_back(0);
    hist2.push_back(0);
    step += bin_width;
  }

  for (int i = 0; i < c.size(); i++)
  {
    double data = c[i];

    if (data > min && data <= max)
    {
      auto index = int((data - min) / bin_width);

      hist0[index]++;
      hist[index] += c2[i];
      hist2[index] += double(c[i] * c2[i]);
    }
  }

  for (int i = 0; i < hist.size() - 1; i++)
  {
    double x = hist2[i] / hist[i];
    auto y = double(hist[i] / hist0[i]);

    if (fabs(y) > 1e-10)
      out << x << "\t" << y << std::endl;
  }
  return 0;
}

int int_histogram(std::vector <int> &c, std::ostream & out)
{
  std::map<int, double> hist;

  double freq = 1 / double(c.size());

  for (int& i : c)
  {
    const auto itf = hist.find(i);
    if (itf == hist.end())
    {
      hist.emplace(i, 1.);
    }
    else
    {
      itf->second++;
    }
  }

  for (auto& it : hist)
  {
    it.second *= freq;
  }

  prints(hist, out);
  return 0;
}

int int_histogram(std::deque <int> &c, std::ostream & out)
{
  std::map<int, double> hist;

  const double freq = 1 / double(c.size());

  for (int& i : c)
  {
    const auto itf = hist.find(i);
    if (itf == hist.end())
      hist.emplace(i, 1.);
    else
      itf->second++;
  }

  for (auto& it : hist)
  {
    it.second *= freq;
  }
  prints(hist, out);
}

int int_histogram(int c, std::map<int, int> & hist)
{
  auto itf = hist.find(c);
  if (itf == hist.end())
    hist.insert(std::make_pair(c, 1));
  else
    itf->second++;
  return 0;
}

int int_histogram(int c, std::map<int, double> & hist, double w)
{
  auto itf = hist.find(c);
  if (itf == hist.end())
    hist.emplace(c, w);
  else
    itf->second += w;
  return 0;
}

int print_cumulative(std::deque<double> & kws, const std::string& file, int number_of_step)
{
  std::ofstream expout(file);
  std::sort(kws.begin(), kws.end());

  int step = (kws.size() - 1) / number_of_step;
  step = std::max(step, 1);

  for (int i = 0; i < kws.size(); i++) if (i%step == 0)
    expout << kws[i] << " " << double(i + 1) / (kws.size()) << std::endl;
  return 0;
}

int print_cumulative(std::deque<int> & kws, const std::string& file, int number_of_step)
{
  std::ofstream expout(file);
  sort(kws.begin(), kws.end());

  int step = (kws.size() - 1) / number_of_step;
  step = std::max(step, 1);

  for (int i = 0; i < kws.size(); i++) if (i%step == 0)
    expout << kws[i] << " " << double(i + 1) / (kws.size()) << std::endl;

  return 0;
}

int print_cumulative(std::vector<double> & kws, const std::string& file, int number_of_step)
{
  std::ofstream expout(file);
  sort(kws.begin(), kws.end());

  int step = (kws.size() - 1) / number_of_step;
  step = std::max(step, 1);

  for (int i = 0; i < kws.size(); i++)
  {
    if (i%step == 0)
    {
      expout << kws[i] << " " << double(i + 1) / (kws.size()) << std::endl;
    }
  }

  return 0;
}

int print_cumulative(std::vector<int> & kws, const std::string& file, int number_of_step)
{
  std::ofstream expout(file);
  sort(kws.begin(), kws.end());

  int step = (kws.size() - 1) / number_of_step;
  step = std::max(step, 1);

  for (int i = 0; i < kws.size(); i++) if (i%step == 0)
    expout << kws[i] << " " << double(i + 1) / (kws.size()) << std::endl;

  return 0;
}

int int_histogram(const std::string& infile, const std::string& outfile)
{
  // this makes a int_histogram of integers from a file

  std::ifstream ing(infile);
  std::deque<int> H;
  int h;
  while (ing >> h)
    H.push_back(h);

  std::ofstream outg(outfile);
  int_histogram(H, outg);

  return 0;
}

#endif
