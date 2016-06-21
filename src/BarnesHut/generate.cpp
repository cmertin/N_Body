#include <iostream>
#include <random>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include "planet.h"

using namespace std;

bool Exists(const string &data_dir);
void WriteBodies(vector<Planet<double> > &planets, string &data_dir, double minX, double maxX);

int main(int argc, char *argv[])
{
  string data_dir = "data/";
  int maxNum = 10000;
  // Random number generator
  int maxDepth = (1<<8)-2;
  double minX = 0;
  double maxX = maxDepth;
  double minMass = 10e-8;
  double maxMass = 10;
  random_device rd;
  default_random_engine generator(rd()); // rd() provides a random seed
  uniform_real_distribution<double> posGen(minX, maxX);
  uniform_real_distribution<double> massGen(minMass, maxMass);

  vector<Planet<double> > planets;

  if(Exists(data_dir) == false)
    {
      cerr << "data directory \"" << data_dir << "\" not found. Creating..." << endl;
      string sys_cmd = "mkdir ";
      sys_cmd.append(data_dir);
      system(sys_cmd.c_str());
    }

  if(argc == 1)
    {
      for(int i = 1; i <= maxNum; ++i)
	{
	  double mass = massGen(generator);
	  double x = posGen(generator);
	  double y = posGen(generator);
	  double z = posGen(generator);
	  Planet<double> temp(x, y, z, mass);
	  planets.push_back(temp);
	  
	  if(i == 10)
	    WriteBodies(planets, data_dir, minX, maxX);
	  if(i == 100)
	    WriteBodies(planets, data_dir, minX, maxX);
	  
	  if(i % 1000 == 0)
	    WriteBodies(planets, data_dir, minX, maxX);
	}
    }
  if(argc == 2)
    {
      maxNum = atoi(argv[1]);
      for(int i = 1; i <= maxNum; ++i)
	{
	  double mass = massGen(generator);
	  double x = posGen(generator);
	  double y = posGen(generator);
	  double z = posGen(generator);
	  Planet<double> temp(x, y, z, mass);
	  planets.push_back(temp);

	  if(i == maxNum)
	    WriteBodies(planets, data_dir, minX, maxX);
	}
    }
  else
    {
      cerr << "Requires 0 or 1 arguments" << endl;
      cerr << "argument: Number_of_bodies" << endl;
      return -1;
    }
  
  return 0;
}

bool Exists(const string &data_dir)
{
  struct stat st;
  if(stat(data_dir.c_str(),&st) == 0)
    {
      if(st.st_mode & S_IFDIR != 0)
	return true;
    }
  else
    return false;
}

void WriteBodies(vector<Planet<double> > &planets, string &data_dir, double minX, double maxX)
{
  string filename = data_dir;
  filename.append(to_string(planets.size()));
  filename.append("_bodies.dat");
  ofstream output;
  output.open(filename.c_str());
  output << planets.size() << endl;
  output << minX << ',' << maxX << ',' << minX << ',' << maxX << ',' << minX << ',' << maxX << endl;

  for(int i = 0; i < planets.size(); ++i)
    output << planets[i] << endl;
  output.close();
}
    
