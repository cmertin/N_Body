#include <iostream>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <random>
//#include "mpi.h"
#include "msd.h"
#include "vector.h"
#include "planet.h"
#include "bh_hermite.h"
#include "EnergyConservation.h"
#include "OctTree.h"
#include "barneshut.h"


using namespace std;

struct Dimensions
{
  double minX;
  double maxX;
  double minY;
  double maxY;
  double minZ;
  double maxZ;
  int numPlanets;
};

void ProgressCheck(int &nodes, int &node_id);
void ReadFile(vector<Planet<double> > &planets, Dimensions &dim);
template <typename T>
bool PlanetInOctant(ot::TreeNode &octant, Planet<T> &planet);
template <typename T>
void CenterOfMass(vector<Planet<T> > &inOctant, Planet<T> &newPlanet, unsigned int &maxDepth);
inline bool NotChild(vector<ot::TreeNode> &octants, ot::TreeNode &newOct);
inline bool NotAncestor(vector<ot::TreeNode> &octants, ot::TreeNode &newOct);
void CoarsenIncomplete(vector<ot::TreeNode> &in, vector<ot::TreeNode> &out);
template <typename T>
void CoarsenParticles(vector<Planet<T> > &in, vector<Planet<T> > &out);//, vector<ot::TreeNode> &octants);
template <typename T>
void RemoveDuplicates(vector<Planet<T> > &in);
template <typename T>
bool CheckCoarsening(vector<Planet<T> > &fine, vector<Planet<T> > &coarse);
template <typename T>
bool operator<(const ot::TreeNode &oct, const Planet<T> &planet);

int main(int argc, char *argv[])
{
  
  MPI::Init(argc, argv);

  chrono::time_point<chrono::system_clock> start, end;
  chrono::duration<double> elapsed_seconds;

  int node_id;
  int nodes;
  int localSize;
  int globalSize;
  bool divisible = true;
  int remainder = 0;
  Dimensions dimension;
  unsigned int dim = 3;
  unsigned int maxDepth = 8;
  unsigned int maxNumPts = 1;
  const int MAX_COARSEN = 10;
  const int MAX_PER_OCT = 1;
  double energyInit = 0;
  double energyFinal = 0;
  vector<int> boundaryIndex;
  double gSize[3] = {1.0, 1.0, 1.0};
  double startTime = 0;
  double endTime = 0;
  double minTime = 0;
  double maxTime = 0;
  int numParticles;
  vector<vector<Planet<double> > > allPlanets;
  allPlanets.resize(1);
  if(argc > 1)
    numParticles = atoi(argv[1]);
  else
    numParticles = 1000;

  if(argc == 1)
    dimension.numPlanets = 1000;
  else
    dimension.numPlanets = atoi(argv[1]);

  double dt = 0.001;
  int years = 2;
  int steps = years/dt;
  int numPts = 0;
  steps = 5;

  _InitializeHcurve(dim);

  nodes = MPI::COMM_WORLD.Get_size();
  node_id = MPI::COMM_WORLD.Get_rank();

  // Generate random particles
  random_device rd;
  default_random_engine generator(rd());
  uniform_real_distribution<double> posGen(1, 250);
  uniform_real_distribution<double> massGen(.3, 10);

  //vector<ot::TreeNode> temp;
  vector<Planet<double> > planets;
  vector<ot::TreeNode> octants;
  vector<vector<ot::TreeNode> > tree;
  vector<ot::TreeNode> temp;
  vector<ot::TreeNode> coarseTemp;
  
  // Build the Particles
  for(int i = 0; i < numParticles; ++i)
    {
      double x = posGen(generator);
      double y = posGen(generator);
      double z = posGen(generator);
      double mass = massGen(generator);
      
      Planet<double> tempPlanet(x, y, z, mass, maxDepth);
      allPlanets[0].push_back(tempPlanet);
      planets.push_back(tempPlanet);
    }

  cout << "Finished Building Particles" << endl;

  // Sort the particles based off of Hilbert/Morton ordering
  #ifndef __SILENT_MODE__
  {
    cout << "Before Sorting, is sorted: " << boolalpha 
	 << seq::test::isSorted(allPlanets[0]) << endl;
    cout << "Before remove: " << allPlanets[0].size() << endl;

    cout << "Sorting..." << endl;
  }
  #endif
  start = chrono::system_clock::now();

  MSD_Sort_rd(allPlanets[0], 0, maxDepth, MAX_PER_OCT);
  //SFC::seqSort::SFC_3D_msd_sort_rd(&(*(planets.begin())), planets.size(), 0, maxDepth, maxDepth);
  #ifndef __SILENT_MODE__
  {
    cout << "\tFinished Initial Sorting" << endl;
 
    cout << "After remove: " << allPlanets[0].size() << endl;

    cout << "After Sorting, is sorted: " << boolalpha 
	 << seq::test::isSorted(allPlanets[0]) << endl;

    cout << "Finished Sorting Planets" << endl;
  
  // Build the octants based on sorted
  ot::TreeNode tempOctant;
  for(int i = 0; i < allPlanets[0].size(); ++i)
    {
      tempOctant = allPlanets[0][i].Octant();
      octants.push_back(tempOctant);
    }

  cout << "Finished Building Octants" << endl;

  treeNodesTovtk(octants, node_id, "finest_octants");
  }
  #endif
  //temp.insert(temp.end(), octants.begin(), octants.begin() + 10);
  //treeNodesTovtk(temp, node_id, "localized_octants-fine");

  // Coarsen - check it works
  //ot::coarsenOctree(completeOctree, coarseTemp);
  //CoarsenIncomplete(temp, coarseTemp);
  
  //treeNodesTovtk(coarseTemp, node_id, "localized_octants-coarse");

  // Coarsens all the octants
  //CoarsenIncomplete(octants, coarseTemp);
  
  //treeNodesTovtk(coarseTemp, node_id, "octants-coarse");
  
  // Coarsens the planets
  vector<Planet<double> > coarsenPlanets;
  string basename = "coarsened_";
  string output;
  int level = maxDepth;
  int itr = 0;
  planets = allPlanets[0];
  do
    {
      allPlanets.resize(allPlanets.size() + 1);
      --level;
      output = basename;
      output.append(to_string(level));
      //CoarsenParticles(planets, coarsenPlanets);//, coarseTemp);
      CoarsenParticles(allPlanets[itr], allPlanets[itr+1]);
      //MSD_Sort_rd(coarsenPlanets, 0, maxDepth);
      //MSD_Sort(coarsenPlanets, maxDepth);
      MSD_Sort(allPlanets[itr+1], maxDepth);
      //SFC::seqSort::SFC_3D_msd_sort_rd(&(*(coarsenPlanets.begin())), coarsenPlanets.size(), 0, maxDepth, maxDepth);
      
      octants.clear();
      for(int i = 0; i < allPlanets[itr+1].size(); ++i)
	{
	  //cout << "i: " << i << '\t' << allPlanets[itr+1][i].Octant() << endl;
	  octants.push_back(allPlanets[itr+1][i].Octant());
	  //octants.push_back(coarsenPlanets[i].Octant());
	}
      
      treeNodesTovtk(octants, node_id, output.c_str());
      
      // Checks Coarsening
      //cout << "Coarsened: " << boolalpha << CheckCoarsening(allPlanets[itr], allPlanets[itr+1]) << endl;
      
      cout << "Level " << level << ": " << allPlanets[itr+1].size() << endl;
      
      ++itr;
    }while(level > 2);//allPlanets[allPlanets.size()-1].size() > MAX_COARSEN);
  
  cout << "Total Levels: " << allPlanets.size() << endl;

  cout << "Last level: " << allPlanets[allPlanets.size()-1].size() << endl;
  
  end = chrono::system_clock::now();
  elapsed_seconds = end - start;
  cout << "Time: " << elapsed_seconds.count() << " seconds" << endl;



  // For the down sweep
  Planet<double> current = allPlanets[0][500];
  ot::TreeNode rootOct(0, 0, 0, 0, dim, maxDepth);
  vector<ot::TreeNode> currOct;
  currOct.push_back(current.Octant(maxDepth - 2));
  vector<Planet<double> > accumulate;
  vector<ot::TreeNode> lower;
  vector<ot::TreeNode> upper;
  vector<ot::TreeNode> neighbors = (current.Octant(2)).getAllNeighbours();
  neighbors.push_back(current.Octant(2));
  sort(neighbors.begin(), neighbors.end());
  //for(int i = 0; i < neighbors.size(); ++i)
  //cout << "i: " << i << '\t' << neighbors[i] << endl;
  neighbors.erase(remove(neighbors.begin(), neighbors.end(), rootOct), neighbors.end());
  ot::TreeNode minLocalOct = neighbors.front();
  ot::TreeNode maxLocalOct = neighbors.back();
  treeNodesTovtk(neighbors, node_id, "local");
  treeNodesTovtk(currOct, node_id, "current");
  cout << minLocalOct << '\t' << maxLocalOct << endl;
  int prevMinIndex = 0;
  int prevMaxIndex = 0;

  auto lowerBound = lower_bound(allPlanets[maxDepth-2].begin(), allPlanets[maxDepth-2].end(), minLocalOct);
  auto upperBound = upper_bound(allPlanets[maxDepth-2].begin(), allPlanets[maxDepth-2].end(), maxLocalOct);

  for(auto i = allPlanets[maxDepth-2].begin(); i < lowerBound; ++i)
    accumulate.push_back(*i);

  for(auto i = upperBound; i < allPlanets[maxDepth-2].end(); ++i)
    accumulate.push_back(*i);
  
  auto allItr = lowerBound;
  auto localItr = neighbors.begin();

  while(allItr != upperBound)
    {
      if(*localItr == (*allItr).Octant())
	{
	  ++localItr;
	  ++allItr;
	}
      else
	{
	  accumulate.push_back(*allItr);
	  ++allItr;
	}
    }

  for(int i = 0; i < accumulate.size(); ++i)
    upper.push_back(accumulate[i].Octant());
  

  treeNodesTovtk(upper, node_id, "upper");
  /*

  auto allItr = allPlanets[maxDepth-2].begin();
  auto localItr = neighbors.begin();

  while(allItr != allPlanets[maxDepth-2].end())
    {
      if(*localItr == (*allItr).Octant())
	{
	  ++localItr;
	  ++allItr;
	}
      else
	{
	  accumulate.push_back((*allItr).Octant());
	  ++allItr;
	}
    }
  treeNodesTovtk(accumulate, node_id, "accumulate");
  /*
  for(auto i = allPlanets[maxDepth-2].begin(); i < allPlanets[maxDepth-2].end(); ++i)
    cout << (*i).Octant() << endl;

  auto lowerBound = lower_bound(allPlanets[maxDepth-2].begin(), allPlanets[maxDepth-2].end(), minLocalOct);
  auto upperBound = upper_bound(allPlanets[maxDepth-2].begin(), allPlanets[maxDepth-2].end(), maxLocalOct);
  cout << "Lower: " << minLocalOct << '\t' << (*(lowerBound-1)).Octant() << endl;
  cout << "\t distance: " << (lowerBound - 1) - allPlanets[maxDepth-2].begin() << endl;
  for(auto i = allPlanets[maxDepth-2].begin(); i < lowerBound; ++i)
    lower.push_back((*i).Octant());

  cout << "Upper: " << maxLocalOct << '\t' << (*upperBound).Octant() << endl;
  cout << "\t distance: " << allPlanets[maxDepth-2].end() - upperBound << endl;
  cout << "\t neighbors: " << neighbors.size() + 1 << endl;
  for(auto i = upperBound-8; i < allPlanets[maxDepth-2].end(); ++i)
    upper.push_back((*i).Octant());
  
  treeNodesTovtk(lower, node_id, "lower");
  treeNodesTovtk(upper, node_id, "upper");
  ot::TreeNode equalCheck;
  for(int i = 0; i < coarsenPlanets.size(); ++i)
    {
      equalCheck = coarsenPlanets[i].GetOctant();
      if(equalCheck != coarseTemp[i])
	{
	  cout << "Failed for i = " << i << '\t' << equalCheck << '\t' << coarseTemp[i] << endl;
	}
    }
  */
  MPI::Finalize();
  return 0;
}

void ProgressCheck(int &nodes, int &node_id)
{
  for(int i = 0; i < nodes; ++i)
    {
      if(node_id == i)
	cout << "Node: " << node_id << endl;
      MPI::COMM_WORLD.Barrier();
    }

  return;
}

void ReadFile(vector<Planet<double> > &planets, Dimensions &dim)
{
  string filename = "data/";
  filename.append(to_string(dim.numPlanets));
  filename.append("_bodies.dat");
  vector<double> dimensions;
  vector<double> bodies;
  string delimeter = ",";
  string token;
  fstream data;
  string line;
  size_t pos = 0;
  bool fileExists = ifstream(filename.c_str()).good();

  if(fileExists == false)
    {
      cerr << "File \"" << filename << "\" does not exist. Exiting..." << endl;
      MPI::COMM_WORLD.Abort(-1);
      MPI::Finalize();
      exit(-1);
    }

  data.open(filename.c_str());
  getline(data, line);
  getline(data, line);
  while((pos = line.find(delimeter)) != string::npos)
    {
      token = line.substr(0, pos);
      dimensions.push_back(stod(token));
      line.erase(0, pos + delimeter.length());
    }

  dim.minX = dimensions[0];
  dim.maxX = dimensions[1];
  dim.minY = dimensions[2];
  dim.maxY = dimensions[3];
  dim.minZ = dimensions[4];
  dim.maxZ = stod(line);

  while(getline(data, line))
    {
      pos = 0;
      bodies.clear();
      while((pos = line.find(delimeter)) != string::npos)
	{
	  token = line.substr(0, pos);
	  bodies.push_back(stod(token));
	  line.erase(0, pos + delimeter.length());
	}
      bodies.push_back(stod(line));
      Planet<double> temp(bodies[1], bodies[2], bodies[3], bodies[0]);
      planets.push_back(temp);
    }

  data.close();
  #ifndef __SILENT_MODE__
  {
    cout << "\tFinished Reading Data" << endl;
  }
  #endif

  return;
}

template <typename T>
bool PlanetInOctant(ot::TreeNode &octant, Planet<T> &planet)
{
  Vector<int> pos = planet.GetIntPos();
  bool inX = false;
  bool inY = false;
  bool inZ = false;
  
  if(octant.minX() <= pos.GetX() && pos.GetX() < octant.maxX())
    inX = true;
  
  if(octant.minY() <= pos.GetY() && pos.GetY() < octant.maxY())
    inY = true;

  if(octant.minZ() <= pos.GetZ() && pos.GetZ() < octant.maxZ())
    inZ = true;
  
  return (inX && inY && inZ);
}

template <typename T>
void CenterOfMass(vector<Planet<T> > &inOctant, Planet<T> &newPlanet, unsigned int &depth)
{
  if(inOctant.size() > 0)
    {
      T mass = 0.0;
      T x = 0.0;
      T y = 0.0;
      T z = 0.0;
      Vector<T> pos;
      T planetMass;
      for(int i = 0; i < inOctant.size(); ++i)
	{
	  pos = inOctant[i].GetPos();
	  planetMass = inOctant[i].GetMass();
	  x += planetMass * pos.GetX();
	  y += planetMass * pos.GetY();
	  z += planetMass * pos.GetZ();
	  mass += planetMass;
	}
      x = x/mass;
      y = y/mass;
      z = z/mass;
      newPlanet.SetPos(x, y, z);
      newPlanet.SetMass(mass);
      newPlanet.SetDepth(depth);
    }
  return;
}

// Checks to see if the newOctant is a child of any of the current octants
inline bool NotChild(vector<ot::TreeNode> &octants, ot::TreeNode &newOct)
{
  for(int i = 0; i < octants.size(); ++i)
    {
      if(octants[i].isAncestor(newOct))
	return false;
    }
  return true;
}

// Checks to see if the newOctant is a parent of any of the current octants
inline bool NotAncestor(vector<ot::TreeNode> &octants, ot::TreeNode &newOct)
{
  for(int i = 0; i < octants.size(); ++i)
    {
      if(newOct.isAncestor(octants[i]))
	return false;
    }
  return true;
}

void CoarsenIncomplete(vector<ot::TreeNode> &in, vector<ot::TreeNode> &out)
{
  ot::TreeNode tempParent;

  out.clear();
  out.push_back(in[0].getParent());

  for(auto itr = in.begin()+1; itr != in.end(); ++itr)
    {
      tempParent = (*itr).getParent();
      if(out[out.size()-1] == tempParent)
	continue;
      else
	out.push_back(tempParent);
    }
  
  return;
}

template <typename T>
void CoarsenParticles(vector<Planet<T> > &in, vector<Planet<T> > &out)//, vector<ot::TreeNode> &octants)
{
  out.clear();

  Planet<T> tempPlanet;
  vector<Planet<T> > planets;
  unsigned int depth = ((*in.begin()).GetDepth()) - 1;
  ot::TreeNode prev = ((*(in.begin())).GetOctant().getAncestor(depth));
  ot::TreeNode curr;
  auto itr = in.begin() + 1;
  planets.push_back(*in.begin());

  while(itr != in.end())
    {
      curr = (*itr).GetOctant().getAncestor(depth);
      // add to vector 
      if(prev == curr)
	{
	  planets.push_back(*itr);
	  ++itr;
	}
      // perform coarsening
      else
	{
	  if(planets.size() == 0)
	    {
	      planets.push_back(*itr);
	    }
	  else
	    {
	      barneshut::Coarsen_CenterOfMass(&(*(planets.begin())), planets.size());
	      out.push_back(planets[0]);
	      planets.clear();
	      --itr;
	    }
	  prev = curr;
	  ++itr;
	}
    }

  if(planets.size() > 0)
    {
      barneshut::Coarsen_CenterOfMass(&(*(planets.begin())), planets.size());
      out.push_back(planets[0]);
      planets.clear();
    }
  
  return;
}

template <typename T>
void RemoveDuplicates(vector<Planet<T> > &in)
{
  // Erases all duplicates
  in.erase(unique(in.begin(), in.end()), in.end());

  return;
}

template <typename T>
bool CheckCoarsening(vector<Planet<T> > &fine, vector<Planet<T> > &coarse)
{
  auto fineItr = fine.begin();
  auto coarseItr = coarse.begin();
  ot::TreeNode fineOct;
  ot::TreeNode coarseOct;
  bool endFine = false;
  bool endCoarse = false;
  int indexFine = 0;
  int indexCoarse = 0;
  
  while(fineItr != fine.end() || coarseItr != coarse.end())
    {
      fineOct = (*fineItr).Octant();
      coarseOct = (*coarseItr).Octant();
      indexFine = fineItr - fine.begin() + 1;
      indexCoarse = coarseItr - coarse.begin() + 1;
      if(coarseOct.isAncestor(fineOct))
	++fineItr;
      else
	++coarseItr;
    }

  return true;
}

template <typename T>
bool operator<(const ot::TreeNode &oct, const Planet<T> &planet)
{
  ot::TreeNode current = planet.Octant();

  return oct < current;
}
