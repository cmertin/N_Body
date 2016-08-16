#ifndef BARNESHUT_H
#define BARNESHUT_H
#include "planet.h"
#include "vector.h"
#include <iostream>
#include <algorithm>
#include <vector>

namespace barneshut
{
  template <typename T>
  void CenterOfMass(T *planets, unsigned int n)
    {
      if(n > 1)
	{
	  T newPlanet;
	  unsigned int depth = 0;
	  double mass = 0.0;
	  double x = 0.0;
	  double y = 0.0;
	  double z = 0.0;
	  double planetMass = 0.0;
	  Vector<double> pos;
	  for(unsigned int i = 0; i < n; ++i)
	    {
	      pos = (planets+i)->GetPos();
	      planetMass = (planets+i)->GetMass();
	      x += planetMass * pos.GetX();
	      y += planetMass * pos.GetY();
	      z += planetMass * pos.GetZ();
	      mass += planetMass;
	    }
	  x = x/mass;
	  y = y/mass;
	  z = z/mass;
	  depth = planets->GetDepth();
	  newPlanet.SetDepth(depth);
	  depth = planets->GetMaxDepth();
	  newPlanet.SetMaxDepth(depth);
	  newPlanet.SetPos(x, y, z);
	  newPlanet.SetMass(mass);
	  newPlanet.SetChildren(n);
	  for(unsigned int i = 0; i < n; ++i)
	    *(planets + i) = newPlanet;
	}
      return;
    }

  template <typename T>
  void Coarsen_CenterOfMass(T *planets, unsigned int n)
    {
      if(n > 0)
	{
	  T newPlanet;
	  unsigned int depth = 0;
	  double mass = 0.0;
	  double x = 0.0;
	  double y = 0.0;
	  double z = 0.0;
	  double planetMass = 0.0;
	  Vector<double> pos;
	  for(unsigned int i = 0; i < n; ++i)
	    {
	      pos = (planets+i)->GetPos();
	      planetMass = (planets+i)->GetMass();
	      x += planetMass * pos.GetX();
	      y += planetMass * pos.GetY();
	      z += planetMass * pos.GetZ();
	      mass += planetMass;
	    }
	  x = x/mass;
	  y = y/mass;
	  z = z/mass;
	  depth = planets->GetDepth() - 1;
	  newPlanet.SetDepth(depth);
	  depth = planets->GetMaxDepth();
	  newPlanet.SetMaxDepth(depth);
	  newPlanet.SetPos(x, y, z);
	  newPlanet.SetMass(mass);
	  newPlanet.SetChildren(n);
	  //newPlanet.SetCount(n);
	  //std::cout << "this count: " << newPlanet.GetCount() << std::endl;
	  for(unsigned int i = 0; i < n; ++i)
	    *(planets + i) = newPlanet;
	  
	  return;
	}
    }
}
#endif
