#ifndef BH_HERMITE_H
#define BH_HERMITE_H

#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>
#include "vector.h"
#include "planet.h"

template <typename T>
void Hermite(std::vector<Planet<T> > &planets, const double dt)
{
  const double G = 4 * M_PI * M_PI;
  Vector<T> zero;
  Vector<T> r;
  Vector<T> v;
  Vector<T> a;
  Vector<T> jk;
  int threads = omp_get_max_threads();
  
  #pragma omp parallel num_threads(threads) private(r, v, a, jk)
  for(int i = 0; i < planets.size(); ++i)
    {
      r = planets[i].GetPos();
      v = planets[i].GetVel();
      a = planets[i].GetAcc();
      jk = planets[i].GetJerk();
      planets[i].SetPosOld(r);
      planets[i].SetVelOld(v);
      planets[i].SetAccOld(a);
      planets[i].SetJerkOld(jk);
      v += a * dt + jk * (dt*dt)/2.0;
      r += v * dt + a * ((dt*dt)/2.0) + jk*((dt*dt*dt)/6.0);
      planets[i].SetPos(r);
      planets[i].SetVel(v);
    }
  
  #pragma omp parallel num_threads(threads) private(r, v, a, jk)
  for(int i = 0; i < planets.size(); ++i)
    {
      a = zero;
      jk = zero;
      for(int j = 0; j < planets.size(); ++j)
	{
	  if(i == j)
	    continue;
	  Vector<T> rji = planets[j].GetPos() - planets[i].GetPos();
	  Vector<T> vji = planets[j].GetVel() - planets[i].GetVel();
	  T rji3 = pow(rji.Magnitude(),3);
	  T rji5 = pow(rji.Magnitude(),5);
	  T rv = rji * vji;
	  a += rji * planets[j].GetMass()/rji3;
	  jk += (vji/rji3 -  rji * ((3/rji5) * rv)) * planets[j].GetMass();
	}
      a = a * G;
      jk = jk * G;
      planets[i].SetAcc(a);
      planets[i].SetJerk(jk);
    }
  
  #pragma omp parallel num_threads(threads) private(r, v, a, jk)
  for(int i = 0; i < planets.size(); ++i)
    {
      v = planets[i].GetVelOld() + (planets[i].GetAcc() + planets[i].GetAccOld()) * (dt/2.0)
	+ (planets[i].GetJerk() - planets[i].GetJerkOld()) * (dt * dt) / 12.0;
      r = planets[i].GetPos() + (planets[i].GetVel() + planets[i].GetVelOld()) * (dt/2.0)
	+ (planets[i].GetAcc() - planets[i].GetAccOld()) * (dt * dt) / 12.0;
      planets[i].SetVel(v);
      planets[i].SetPos(r);
    }
  return;
}


#endif
