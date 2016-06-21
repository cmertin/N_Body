#ifndef ENERGYCONSERVATION_H
#define ENERGYCONSERVATION_H
#include "planet.h"
#include "vector.h"

template <typename T>
T EnergyConservation(std::vector<Planet<T> > &planets)
{
  const double G = 6.67408e-11;
  const double solarMass = 1.98855e30;
  const double AUMeters = 1.496e11;
  const double AUyr2ms = 4744;
  T energy = 0;
  T potentialEnergy = 0;
  T kineticEnergy = 0;

  // Calculate potential energy
  for(int i = 0; i < planets.size(); ++i)
    {
      for(int j = 0; j < planets.size(); ++j)
	{
	  if(i == j)
	    continue;
	  Vector<T> r1 = planets[i].GetPos() * AUMeters;
	  Vector<T> r2 = planets[j].GetPos() * AUMeters; 
	  T distance = (r2 - r1).Magnitude();
	  T m1 = planets[i].GetMass() * solarMass;
	  T m2 = planets[j].GetMass() * solarMass;
	  potentialEnergy += - (G * m1 * m2)/distance;
	}
    }

  // Calculate kinetic energy
  for(int i = 0; i < planets.size(); ++i)
    {
      Vector<T> vel = planets[i].GetVel() * AUyr2ms;
      T vx2 = vel.GetX() * vel.GetX();
      T vy2 = vel.GetY() * vel.GetY();
      T vz2 = vel.GetZ() * vel.GetZ();
      T mass = planets[i].GetMass() * solarMass;
      kineticEnergy += 0.5 * mass * (vx2 + vy2 + vz2);
    }

  energy = kineticEnergy + potentialEnergy;

  return energy;
}


#endif
