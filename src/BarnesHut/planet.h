#ifndef PLANET_H
#define PLANET_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include "vector.h"
#include "OctTree.h"

template <typename T>
class Planet
{
 public:
  Planet();
  Planet(T x, T y, T z, T mass);
  Planet(T x, T y, T z, T mass, unsigned int level);
  Planet(T x, T y, T z, T mass, unsigned int level, unsigned int maxDepth);
  ~Planet();

  friend std::ostream &operator<<(std::ostream &out, const Planet<T> &planet)
  {
    out << std::setprecision(5) << std::scientific << planet.mass << ',' 
	<< std::fixed << planet.pos.GetX() << ',' << planet.pos.GetY() 
	<< ',' << planet.pos.GetZ();
    return out;
  }

  bool operator<(const Planet<T> &p2) const;
  bool operator>(const Planet<T> &p2) const;
  bool operator==(const Planet<T> &p2) const;
  bool operator<=(const Planet<T> &p2) const;
  bool operator!=(const Planet<T> &p2) const;
  bool operator<(const ot::TreeNode &oct) const;
  bool operator>(const ot::TreeNode &oct) const;
  Planet<T> &operator=(const Planet<T> &rhs);

  Vector<T> GetPos() const;
  //Vector<int> GetIntPos();
  Vector<T> GetVel() const;
  Vector<T> GetAcc() const;
  Vector<T> GetJerk() const;
  Vector<T> GetPosOld() const;
  Vector<T> GetVelOld() const;
  Vector<T> GetAccOld() const;
  Vector<T> GetJerkOld() const;
  unsigned int GetDepth() const;
  unsigned int getLevel() const;
  unsigned int getMaxDepth() const;
  unsigned int GetMaxDepth() const;
  unsigned int getX() const;
  unsigned int getY() const;
  unsigned int getZ() const;
  T GetMass() const;
  ot::TreeNode GetOctant() const;
  ot::TreeNode Octant() const;
  ot::TreeNode Octant(const unsigned int &level) const;

  void SetPos(Vector<T> &rhs);
  void SetPos(T &x, T &y, T &z);
  void SetPosOld(Vector<T> &rhs);
  void SetAcc(Vector<T> &rhs);
  void SetAccOld(Vector<T> &rhs);
  void SetVel(Vector<T> &rhs);
  void SetVelOld(Vector<T> &rhs);
  void SetJerk(Vector<T> &rhs);
  void SetJerkOld(Vector<T> &rhs);
  void SetMass(T mass);
  void SetDepth(unsigned int &level);
  void SetMaxDepth(unsigned int &maxDepth);
  
 private:
  unsigned int level;
  unsigned int maxDepth;
  Vector<T> pos;
  Vector<int> intPos;
  Vector<T> pos_old;
  Vector<T> vel;
  Vector<T> vel_old;
  Vector<T> acc;
  Vector<T> acc_old;
  Vector<T> jerk;
  Vector<T> jerk_old;
  T mass;
};

template <typename T>
Planet<T>::Planet()
{
  Vector<T> zero;
  this->pos = zero;
  this->pos_old = zero;
  this->vel = zero;
  this->vel_old = zero;
  this->acc = zero;
  this->acc_old = zero;
  this->jerk = zero;
  this->jerk_old = zero;
  this->mass = 0;
  this->level = 0;
}

template <typename T>
Planet<T>::Planet(T x, T y, T z, T mass)
{
  Vector<T> A(x, y, z);
  Vector<T> zero;
  this->pos = A;
  this->pos_old = zero;
  this->vel = zero;
  this->vel_old = zero;
  this->acc = zero;
  this->acc_old = zero;
  this->jerk = zero;
  this->jerk_old = zero;
  this->mass = mass;
  this->level = 0;
}

template <typename T>
Planet<T>::Planet(T x, T y, T z, T mass, unsigned int level)
{
  Vector<T> A(x, y, z);
  Vector<T> zero;
  this->pos = A;
  this->pos_old = zero;
  this->vel = zero;
  this->vel_old = zero;
  this->acc = zero;
  this->acc_old = zero;
  this->acc_old = zero;
  this->jerk = zero;
  this->jerk_old = zero;
  this->mass = mass;
  this->level = level;
  this->maxDepth = level;
}

template <typename T>
Planet<T>::Planet(T x, T y, T z, T mass, unsigned int level, unsigned int maxDepth)
{
  Vector<T> A(x, y, z);
  Vector<T> zero;
  this->pos = A;
  this->pos_old = zero;
  this->vel = zero;
  this->vel_old = zero;
  this->acc = zero;
  this->acc_old = zero;
  this->acc_old = zero;
  this->jerk = zero;
  this->jerk_old = zero;
  this->mass = mass;
  this->level = level;
  this->maxDepth = maxDepth;
}

template <typename T>
Planet<T>::~Planet()
{

}

template <typename T>
bool Planet<T>::operator<(const Planet<T> &p2) const
{
  unsigned int x1 = getX();
  unsigned int x2 = p2.getX();
  unsigned int y1 = getY();
  unsigned int y2 = p2.getY();
  unsigned int z1 = getZ();
  unsigned int z2 = p2.getZ();
  unsigned int depth1 = GetDepth();
  unsigned int depth2 = p2.GetDepth();
  unsigned int maxDepth1 = GetMaxDepth();
  unsigned int maxDepth2 = p2.GetMaxDepth();

  ot::TreeNode lhs(x1, y1, z1, depth1, 3, maxDepth1);
  ot::TreeNode rhs(x2, y2, z2, depth2, 3, maxDepth2);

  return lhs < rhs;
}

template <typename T>
bool Planet<T>::operator>(const Planet<T> &p2) const
{
  unsigned int x1 = getX();
  unsigned int x2 = p2.getX();
  unsigned int y1 = getY();
  unsigned int y2 = p2.getY();
  unsigned int z1 = getZ();
  unsigned int z2 = p2.getZ();
  unsigned int depth1 = GetDepth();
  unsigned int depth2 = p2.GetDepth();
  unsigned int maxDepth1 = GetMaxDepth();
  unsigned int maxDepth2 = p2.GetMaxDepth();

  ot::TreeNode lhs(x1, y1, z1, depth1, 3, maxDepth1);
  ot::TreeNode rhs(x2, y2, z2, depth2, 3, maxDepth2);

  return lhs > rhs;
}

template <typename T>
bool Planet<T>::operator==(const Planet<T> &p2) const
{
  unsigned int x1 = getX();
  unsigned int x2 = p2.getX();
  unsigned int y1 = getY();
  unsigned int y2 = p2.getY();
  unsigned int z1 = getZ();
  unsigned int z2 = p2.getZ();
  unsigned int depth1 = GetDepth();
  unsigned int depth2 = p2.GetDepth();
  unsigned int maxDepth1 = GetMaxDepth();
  unsigned int maxDepth2 = p2.GetMaxDepth();

  ot::TreeNode lhs(x1, y1, z1, depth1, 3, maxDepth1);
  ot::TreeNode rhs(x2, y2, z2, depth2, 3, maxDepth2);

  return lhs == rhs;
}

template <typename T>
bool Planet<T>::operator<=(const Planet<T> &p2) const
{
  unsigned int x1 = getX();
  unsigned int x2 = p2.getX();
  unsigned int y1 = getY();
  unsigned int y2 = p2.getY();
  unsigned int z1 = getZ();
  unsigned int z2 = p2.getZ();
  unsigned int depth1 = GetDepth();
  unsigned int depth2 = p2.GetDepth();
  unsigned int maxDepth1 = GetMaxDepth();
  unsigned int maxDepth2 = p2.GetMaxDepth();

  ot::TreeNode lhs(x1, y1, z1, depth1, 3, maxDepth1);
  ot::TreeNode rhs(x2, y2, z2, depth2, 3, maxDepth2);

  return lhs <= rhs;
}

template <typename T>
bool Planet<T>::operator!=(const Planet<T> &p2) const
{
  unsigned int x1 = getX();
  unsigned int x2 = p2.getX();
  unsigned int y1 = getY();
  unsigned int y2 = p2.getY();
  unsigned int z1 = getZ();
  unsigned int z2 = p2.getZ();
  unsigned int depth1 = GetDepth();
  unsigned int depth2 = p2.GetDepth();
  unsigned int maxDepth1 = GetMaxDepth();
  unsigned int maxDepth2 = p2.GetMaxDepth();

  ot::TreeNode lhs(x1, y1, z1, depth1, 3, maxDepth1);
  ot::TreeNode rhs(x2, y2, z2, depth2, 3, maxDepth2);

  return lhs != rhs;
}

template <typename T>
bool Planet<T>::operator<(const ot::TreeNode &oct) const
{
  ot::TreeNode current = Octant();
  
  return current < oct;
}

template <typename T>
bool Planet<T>::operator>(const ot::TreeNode &oct) const
{
  ot::TreeNode current = Octant();
  
  return oct < current;
}

template <typename T>
Planet<T> &Planet<T>::operator=(const Planet<T> &rhs)
{
  Vector<T> tempVec = rhs.GetPos();
  unsigned int tempInt = 0;
  T tempMass = 0;
  
  SetPos(tempVec);
  tempVec = rhs.GetPosOld();
  SetPosOld(tempVec);
  tempVec = rhs.GetVel();
  SetVel(tempVec);
  tempVec = rhs.GetVelOld();
  SetVelOld(tempVec);
  tempVec = rhs.GetAcc();
  SetAcc(tempVec);
  tempVec = rhs.GetAccOld();
  SetAccOld(tempVec);
  tempVec = rhs.GetJerk();
  SetJerk(tempVec);
  tempVec = rhs.GetJerkOld();
  SetJerkOld(tempVec);
  tempMass = rhs.GetMass();
  SetMass(tempMass);
  tempInt = rhs.GetDepth();
  SetDepth(tempInt);
  tempInt = rhs.GetMaxDepth();
  SetMaxDepth(tempInt);

  return *this;
}

template <typename T>
Vector<T> Planet<T>::GetPos() const
{
  return this->pos;
}

/*
template <typename T>
Vector<int> Planet<T>::GetIntPos()
{
  return this->intPos;
}
*/

template <typename T>
Vector<T> Planet<T>::GetPosOld() const
{
  return this->pos_old;
}

template <typename T>
Vector<T> Planet<T>::GetVel() const
{
  return this->vel;
}

template <typename T>
Vector<T> Planet<T>::GetVelOld() const
{
  return this->vel_old;
}

template <typename T>
Vector<T> Planet<T>::GetAcc() const
{
  return this->acc;
}

template <typename T>
Vector<T> Planet<T>::GetAccOld() const
{
  return this->acc_old;
}

template <typename T>
Vector<T> Planet<T>::GetJerk() const
{
  return this->jerk;
}

template <typename T>
Vector<T> Planet<T>::GetJerkOld() const
{
  return this->jerk_old;
}

template <typename T>
T Planet<T>::GetMass() const
{
  return this->mass;
}

template <typename T>
ot::TreeNode Planet<T>::GetOctant() const
{
  Vector<T> currPos = this->pos;
  unsigned int x = (unsigned int)currPos.GetX();//getX();
  unsigned int y = (unsigned int)currPos.GetY();//getY();
  unsigned int z = (unsigned int)currPos.GetZ();//getZ();
  unsigned int depth = GetDepth();
  unsigned int maxDepth = GetMaxDepth();
  
  ot::TreeNode octant(x, y, z, maxDepth, 3, maxDepth);
  //std::cout << '\t' << currPos << '\t' << depth << '\t' << maxDepth << '\t' << octant << std::endl;
  return octant;
}

template <typename T>
ot::TreeNode Planet<T>::Octant() const
{
  ot::TreeNode oct = GetOctant();
  unsigned int depth = GetDepth();
  return oct.getAncestor(depth);
  /*
  ot::TreeNode oct = GetOctant();
  unsigned int x = oct.getX();
  unsigned int y = oct.getY();
  unsigned int z = oct.getZ();
  unsigned int depth = GetDepth();
  unsigned int maxDepth = GetMaxDepth();
  
  ot::TreeNode octant(x, y, z, depth, 3, maxDepth);
  return octant;
  */
}

template <typename T>
ot::TreeNode Planet<T>::Octant(const unsigned int &level) const
{
  ot::TreeNode oct = GetOctant();
  return oct.getAncestor(level);
}

template <typename T>
unsigned int Planet<T>::GetDepth() const
{
  return this->level;
}

template <typename T>
unsigned int Planet<T>::getLevel() const
{
  return this->level;
}

template <typename T>
unsigned int Planet<T>::getMaxDepth() const
{
  return this->maxDepth;
}

template <typename T>
unsigned int Planet<T>::GetMaxDepth() const
{
  return this->maxDepth;
}

template <typename T>
unsigned int Planet<T>::getX() const
{

  /*
  Vector<T> pos = GetPos();
  unsigned int level = GetDepth();
  unsigned int maxDepth = GetMaxDepth();
  unsigned int x = pos.GetX()/(1u<<maxDepth) * ((1u<<level)-1);
  */
  ot::TreeNode octant = GetOctant();
  return octant.getAncestor(level).getX();
  //return x;
}

template <typename T>
unsigned int Planet<T>::getY() const
{
  /*
  Vector<T> pos = GetPos();
  unsigned int level = GetDepth();
  unsigned int maxDepth = GetMaxDepth();
  unsigned int y = pos.GetY()/(1u<<maxDepth) * ((1u<<level)-1);
  return y;
  */
  ot::TreeNode octant = GetOctant();
  return octant.getAncestor(level).getY();
}

template <typename T>
unsigned int Planet<T>::getZ() const
{
  /*
  Vector<T> pos = GetPos();
  unsigned int level = GetDepth();
  unsigned int maxDepth = GetMaxDepth();
  unsigned int z = pos.GetZ()/(1u<<maxDepth) * ((1u<<level)-1);
  return z;
  */
  ot::TreeNode octant = GetOctant();
  return octant.getAncestor(level).getZ();
}

template <typename T>
void Planet<T>::SetPos(Vector<T> &rhs)
{
  this->pos = rhs;
}

template <typename T>
void Planet<T>::SetPos(T &x, T &y, T &z)
{
  Vector<T> newPos(x, y, z);
  this->pos = newPos;
}

template <typename T>
void Planet<T>::SetPosOld(Vector<T> &rhs)
{
  this->pos_old = rhs;
}

template <typename T>
void Planet<T>::SetVel(Vector<T> &rhs)
{
  this->vel = rhs;
}

template <typename T>
void Planet<T>::SetVelOld(Vector<T> &rhs)
{
  this->vel_old = rhs;
}

template <typename T>
void Planet<T>::SetAcc(Vector<T> &rhs)
{
  this->acc = rhs;
}

template <typename T>
void Planet<T>::SetAccOld(Vector<T> &rhs)
{
  this->acc_old = rhs;
}

template <typename T>
void Planet<T>::SetJerk(Vector<T> &rhs)
{
  this->jerk = rhs;
}

template <typename T>
void Planet<T>::SetJerkOld(Vector<T> &rhs)
{
  this->jerk_old = rhs;
}

template <typename T>
void Planet<T>::SetMass(T mass)
{
  this->mass = mass;
}

template <typename T>
void Planet<T>::SetDepth(unsigned int &level)
{
  this->level = level;
}

template <typename T>
void Planet<T>::SetMaxDepth(unsigned int &maxDepth)
{
  this->maxDepth = maxDepth;
}

#endif
