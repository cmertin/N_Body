#ifndef MSD_H
#define MSD_H

#include "OctTree.h"
#include "barneshut.h"
#include <cassert>


template<typename T>
void MSD_Sort(T *pNodes, DendroIntL n, unsigned int rot_id, unsigned int pMaxDepthBit, unsigned int pMaxDepth, int numPer = 1)
{
  register unsigned int cnum;
  register unsigned int cnum_prev = 0;
  //register unsigned int n=0;
  unsigned int rotation = 0;
  DendroIntL count[NUM_CHILDREN + 2] = {};
  unsigned int lev = pMaxDepth-pMaxDepthBit;
  pMaxDepthBit--;
  
  count[0] = 0;
  for(DendroIntL i = 0; i < n; ++i) 
    {
      if(lev==pNodes[i].getLevel())
	count[1]++;
      else
	{
	  cnum = (((pNodes[i].getZ() >> pMaxDepthBit) & 1u) << 2u) | (((pNodes[i].getY() >> pMaxDepthBit) & 1u) << 1u) | ((pNodes[i].getX() >>pMaxDepthBit) & 1u);
	  count[cnum + 2]++;
	}
    }
  
  DendroIntL loc[NUM_CHILDREN + 1];
  T unsorted[NUM_CHILDREN + 1];
  unsigned int live = 0;
  
  
  for(unsigned int i = 0; i < NUM_CHILDREN + 1; ++i) 
    {
      if(i==0)
	{
	  loc[0] = count[0];
	  count[1] += count[0];
	  unsorted[live] = pNodes[loc[0]];
	  if(loc[0] < count[1]) 
	    live++;
	}
      else
	{
	  cnum = (rotations[ROTATION_OFFSET * rot_id+ i-1] - '0');
	  (i>1) ? cnum_prev = ((rotations[ROTATION_OFFSET * rot_id+i-2] - '0')+2): cnum_prev=1;
	  loc[cnum+1] = count[cnum_prev];
	  count[cnum+2] += count[cnum_prev];
	  unsorted[live] = pNodes[loc[cnum+1]];
	  if(loc[cnum+1] < count[cnum+2])
	    live++;
	}
    }
  live--;
  
  for (DendroIntL i=0; i < n ; ++i) 
    {
      if(lev==unsorted[live].getLevel())
	{
	  pNodes[loc[0]++] = unsorted[live];
	  unsorted[live] = pNodes[loc[0]];
	  if((loc[0] == count[1])) 
	    live--;
	}
      else
	{
	  cnum = (((unsorted[live].getZ() >> pMaxDepthBit) & 1u) << 2u) | (((unsorted[live].getY() >> pMaxDepthBit) & 1u) << 1u) | ((unsorted[live].getX() >> pMaxDepthBit) & 1u);
	  pNodes[loc[(cnum + 1)]++] = unsorted[live];
	  unsorted[live] = pNodes[loc[cnum+1]];
	  if((loc[cnum+1] == count[cnum + 2])) 
	    live--;
	}
    }

  if(pMaxDepthBit > 0) 
    {
      for(unsigned int i = 1; i < NUM_CHILDREN + 1; i++) 
	{
	  cnum=(rotations[ROTATION_OFFSET*rot_id+i-1]-'0');
	  (i>1)? cnum_prev = ((rotations[ROTATION_OFFSET * rot_id+i-2] - '0')+2) : cnum_prev=1;
	  n = count[cnum+2] - count[cnum_prev];
	  if (n > numPer) 
	    {
	      rotation=HILBERT_TABLE[NUM_CHILDREN * rot_id + cnum];
	      MSD_Sort(pNodes + count[cnum_prev], n, rotation, (pMaxDepthBit), pMaxDepth, numPer);
	    }
	  else
	    {
	      // Perform the coarsening/averaging
	      barneshut::CenterOfMass(pNodes + count[cnum_prev], n);
	    }
	}
    }
}

template <typename T>
void MSD_Sort_rd(std::vector<T>&pNodes,unsigned int rot_id,unsigned int pMaxDepth, int numPer = 1, bool isSorted = false)
{
  if(!isSorted)
    {
      MSD_Sort((&(*(pNodes.begin()))),pNodes.size(),rot_id,pMaxDepth,pMaxDepth, numPer);
    }
  
  SFC::seqSort::makeVectorUnique(pNodes,true);
}


template <typename T>
void MSD_Sort(std::vector<T> &pNodes, unsigned int maxDepth, int NUM_PER_OCT = 1)
{
  MSD_Sort(&(*(pNodes.begin())), pNodes.size(), 0, maxDepth, maxDepth);
}

/*

template<typename T>
void SFC_3D_msd_sort(T *pNodes, unsigned int n, unsigned int rot_id, unsigned int pMaxDepth, int numPer = 1)  
{

  assert(numPer >= 1);

  register unsigned int cnum;
  register unsigned int cnum_prev = 0;
//const int NUM_CHILDREN = 8;
  unsigned int rotation = 0;
  DendroIntL count[NUM_CHILDREN+1] = {};
  --pMaxDepth;
  
  for(DendroIntL i = 0; i < n; ++i) 
    {
      cnum = (((pNodes[i].getZ() >> pMaxDepth) & 1u) << 2u) | (((pNodes[i].getY() >> pMaxDepth) & 1u) << 1u) | ((pNodes[i].getX() >> pMaxDepth) & 1u);
      ++count[cnum+1];
    }
  
  DendroIntL loc[NUM_CHILDREN];
  T unsorted[NUM_CHILDREN];
  unsigned int live = 0;
  
  for(unsigned int i = 0; i < NUM_CHILDREN; ++i) 
    {
      cnum=rotations[ROTATION_OFFSET * rot_id+i] - '0';
      (i>0) ? cnum_prev = ((rotations[ROTATION_OFFSET * rot_id+i-1] - '0')+1) : cnum_prev = 0;
      loc[cnum] = count[cnum_prev];
      count[cnum+1] += count[cnum_prev];
      unsorted[live] = pNodes[loc[cnum]];
      if(loc[cnum] < count[cnum+1])
	++live;
    }
  --live;
  
  for(DendroIntL i = 0; i < n; ++i) 
    {
      cnum = (((unsorted[live].getZ() >> pMaxDepth) & 1u) << 2u) | (((unsorted[live].getY() >> pMaxDepth) & 1u) << 1u) | ((unsorted[live].getX() >> pMaxDepth) & 1u);
      pNodes[loc[cnum]++] = unsorted[live];
      unsorted[live] = pNodes[loc[cnum]];
      if(loc[cnum] == count[cnum + 1]) 
	--live;
    }
  
  if(pMaxDepth > 0) 
    {
      for(int i = 0; i < NUM_CHILDREN; ++i) 
	{
	  cnum = rotations[ROTATION_OFFSET*rot_id + i] - '0';
	  (i>0) ? cnum_prev = ((rotations[ROTATION_OFFSET * rot_id+i-1] - '0') + 1) : cnum_prev = 0;
	  n = count[cnum+1] - count[cnum_prev];
	  if(n > numPer) 
	    {
	      rotation = HILBERT_TABLE[NUM_CHILDREN * rot_id + cnum];
	      SFC_3D_msd_sort(pNodes + count[cnum_prev],n,rotation,pMaxDepth,numPer);
	    } 
	  else
	    {
	      // Perform the coarsening/averaging
	      barneshut::CenterOfMass(pNodes + count[cnum_prev], n);
	    }
	}
    }
}
*/
#endif
