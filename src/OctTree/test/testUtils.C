
/**
  @file testUtils.C
  @author Rahul Sampath, rahul.sampath@gmail.com
  */

#include "TreeNode.h"
#include "testUtils.h"
#include "parUtils.h"
#include "seqUtils.h"
#include <cstring>
#include "hcurvedata.h"

namespace ot {
  namespace test {

    //Assumes nodes is sorted and unique
    bool isLinear(const std::vector<ot::TreeNode > &nodes) {
      for(int i=0; i < (nodes.size() -1); i++) {
        if(nodes[i].isAncestor(nodes[i+1])){
          return false;
        }
      }
      return true;       
    }

//Assumes nodes is sorted and unique.
bool isComplete(const std::vector<ot::TreeNode >& nodes) {
  assert(!nodes.empty());

  unsigned int dim = nodes[0].getDim();
  unsigned int maxDepth = nodes[0].getMaxDepth();

  std::vector<ot::TreeNode > tmp1 = nodes;
  for (unsigned int lev = maxDepth; lev > 0; lev--) {
    //1. Select all octants at lev
    std::vector<ot::TreeNode > tmp2;
    for(unsigned int i=0; i < tmp1.size(); i++) {
      if(tmp1[i].getLevel() == lev) {
        tmp2.push_back(tmp1[i]);
      }
    }//end for i

    if( (tmp2.size()%8) != 0) {
      std::cout<<"octants in "<<lev<<" not divisible by 8"<<std::endl;
      return false;
    }

    std::vector<ot::TreeNode > tmp3;
    unsigned int tmp2Cnt = 0;
    while(tmp2Cnt < tmp2.size()) {
      ot::TreeNode currParent = tmp2[tmp2Cnt].getParent();
      for(int i = 1; i < (1 << dim); i++) {
        assert((tmp2Cnt + i) < tmp2.size()); 
        if(tmp2[tmp2Cnt + i].getParent() != currParent) {
          std::cout<<"Treenodes at lev:"<<lev<<" do not have the same parents"<<std::endl;
          return false; 
        }
      }
      tmp3.push_back(currParent);
      tmp2Cnt +=8;
    }//end while

    tmp2.clear();

    for(unsigned int i=0; i < tmp1.size(); i++) {
      if(tmp1[i].getLevel() != lev) {
        tmp3.push_back(tmp1[i]);
      }
    }//end for i

    std::sort(tmp3.begin(),tmp3.end());
    tmp1 = tmp3;
    tmp3.clear();

  }//end for lev
  return true;
}



bool isBalanced(unsigned int dim, unsigned int maxDepth, char*failFileName,
    const std::vector<ot::TreeNode > &nodes, bool incCorn, unsigned int maxLevDiff) {
  TreeNode root (dim, maxDepth);
  return	isBalancedInternal(dim,maxDepth,failFileName,nodes,root, incCorn, maxLevDiff);
}//end function

bool isBalancedInternal(unsigned int dim, unsigned int maxDepth,char*failFileName,

  const  std::vector<ot::TreeNode > & nodes, ot::TreeNode holder,
  bool incCorn, unsigned int maxLevDiff) {
  bool yesBalanced = true;
  std::vector<TreeNode> failedCorners;
  std::vector<TreeNode> failedEdges;
  std::vector<TreeNode> failedFaces;
  TreeNode root (dim, maxDepth);
  unsigned int retIdx;
  double nxz = (double) nodes.size();
  for (unsigned int i=0; i<nodes.size(); i++) {
    /*
       if (!(i%100)) {
       printf("%4.2f\r", 100.0*((double)i)/nxz);
       }
       */

    // cout<<RED<<"R"<<NRM<<endl;
    {
      TreeNode it = nodes[i].getRight();

      if(holder.isAncestor(it)){

        bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
        if(found) {
          unsigned int retLev = nodes[retIdx].getLevel();
          unsigned int myLev = nodes[i].getLevel();					
          if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
            found = false;
          }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
            //Found a very big Neighbour
            found = false;
          }						
        }
        if(!found) {
          yesBalanced = false;
          std::cout<<nodes[i]<<": (R) ->"<<nodes[retIdx]<<std::endl;
          failedFaces.push_back(nodes[i]);					
          //break;
        }
      }

    }

    // cout<<RED<<"L"<<NRM<<endl;
    {
      TreeNode it = nodes[i].getLeft();

      if(holder.isAncestor(it)){

        bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
        if(found) {
          unsigned int retLev = nodes[retIdx].getLevel();
          unsigned int myLev = nodes[i].getLevel();					
          if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
            found = false;
          }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
            //Found a very big Neighbour
            found = false;
          }						
        }
        if(!found) {
          yesBalanced = false;
          std::cout<<nodes[i]<<": (L) ->"<<nodes[retIdx]<<std::endl;
          failedFaces.push_back(nodes[i]);					
          //break;
        }
      }

    }


    // cout<<RED<<"T"<<NRM<<endl;
    {
      TreeNode it = nodes[i].getTop();

      if(holder.isAncestor(it)){

        bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
        if(found) {
          unsigned int retLev = nodes[retIdx].getLevel();
          unsigned int myLev = nodes[i].getLevel();					
          if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
            found = false;
          }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
            //Found a very big Neighbour
            found = false;
          }						
        }
        if(!found) {
          yesBalanced = false;
          std::cout<<nodes[i]<<": (T) ->"<<nodes[retIdx]<<std::endl;
          failedFaces.push_back(nodes[i]);					
          //break;
        }
      }
    }

    // cout<<RED<<"Bo"<<NRM<<endl;
    {
      TreeNode it = nodes[i].getBottom();


      if(holder.isAncestor(it)){

        bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
        if(found) {
          unsigned int retLev = nodes[retIdx].getLevel();
          unsigned int myLev = nodes[i].getLevel();					
          if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
            found = false;
          }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
            //Found a very big Neighbour
            found = false;
          }						
        }
        if(!found) {
          yesBalanced = false;
          std::cout<<nodes[i]<<": (Bo) ->"<<nodes[retIdx]<<std::endl;
          failedFaces.push_back(nodes[i]);					
         // break;
        }
      }

    }

    // cout<<RED<<"F"<<NRM<<endl;
    {
      TreeNode it = nodes[i].getFront();


      if(holder.isAncestor(it)){

        bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
        if(found) {
          unsigned int retLev = nodes[retIdx].getLevel();
          unsigned int myLev = nodes[i].getLevel();					
          if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
            found = false;
          }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
            //Found a very big Neighbour
            found = false;
          }						
        }
        if(!found) {
          yesBalanced = false;
          std::cout<<nodes[i]<<": (F) ->"<<nodes[retIdx]<<std::endl;
          failedFaces.push_back(nodes[i]);					
         // break;
        }
      }

    }

    // cout<<RED<<"Bk"<<NRM<<endl;
    {
      TreeNode it = nodes[i].getBack();


      if(holder.isAncestor(it)){

        bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
        if(found) {
          unsigned int retLev = nodes[retIdx].getLevel();
          unsigned int myLev = nodes[i].getLevel();					
          if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
            found = false;
          }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
            //Found a very big Neighbour
            found = false;
          }						
        }
        if(!found) {
          yesBalanced = false;
          std::cout<<nodes[i]<<": (Bk) ->"<<nodes[retIdx]<<std::endl;
          failedFaces.push_back(nodes[i]);					
          //break;
        }
      }

    }

    if(dim == 3 || incCorn) {
      // cout<<RED<<"TR"<<NRM<<endl;
      {
        TreeNode it = nodes[i].getTopRight();


        if(holder.isAncestor(it)){

          bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
          if(found) {
            unsigned int retLev = nodes[retIdx].getLevel();
            unsigned int myLev = nodes[i].getLevel();					
            if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
              found = false;
            }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
              //Found a very big Neighbour
              found = false;
            }						
          }
          if(!found) {
            yesBalanced = false;
            std::cout<<nodes[i]<<": (TR) ->"<<nodes[retIdx]<<std::endl;
            failedEdges.push_back(nodes[i]);					
            //break;
          }
        }

      }

      // cout<<RED<<"TL"<<NRM<<endl;
      {
        TreeNode it = nodes[i].getTopLeft();


        if(holder.isAncestor(it)){

          bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
          if(found) {
            unsigned int retLev = nodes[retIdx].getLevel();
            unsigned int myLev = nodes[i].getLevel();					
            if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
              found = false;
            }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
              //Found a very big Neighbour
              found = false;
            }						
          }
          if(!found) {
            yesBalanced = false;
            std::cout<<nodes[i]<<": (TL) ->"<<nodes[retIdx]<<std::endl;
            failedEdges.push_back(nodes[i]);					
           // break;
          }
        }

      }

      // cout<<RED<<"BoL"<<NRM<<endl;
      {
        TreeNode it = nodes[i].getBottomLeft();


        if(holder.isAncestor(it)){

          bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
          if(found) {
            unsigned int retLev = nodes[retIdx].getLevel();
            unsigned int myLev = nodes[i].getLevel();					
            if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
              found = false;
            }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
              //Found a very big Neighbour
              found = false;
            }						
          }
          if(!found) {
            yesBalanced = false;
            std::cout<<nodes[i]<<": (BoL) ->"<<nodes[retIdx]<<std::endl;
            failedEdges.push_back(nodes[i]);					
            //break;
          }
        }

      }

      // cout<<RED<<"BoR"<<NRM<<endl;
      {
        TreeNode it = nodes[i].getBottomRight();


        if(holder.isAncestor(it)){

          bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
          if(found) {
            unsigned int retLev = nodes[retIdx].getLevel();
            unsigned int myLev = nodes[i].getLevel();					
            if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
              found = false;
            }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
              //Found a very big Neighbour
              found = false;
            }						
          }
          if(!found) {
            yesBalanced = false;
            std::cout<<nodes[i]<<": (BoR) ->"<<nodes[retIdx]<<std::endl;
            failedEdges.push_back(nodes[i]);					
           // break;
          }
        }

      }

      // cout<<RED<<"RBk"<<NRM<<endl;
      {
        TreeNode it = nodes[i].getRightBack();


        if(holder.isAncestor(it)){

          bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
          if(found) {
            unsigned int retLev = nodes[retIdx].getLevel();
            unsigned int myLev = nodes[i].getLevel();					
            if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
              found = false;
            }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
              //Found a very big Neighbour
              found = false;
            }						
          }
          if(!found) {
            yesBalanced = false;
            std::cout<<nodes[i]<<": (RBk) ->"<<nodes[retIdx]<<std::endl;
            failedEdges.push_back(nodes[i]);					
            //break;
          }
        }

      }

      // cout<<RED<<"TBk"<<NRM<<endl;
      {
        TreeNode it = nodes[i].getTopBack();


        if(holder.isAncestor(it)){

          bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
          if(found) {
            unsigned int retLev = nodes[retIdx].getLevel();
            unsigned int myLev = nodes[i].getLevel();					
            if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
              found = false;
            }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
              //Found a very big Neighbour
              found = false;
            }						
          }
          if(!found) {
            yesBalanced = false;
            std::cout<<nodes[i]<<": (TBk) ->"<<nodes[retIdx]<<std::endl;
            failedEdges.push_back(nodes[i]);					
           // break;
          }
        }

      }

      // cout<<RED<<"LF"<<NRM<<endl;
      {
        TreeNode it = nodes[i].getLeftFront();


        if(holder.isAncestor(it)){

          bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
          if(found) {
            unsigned int retLev = nodes[retIdx].getLevel();
            unsigned int myLev = nodes[i].getLevel();					
            if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
              found = false;
            }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
              //Found a very big Neighbour
              found = false;
            }						
          }
          if(!found) {
            yesBalanced = false;
            std::cout<<nodes[i]<<": (LF) ->"<<nodes[retIdx]<<std::endl;
            failedEdges.push_back(nodes[i]);					
            //break;
          }
        }

      }

      // cout<<RED<<"LBk"<<NRM<<endl;
      {
        TreeNode it = nodes[i].getLeftBack();


        if(holder.isAncestor(it)){

          bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
          if(found) {
            unsigned int retLev = nodes[retIdx].getLevel();
            unsigned int myLev = nodes[i].getLevel();					
            if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
              found = false;
            }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
              //Found a very big Neighbour
              found = false;
            }						
          }
          if(!found) {
            yesBalanced = false;
            std::cout<<nodes[i]<<": (LBk) ->"<<nodes[retIdx]<<std::endl;
            failedEdges.push_back(nodes[i]);					
            //break;
          }
        }

      }

      // cout<<RED<<"RF"<<NRM<<endl;
      {
        TreeNode it = nodes[i].getRightFront();


        if(holder.isAncestor(it)){

          bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
          if(found) {
            unsigned int retLev = nodes[retIdx].getLevel();
            unsigned int myLev = nodes[i].getLevel();					
            if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
              found = false;
            }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
              //Found a very big Neighbour
              found = false;
            }						
          }
          if(!found) {
            yesBalanced = false;
            std::cout<<nodes[i]<<": (RF) ->"<<nodes[retIdx]<<std::endl;
            failedEdges.push_back(nodes[i]);					
            //break;
          }
        }

      }

      // cout<<RED<<"TF"<<NRM<<endl;
      {
        TreeNode it = nodes[i].getTopFront();


        if(holder.isAncestor(it)){

          bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
          if(found) {
            unsigned int retLev = nodes[retIdx].getLevel();
            unsigned int myLev = nodes[i].getLevel();					
            if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
              found = false;
            }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
              //Found a very big Neighbour
              found = false;
            }						
          }
          if(!found) {
            yesBalanced = false;
            std::cout<<nodes[i]<<": (TF) ->"<<nodes[retIdx]<<std::endl;
            failedEdges.push_back(nodes[i]);					
            //break;
          }
        }

      }

      // cout<<RED<<"BoBk"<<NRM<<endl;
      {
        TreeNode it = nodes[i].getBottomBack();


        if(holder.isAncestor(it)){

          bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
          if(found) {
            unsigned int retLev = nodes[retIdx].getLevel();
            unsigned int myLev = nodes[i].getLevel();					
            if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
              found = false;
            }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
              //Found a very big Neighbour
              found = false;
            }						
          }
          if(!found) {
            yesBalanced = false;
            std::cout<<nodes[i]<<": (BoBk) ->"<<nodes[retIdx]<<std::endl;
            failedEdges.push_back(nodes[i]);					
            //break;
          }
        }

      }

      // cout<<RED<<"BoF"<<NRM<<endl;
      {
        TreeNode it = nodes[i].getBottomFront();


        if(holder.isAncestor(it)){

          bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
          if(found) {
            unsigned int retLev = nodes[retIdx].getLevel();
            unsigned int myLev = nodes[i].getLevel();					
            if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
              found = false;
            }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
              //Found a very big Neighbour
              found = false;
            }						
          }
          if(!found) {
            yesBalanced = false;
            std::cout<<nodes[i]<<": (BoF) ->"<<nodes[retIdx]<<std::endl;
            failedEdges.push_back(nodes[i]);					
            //break;
          }
        }

      }

    }//end if dim=3 or incCorn

    if(incCorn) {
      // cout<<RED<<"TRBk"<<NRM<<endl;
      {
        TreeNode it = nodes[i].getTopRightBack();


        if(holder.isAncestor(it)){

          bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
          if(found) {
            unsigned int retLev = nodes[retIdx].getLevel();
            unsigned int myLev = nodes[i].getLevel();					
            if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
              found = false;
            }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
              //Found a very big Neighbour
              found = false;
            }						
          }
          if(!found) {
            yesBalanced = false;
            std::cout<<nodes[i]<<": (TRBk) ->"<<nodes[retIdx]<<std::endl;
            failedCorners.push_back(nodes[i]);					
            //break;
          }
        }

      }

      // cout<<RED<<"TLF"<<NRM<<endl;
      {
        TreeNode it = nodes[i].getTopLeftFront();


        if(holder.isAncestor(it)){

          bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
          if(found) {
            unsigned int retLev = nodes[retIdx].getLevel();
            unsigned int myLev = nodes[i].getLevel();					
            if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
              found = false;
            }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
              //Found a very big Neighbour
              found = false;
            }						
          }
          if(!found) {
            yesBalanced = false;
            std::cout<<nodes[i]<<": (TLF) ->"<<nodes[retIdx]<<std::endl;
            failedCorners.push_back(nodes[i]);					
            //break;
          }
        }

      }

      // cout<<RED<<"BoRBk"<<NRM<<endl;
      {
        TreeNode it = nodes[i].getBottomRightBack();


        if(holder.isAncestor(it)){

          bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
          if(found) {
            unsigned int retLev = nodes[retIdx].getLevel();
            unsigned int myLev = nodes[i].getLevel();					
            if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
              found = false;
            }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
              //Found a very big Neighbour
              found = false;
            }						
          }
          if(!found) {
            yesBalanced = false;
            std::cout<<nodes[i]<<": (BoRBk) ->"<<nodes[retIdx]<<std::endl;
            failedCorners.push_back(nodes[i]);					
            //break;
          }
        }

      }

      // cout<<RED<<"TLBk"<<NRM<<endl;
      {
        TreeNode it = nodes[i].getTopLeftBack();


        if(holder.isAncestor(it)){

          bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
          if(found) {
            unsigned int retLev = nodes[retIdx].getLevel();
            unsigned int myLev = nodes[i].getLevel();					
            if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
              found = false;
            }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
              //Found a very big Neighbour
              found = false;
            }						
          }
          if(!found) {
            yesBalanced = false;
            std::cout<<nodes[i]<<": (TLBk) ->"<<nodes[retIdx]<<std::endl;
            failedCorners.push_back(nodes[i]);					
           // break;
          }
        }

      }

      // cout<<RED<<"BoRF"<<NRM<<endl;
      {
        TreeNode it = nodes[i].getBottomRightFront();


        if(holder.isAncestor(it)){

          bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
          if(found) {
            unsigned int retLev = nodes[retIdx].getLevel();
            unsigned int myLev = nodes[i].getLevel();					
            if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
              found = false;
            }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
              //Found a very big Neighbour
              found = false;
            }						
          }
          if(!found) {
            yesBalanced = false;
            std::cout<<nodes[i]<<": (BoRF) ->"<<nodes[retIdx]<<std::endl;
            failedCorners.push_back(nodes[i]);					
            //break;
          }
        }

      }

      // cout<<RED<<"TRF"<<NRM<<endl;
      {
        TreeNode it = nodes[i].getTopRightFront();


        if(holder.isAncestor(it)){

          bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
          if(found) {
            unsigned int retLev = nodes[retIdx].getLevel();
            unsigned int myLev = nodes[i].getLevel();					
            if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
              found = false;
            }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
              //Found a very big Neighbour
              found = false;
            }						
          }
          if(!found) {
            yesBalanced = false;
            std::cout<<nodes[i]<<": (TRF) ->"<<nodes[retIdx]<<std::endl;
            failedCorners.push_back(nodes[i]);					
            //break;
          }
        }

      }


      // cout<<RED<<"BoLBk"<<NRM<<endl;
      {
        TreeNode it = nodes[i].getBottomLeftBack();


        if(holder.isAncestor(it)){

          bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
          if(found) {
            unsigned int retLev = nodes[retIdx].getLevel();
            unsigned int myLev = nodes[i].getLevel();					
            if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
              found = false;
            }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
              //Found a very big Neighbour
              found = false;
            }						
          }
          if(!found) {
            yesBalanced = false;
            std::cout<<nodes[i]<<": (BoLBk) ->"<<nodes[retIdx]<<std::endl;
            failedCorners.push_back(nodes[i]);					
            //break;
          }
        }


      }


      // cout<<RED<<"BoLF"<<NRM<<endl;
      {
        TreeNode it = nodes[i].getBottomLeftFront();


        if(holder.isAncestor(it)){


          bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
          if(found) {
            unsigned int retLev = nodes[retIdx].getLevel();
            unsigned int myLev = nodes[i].getLevel();					
            if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
              found = false;
            }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
              //Found a very big Neighbour
              found = false;
            }						
          }
          if(!found) {
            yesBalanced = false;
            std::cout<<nodes[i]<<": (BoLF) ->"<<nodes[retIdx]<<std::endl;
            failedCorners.push_back(nodes[i]);					
           // break;
          }
        }

      }

    }//end if incCorn

  }//end for i

  seq::makeVectorUnique(failedFaces,false);
  seq::makeVectorUnique(failedEdges,false);
  seq::makeVectorUnique(failedCorners,false);

  treeNodesTovtk(failedFaces,0,"failedFaces");
  treeNodesTovtk(failedEdges,0,"failedEdges");
  treeNodesTovtk(failedCorners,0,"failedCorners");


  if(!yesBalanced)
  {
    std::cout << BLU << "===============================================" << NRM << std::endl;
    std::cout << RED " Balance test failed. " NRM << std::endl;
    std::cout << YLW << "Failed Face Cases:"<<failedFaces.size()<<NRM<<std::endl;
    std::cout << YLW << "Failed Edges Cases:"<<failedEdges.size()<<NRM<<std::endl;
    std::cout << YLW << "Failed Corners Cases:"<<failedCorners.size()<<NRM<<std::endl;
    std::cout << BLU << "===============================================" << NRM << std::endl;

  }


  char failCornerFileName[100];
  char failEdgeFileName[100];
  char failFaceFileName[100];
  strcpy(failFaceFileName,failFileName );
  strcpy(failEdgeFileName,failFileName);
  strcpy(failCornerFileName,failFileName);
  strcat(failFaceFileName,"_Faces.ot\0");	
  strcat(failEdgeFileName,"_Edges.ot\0");	
  strcat(failCornerFileName,"_Corners.ot\0");	
  writeNodesToFile(failCornerFileName,failedCorners);
  writeNodesToFile(failEdgeFileName,failedEdges);
  writeNodesToFile(failFaceFileName,failedFaces);
  return yesBalanced;
}//end function

}//end namespace
}//end namespace

namespace oda {
    namespace test {

        bool odaTest(std::vector<ot::TreeNode> &in, ot::DA &da, MPI_Comm comm) {

          int rank, size;
          MPI_Comm_rank(comm, &rank);
          MPI_Comm_size(comm, &size);
          std::vector<ot::TreeNode> DNodes;
          std::vector<ot::TreeNode> DKeys;

          if (!rank)
            std::cout << GRN << "============ODA TEST START==============" << NRM << std::endl;

          if (da.iAmActive()) {
            std::vector<unsigned int> nodeIdx;
            unsigned int *nodeIdx_ptr;
            ot::DA da_cpy(da);

            for (da.init<ot::DA_FLAGS::ALL>(); da.curr() < (da.end<ot::DA_FLAGS::ALL>());da.next<ot::DA_FLAGS::ALL>()) {
              nodeIdx.resize(8);
              nodeIdx_ptr = &(*nodeIdx.begin());
              da.getNodeIndices(nodeIdx_ptr);
              Point currentPoint;
              currentPoint=da.getCurrentOffset();
              unsigned int currentLev=da.getLevel(da.curr());
              std::vector<ot::TreeNode> child;

              ot::TreeNode currentNode(1,currentPoint.xint(),currentPoint.yint(),currentPoint.zint(),currentLev,da.getDimension(),da.getMaxDepth());

              if(!da.isHanging(da.curr())) {
                currentNode.getParent().addChildrenMorton(child);
                for (int j = 0; j < child.size(); j++) {

                  for(da_cpy.init<ot::DA_FLAGS::ALL>();da_cpy.curr()<da_cpy.end<ot::DA_FLAGS::ALL>();da_cpy.next<ot::DA_FLAGS::ALL>())
                  {
                    if(da_cpy.curr()==nodeIdx[j])
                    {
                      assert(da_cpy.getCurrentOffset()==currentNode.getAnchor());
                      break;
                    }
                  }

                }

                nodeIdx.clear();
                child.clear();

              }

            }


          }
        }

        bool odaLoopTestAll(ot::DA& da, MPI_Comm comm)
        {
          std::vector<ot::TreeNode> allNodes;
          ot::TreeNode tmp;
          int rank,size;
          MPI_Comm_rank(comm,&rank);
          MPI_Comm_size(comm,&size);

          for(da.init<ot::DA_FLAGS::ALL>();da.curr()<da.end<ot::DA_FLAGS::ALL>();da.next<ot::DA_FLAGS::ALL>())
          {

            Point p=da.getCurrentOffset();

            tmp=ot::TreeNode(1,p.xint(),p.yint(),p.zint(),da.getLevel(da.curr()),da.getDimension(),da.getMaxDepth());
            allNodes.push_back(tmp);
//            if(!rank)
//              treeNodesTovtk(allNodes,da.curr(),"oda_loop_all0");
//            else
//              treeNodesTovtk(allNodes,da.curr(),"oda_loop_all1");
//
//            allNodes.clear();

          }

          treeNodesTovtk(allNodes,rank,"oda_loop_all");
          //assert(seq::test::isUniqueAndSorted(allNodes));

          return true;
          // this won't be globally sorted since we have the ghost octants.

        }



        bool odaLoopTestWritable(ot::DA & da, MPI_Comm comm)
        {
          std::vector<ot::TreeNode> allWNodes;
          ot::TreeNode tmp;
          int rank,size;
          MPI_Comm_rank(comm,&rank);
          MPI_Comm_size(comm,&size);

          for(da.init<ot::DA_FLAGS::WRITABLE>();da.curr()<0.25*da.end<ot::DA_FLAGS::WRITABLE>();da.next<ot::DA_FLAGS::WRITABLE>())
          {

            Point p=da.getCurrentOffset();

            tmp=ot::TreeNode(1,p.xint(),p.yint(),p.zint(),da.getLevel(da.curr()),da.getDimension(),da.getMaxDepth());
            allWNodes.push_back(tmp);

//            if(!rank)
//              treeNodesTovtk(allWNodes,da.curr(),"oda_loop_W0");
//            else
//              treeNodesTovtk(allWNodes,da.curr(),"oda_loop_W1");
//
//            allWNodes.clear();


          }

          treeNodesTovtk(allWNodes,rank,"oda_loop_writable");
          assert(seq::test::isUniqueAndSorted(allWNodes));

          return true;
          // this won't be globally sorted since we have the ghost octants.

        }



    }
}