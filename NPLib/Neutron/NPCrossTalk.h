#ifndef NPCROSSTALK_H
#define NPCROSSTALK_H
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 *                                                                           *
 * Original Author :  Cyril LENAIN   contact address: lenain@lpcaen.in2p3.fr *
 *                                                                           *
 * Creation Date   : November 2020                                           *
 * Last update     : November 2020                                           *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class find multi-neutrons in an neutron detection array by          *
 *    rejecting "crosstalk"                                                  *
 *  Notation:                                                                *
 *****************************************************************************/
#include <vector>
#include "TVector3.h"
#include "TGraph2DErrors.h"
#include "Math/Functor.h"
#include <map>

class TGraph;

namespace NPL{

  class CrossTalk{
    public:
      CrossTalk();
      ~CrossTalk();

    public:

      void AddHitVector(const std::vector<double>& X, const std::vector<double>& Y,const std::vector<double>& Z, const std::vector<double>& dX, const std::vector<double>& dY, const std::vector<double>& dZ, const std::vector<double>& T);

      std::vector<int> ComputeCrossTalk();
      std::vector<int> ComputeCrossTalk_HyperSphere();
      std::vector<int> GetSortedHits();
      std::vector<TGraph> GetClusters();
      TGraph2DErrors *Get3DClusters();
      std::vector<int> GetHeadClust();

    private: // private member used by
      ROOT::Math::Functor    m_func;
      const std::vector<double>* HitX;
      const std::vector<double>* HitY;
      const std::vector<double>* HitZ;
      const std::vector<double>* HitdX;
      const std::vector<double>* HitdY;
      const std::vector<double>* HitdZ;
      const std::vector<double>* HitT;
      
      std::vector<int> m_SortedID;
      std::map<unsigned int, std::vector<unsigned int>> mapOfClust;
      std::vector<int> m_HeadClust;
      std::vector<TGraph> Clusters_2D;
      TGraph2DErrors *Clusters_3D = new TGraph2DErrors();
      std::vector<int> ClustHit;
      std::vector<int> m_Neutrons;
      double coef;
      unsigned int sizeHit ;
      
  };
}

#endif
