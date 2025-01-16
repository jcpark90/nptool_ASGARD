/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : November 2012                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold WAS3ABi Raw data                                           *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

#include "TWAS3ABiData.h"

ClassImp(TWAS3ABiData)

/////////////////////////
TWAS3ABiData::TWAS3ABiData(){
}

/////////////////////////
TWAS3ABiData::~TWAS3ABiData(){
}

/////////////////////////
void TWAS3ABiData::Clear(){
  fWAS3ABi_StripFront_DetectorNbr.clear();
  fWAS3ABi_StripFront_StripNbr.clear();
  fWAS3ABi_StripFront_Energy.clear();
  fWAS3ABi_StripFront_TimeCFD.clear();
  fWAS3ABi_StripFront_TimeLED.clear();

  fWAS3ABi_StripBack_DetectorNbr.clear();
  fWAS3ABi_StripBack_StripNbr.clear();
  fWAS3ABi_StripBack_Energy.clear();
  fWAS3ABi_StripBack_TimeCFD.clear();
  fWAS3ABi_StripBack_TimeLED.clear();
  

}

/////////////////////////
void TWAS3ABiData::Dump() const{
   
  // Front
  cout << "WAS3ABi Strip Front Mult = " << fWAS3ABi_StripFront_DetectorNbr.size() << endl;
  for (UShort_t i = 0; i < fWAS3ABi_StripFront_DetectorNbr.size(); i++){
    cout << "DetNbr (Front): " << fWAS3ABi_StripFront_DetectorNbr[i]
         << "   Strip: " << fWAS3ABi_StripFront_StripNbr[i]
         << "   Energy: " << fWAS3ABi_StripFront_Energy[i]
         << "   Time CFD: " << fWAS3ABi_StripFront_TimeCFD[i]
         << "   Time LED: " << fWAS3ABi_StripFront_TimeLED[i] << endl;
  }

  // Back  
  cout << "WAS3ABi Strip Back Mult  = " << fWAS3ABi_StripBack_DetectorNbr.size() << endl;
  for (UShort_t i = 0; i < fWAS3ABi_StripBack_DetectorNbr.size(); i++){
    cout << "DetNbr (Back): " << fWAS3ABi_StripBack_DetectorNbr[i]
    << "   Strip: " << fWAS3ABi_StripBack_StripNbr[i]
    << "   Energy: " << fWAS3ABi_StripBack_Energy[i]
    << "   Time CFD: " << fWAS3ABi_StripBack_TimeCFD[i]
    << "   Time LED: " << fWAS3ABi_StripBack_TimeLED[i] << endl;
  }
  
}
