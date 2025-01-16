/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 ****************************************************************************/

/*****************************************************************************
 * Original Author: A. Matta                 address: a.matta@surrey.ac.uk   *
 *                                                                           *
 * Creation Date  : June 2015                                                *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class holds all the online spectra needed for WAS3ABi                 *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

// NPL
#include "TWAS3ABiSpectra.h"
#include "NPOptionManager.h"
#include "NPGlobalSystemOfUnits.h"
#include "NPPhysicalConstants.h"
#ifdef NP_SYSTEM_OF_UNITS_H
using namespace NPUNITS;
#endif

// STL
#include <stdexcept>
#include <iostream>  
#include <cstdlib>
#include <string>
using namespace std;


////////////////////////////////////////////////////////////////////////////////
TWAS3ABiSpectra::TWAS3ABiSpectra(){
  SetName("WAS3ABi");
  fNumberOfDetector = 0;
  fStripFront=24;
  fStripBack=48;
}

////////////////////////////////////////////////////////////////////////////////
TWAS3ABiSpectra::TWAS3ABiSpectra(unsigned int NumberOfDetector){
  if(NPOptionManager::getInstance()->GetVerboseLevel()>0)
    cout << "************************************************" << endl
      << "TWAS3ABiSpectra : Initalising control spectra for " 
      << NumberOfDetector << " Detectors" << endl
      << "************************************************" << endl ;
  SetName("WAS3ABi");
  fNumberOfDetector = NumberOfDetector;
  fStripFront=24;
  fStripBack=48;

  InitRawSpectra();
  InitPreTreatedSpectra();
  InitPhysicsSpectra();
}

////////////////////////////////////////////////////////////////////////////////
TWAS3ABiSpectra::~TWAS3ABiSpectra(){
}

////////////////////////////////////////////////////////////////////////////////
void TWAS3ABiSpectra::InitRawSpectra(){

  static string name;
  for (unsigned int i = 0; i < fNumberOfDetector; i++) { // loop on number of detectors
    name = "WAS3ABiRaw"+NPL::itoa(i+1);
    // STR_FRONT_E_RAW
    name = "WAS3ABI"+NPL::itoa(i+1)+"_STR_FRONT_E_RAW";
    AddHisto2D(name, name, fStripFront, 1, fStripFront+1, 5000, 0, 1.5e6, "WAS3ABI/RAW/STR_FRONT_E")->Draw("colz");

    // STR_BACK_E_RAW
    name = "WAS3ABI"+NPL::itoa(i+1)+"_STR_BACK_E_RAW";
    AddHisto2D(name, name, fStripBack, 1, fStripBack+1, 5000, 0, 1.5e6, "WAS3ABI/RAW/STR_BACK_E")->Draw("colz");

    // STR_FRONT_EMAX_RAW
    name = "WAS3ABI"+NPL::itoa(i+1)+"_STR_FRONT_EMAX_RAW";
    AddHisto2D(name, name, fStripFront, 1, fStripFront+1, 5000, 0, 1.5e6, "WAS3ABI/RAW/STR_FRONT_EMAX");

    // STR_BACK_EMAX_Raw
    name = "WAS3ABI"+NPL::itoa(i+1)+"_STR_BACK_EMAX_RAW";
    AddHisto2D(name, name, fStripBack, 1, fStripBack+1, 5000, 0, 1.5e6, "WAS3ABI/RAW/STR_BACK_EMAX");


    // STR_FRONT_RAW_MULT
    name = "WAS3ABI"+NPL::itoa(i+1)+"_STR_FRONT_RAW_MULT";
    AddHisto1D(name, name, fStripFront, 1, fStripFront+1, "WAS3ABI/RAW/MULT")->Draw("");
    gPad->SetLogy();

    // STR_BACK_RAW_MULT
    name = "WAS3ABI"+NPL::itoa(i+1)+"_STR_BACK_RAW_MULT";
    AddHisto1D(name, name, fStripFront, 1, fStripFront+1, "WAS3ABI/RAW/MULT")->Draw("");
    gPad->SetLogy();


  } // end loop on number of detectors


}

////////////////////////////////////////////////////////////////////////////////
void TWAS3ABiSpectra::InitPreTreatedSpectra(){
  static string name;
  for (unsigned int i = 0; i < fNumberOfDetector; i++) { // loop on number of detectors
    // STR_FRONT_E_CAL
    name = "WAS3ABI"+NPL::itoa(i+1)+"_STR_FRONT_E_CAL";
    AddHisto2D(name, name, fStripFront, 1, fStripFront+1, 500, 0, 25, "WAS3ABI/CAL/STR_FRONT_E");

    // STR_BACK_E_CAL
    name = "WAS3ABI"+NPL::itoa(i+1)+"_STR_BACK_E_CAL";
    AddHisto2D(name, name, fStripBack, 1, fStripBack+1, 500, 0, 25, "WAS3ABI/CAL/STR_BACK_E");


    // STR_FRONT_CAL_MULT
    name = "WAS3ABI"+NPL::itoa(i+1)+"_STR_FRONT_CAL_MULT";
    AddHisto1D(name, name, fStripFront, 1, fStripFront+1, "WAS3ABI/CAL/MULT");

    // STR_BACK_CAL_MULT
    name = "WAS3ABI"+NPL::itoa(i+1)+"_STR_BACK_CAL_MULT";
    AddHisto1D(name, name, fStripFront, 1, fStripFront+1, "WAS3ABI/CAL/MULT");


    // Front-Back Energy Correlation
      name = "WAS3ABI"+NPL::itoa(i+1)+"_FB_COR";
      AddHisto2D(name, name,500,0,25,500,0,25, "WAS3ABI/CAL/FB"); 

  }  // end loop on number of detectors


}

////////////////////////////////////////////////////////////////////////////////
void TWAS3ABiSpectra::InitPhysicsSpectra(){
  static string name;
  // Kinematic Plot 
  name = "WAS3ABI_THETA_E";
  AddHisto2D(name, name,360,0,180,500,0,50,"WAS3ABI/PHY");


}



////////////////////////////////////////////////////////////////////////////////
void TWAS3ABiSpectra::FillRawSpectra(TWAS3ABiData* RawData){
  static string index;
  static string name;

  // STR_FRONT_E 
  unsigned int mysize = RawData->GetMultiplicityFront();
  double EFMAX = 0 ;
  int SFMAX = 0;
  int DFMAX = 0 ;
  for (unsigned int i = 0; i < mysize; i++) {
    index = "WAS3ABI/RAW/STR_FRONT_E/WAS3ABI"+NPL::itoa(RawData->GetFront_DetectorNbr(i))+"_STR_FRONT_E_RAW";
    if(RawData->GetFront_Energy(i) > EFMAX){
      EFMAX = RawData->GetFront_Energy(i);
      SFMAX = RawData->GetFront_StripNbr(i);
      DFMAX = RawData->GetFront_DetectorNbr(i);
    }
    
    FillSpectra(index
      ,RawData->GetFront_StripNbr(i), 
          RawData->GetFront_Energy(i));
  }
 
  if(DFMAX!=0){
    index = "WAS3ABI/RAW/STR_FRONT_EMAX/WAS3ABI"+NPL::itoa(DFMAX)+"_STR_FRONT_EMAX_RAW";
    FillSpectra(index,SFMAX, EFMAX);
  }
 
  // STR_BACK_E
  mysize = RawData->GetMultiplicityBack();
  double EBMAX = 0 ;
  int SBMAX = 0;
  int DBMAX = 0 ;
 
  for (unsigned int i = 0; i < mysize; i++) {
     index = "WAS3ABI/RAW/STR_BACK_E/WAS3ABI"+NPL::itoa( RawData->GetBack_DetectorNbr(i) )+"_STR_BACK_E_RAW";
    if(RawData->GetBack_Energy(i) > EBMAX){
      EBMAX = RawData->GetBack_Energy(i);
      SBMAX = RawData->GetBack_StripNbr(i);
      DBMAX = RawData->GetBack_DetectorNbr(i);
    }
   
    FillSpectra(index
      ,RawData->GetBack_StripNbr(i),
          RawData->GetBack_Energy(i));
  }
 
  if(DBMAX!=0){
    index = "WAS3ABI/RAW/STR_BACK_EMAX/WAS3ABI"+NPL::itoa(DBMAX)+"_STR_BACK_EMAX_RAW";
    FillSpectra(index,SBMAX, EBMAX);
  }


  // STR_FRONT MULT
  int myMULT[fNumberOfDetector];
  for( unsigned int i = 0; i < fNumberOfDetector; i++)
    myMULT[i] = 0 ; 

  for(unsigned int i = 0 ; i < RawData->GetMultiplicityFront();i++){
    myMULT[RawData->GetFront_DetectorNbr(i)-1] += 1;  
  }

  for( unsigned int i = 0; i < fNumberOfDetector; i++){
    index = "WAS3ABI/RAW/MULT/WAS3ABI"+NPL::itoa(i+1)+"_STR_FRONT_RAW_MULT";
    FillSpectra(index
      ,myMULT[i]);
  }

  // STR_BACK MULT
  for( unsigned int i = 0; i < fNumberOfDetector; i++)
    myMULT[i] = 0 ; 

  mysize = RawData->GetMultiplicityBack();
  for(unsigned int i = 0 ; i < mysize;i++){
    myMULT[RawData->GetBack_DetectorNbr(i)-1] += 1;  
  }

  for( unsigned int i = 0; i < fNumberOfDetector; i++){
    index= "WAS3ABI/RAW/MULT/WAS3ABI"+NPL::itoa(i+1)+"_STR_BACK_RAW_MULT";
    
    FillSpectra(index
      ,myMULT[i]);
  }



}

////////////////////////////////////////////////////////////////////////////////
void TWAS3ABiSpectra::FillPreTreatedSpectra(TWAS3ABiData* PreTreatedData){
  static string index;
  static string name;

  // Front-Back
  unsigned int mysizeF = PreTreatedData->GetMultiplicityFront();
  unsigned int mysizeB = PreTreatedData->GetMultiplicityBack();

  for (unsigned int i = 0; i < mysizeF; i++) {
    for (unsigned int j = 0; j < mysizeB; j++) {
      if(PreTreatedData->GetFront_DetectorNbr(i)==PreTreatedData->GetBack_DetectorNbr(j)){
        index="WAS3ABI/CAL/FB/";
        name="WAS3ABI"+NPL::itoa(PreTreatedData->GetFront_DetectorNbr(i))+"_FB_COR";

      FillSpectra(index, name
        ,PreTreatedData->GetFront_Energy(i),
                PreTreatedData->GetBack_Energy(j) );
      }
    }
  } 

  // STR_FRONT_E
  unsigned int mysize = PreTreatedData->GetMultiplicityFront();
  for (unsigned int i = 0; i < mysize; i++) {
    index = "WAS3ABI/CAL/STR_FRONT_E";
    name="WAS3ABI"+NPL::itoa(PreTreatedData->GetFront_DetectorNbr(i))+"_STR_FRONT_E_CAL";

    FillSpectra(index,name
      ,PreTreatedData->GetFront_StripNbr(i), 
          PreTreatedData->GetFront_Energy(i));
  }
  // STR_BACK_E
  mysize = PreTreatedData->GetMultiplicityBack();
  for (unsigned int i = 0; i < mysize; i++) {
   index = "WAS3ABI/CAL/STR_BACK_E";
   string name = "WAS3ABI"+NPL::itoa( PreTreatedData->GetBack_DetectorNbr(i))+"_STR_BACK_E_CAL";

    FillSpectra(index,name
      ,PreTreatedData->GetBack_StripNbr(i), 
          PreTreatedData->GetBack_Energy(i));
  }

  // STR_FRONT MULT
  int myMULT[fNumberOfDetector];
  for( unsigned int i = 0; i < fNumberOfDetector; i++)
    myMULT[i] = 0 ; 

  mysize = PreTreatedData->GetMultiplicityFront(); 
  for(unsigned int i = 0 ; i < mysize ;i++){
    myMULT[PreTreatedData->GetFront_DetectorNbr(i)-1] += 1;  
  }

  for( unsigned int i = 0; i < fNumberOfDetector; i++){
    index= "WAS3ABI/CAL/MULT";
    name= "WAS3ABI"+NPL::itoa(i+1)+"_STR_FRONT_CAL_MULT";
    FillSpectra(index,name,myMULT[i]);
  }

  // STR_BACK MULT
  for( unsigned int i = 0; i < fNumberOfDetector; i++)
    myMULT[i] = 0 ; 

  mysize = PreTreatedData->GetMultiplicityBack();
  for(unsigned int i = 0 ; i < mysize ;i++){
    myMULT[PreTreatedData->GetBack_DetectorNbr(i)-1] += 1;  
  }

  for( unsigned int i = 0; i < fNumberOfDetector; i++){
    index= "WAS3ABI/CAL/MULT";
    name = "WAS3ABI"+NPL::itoa(i+1)+"_STR_BACK_CAL_MULT";
    FillSpectra(index,name
      ,myMULT[i]);
  }




}

////////////////////////////////////////////////////////////////////////////////
void TWAS3ABiSpectra::FillPhysicsSpectra(TWAS3ABiPhysics* Physics){
  static string index="WAS3ABI/PHY";
  static string name;
  // Kine plot
  unsigned int mysize = Physics->Strip_E.size();
  for(unsigned int i = 0 ; i < mysize ; i++){
    double Theta = Physics->GetPositionOfInteraction(i).Angle(TVector3(0,0,1));
    Theta = Theta/deg;
    double Etot=Physics->Strip_E[i];

    name = "WAS3ABI_THETA_E";
    FillSpectra(index,name,Theta,Etot);
  }
}

