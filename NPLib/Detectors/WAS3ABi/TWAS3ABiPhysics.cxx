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
 *  This class hold WAS3ABi treated data                                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/
#include "TWAS3ABiPhysics.h"
using namespace WAS3ABi_LOCAL;

//   STL
#include <sstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <limits>

//   NPL
#include "RootInput.h"
#include "RootOutput.h"
#include "TAsciiFile.h"
#include "NPOptionManager.h"
#include "NPDetectorFactory.h"
//   ROOT
#include "TChain.h"
///////////////////////////////////////////////////////////////////////////

ClassImp(TWAS3ABiPhysics)
  ///////////////////////////////////////////////////////////////////////////
  TWAS3ABiPhysics::TWAS3ABiPhysics(){
    m_Rand= new TRandom3();
    EventMultiplicity   = 0 ;
    m_EventData         = new TWAS3ABiData ;
    m_PreTreatedData    = new TWAS3ABiData ;
    
    m_EventPhysics      = this ;
    m_Spectra           = NULL;
    m_NumberOfDetector = 0 ;
    m_MaximumStripMultiplicityAllowed = 10;
    m_StripEnergyMatchingSigma = 0.060 ; //keV
    m_StripEnergyMatchingNumberOfSigma = 5;

    // Threshold
    m_StripFront_E_RAW_Threshold = 0 ;
    m_StripFront_E_Threshold = 0 ;

    m_StripBack_E_RAW_Threshold = 0 ;
    m_StripBack_E_Threshold = 0 ;
    
    m_Take_E_Front=true;
    m_Take_T_Back=true;
  }

///////////////////////////////////////////////////////////////////////////
void TWAS3ABiPhysics::BuildSimplePhysicalEvent(){
  BuildPhysicalEvent();
}

///////////////////////////////////////////////////////////////////////////
void TWAS3ABiPhysics::BuildPhysicalEvent(){
  PreTreat();

  vector< TVector2 > couple = Match_Front_Back() ;
  EventMultiplicity = couple.size();
  
  /*
  for(unsigned int i = 0 ; i < EventMultiplicity ; ++i){
    cout << " iX= " << couple[i].X() << "  iY=" << couple[i].Y() << endl;
  }
  //CHECK
  */

  unsigned int size = couple.size();
  for(unsigned int i = 0 ; i < size ; ++i){
    
    int N = m_PreTreatedData->GetFront_DetectorNbr(couple[i].X()) ;

    int Front = m_PreTreatedData->GetFront_StripNbr(couple[i].X()) ;
    int Back  = m_PreTreatedData->GetBack_StripNbr(couple[i].Y()) ;

    double Front_E = m_PreTreatedData->GetFront_Energy( couple[i].X() ) ;
    double Back_E  = m_PreTreatedData->GetBack_Energy( couple[i].Y() ) ;

    double Front_T = m_PreTreatedData->GetFront_TimeCFD( couple[i].X() ) ;
    double Back_T  = m_PreTreatedData->GetBack_TimeCFD ( couple[i].Y() ) ;

    DetectorNumber.push_back(N);
    StripFront_E.push_back(Front_E);
    StripFront_T.push_back(Front_T) ;
    StripBack_E.push_back(Back_E) ;
    StripBack_T.push_back(Back_T) ;
    Strip_Front.push_back(Front) ;
    Strip_Back.push_back(Back) ;

    // Try to obtain Pixel Calibration
    static CalibrationManager* Cal = CalibrationManager::getInstance();
    static string name;
    //Store for calibration purposes
    Strip_Front_RawE.push_back(StripFront_OriginalE[couple[i].X()]);
    Strip_Back_RawE.push_back(StripBack_OriginalE[couple[i].Y()]);
    //Proceed for Pixel Calibration
    name = "WAS3ABI/D"+ NPL::itoa(N)+"_STRIP_FRONT"+ NPL::itoa(Front)+"_BACK"+ NPL::itoa(Back)+"_E";
    double Pixel_E = Cal->ApplyCalibration(name,StripFront_OriginalE[couple[i].X()] );
    if(Pixel_E != StripFront_OriginalE[couple[i].X()]){
      Strip_E.push_back(Pixel_E);
      name = "WAS3ABI/D"+ NPL::itoa(N)+"_STRIP_FRONT"+ NPL::itoa(Front)+"_BACK"+ NPL::itoa(Back)+"_DEADLAYER";
      DeadLayer.push_back(Cal->GetValue(name,0));
    }
    // Fall Back option, take the Strip Calibration
    else if(m_Take_E_Front){
      Strip_E.push_back(Front_E) ;
      name = "WAS3ABI/D"+ NPL::itoa(N)+"_STRIP_FRONT"+ NPL::itoa(Front)+"_DEADLAYER";
      DeadLayer.push_back(Cal->GetValue(name,0));
    }
    else{
      Strip_E.push_back(Back_E) ;
      name = "WAS3ABI/D"+ NPL::itoa(N)+"_STRIP_BACK"+ NPL::itoa(Back)+"_DEADLAYER";
      DeadLayer.push_back(Cal->GetValue(name,0));
    }

    if(m_Take_T_Back)
      Strip_T.push_back(Back_T) ;
    else
      Strip_T.push_back(Front_T) ;


//CHECK
   /* 
    for(unsigned int j = 0 ; j < DetectorNumber.size() ; ++j){
      cout << "DetectorNumber " <<(N)<<endl;
      cout << "StripFront_E " <<(Front_E)<<endl;
      cout << "StripFront_T " <<(Front_T) <<endl;
      cout << "StripBack_E " <<(Back_E) <<endl;
      cout << "StripBack_T " <<(Back_T) <<endl;
      cout << "Strip_Front " <<(Front) <<endl;
      cout << "Strip_Back " <<(Back) <<endl;
    }

    */
  }

  if(DetectorNumber.size()==1)
    return;
}

///////////////////////////////////////////////////////////////////////////
void TWAS3ABiPhysics::PreTreat(){
  ClearPreTreatedData();
  //   Front
  unsigned int sizeFront = m_EventData->GetMultiplicityFront();
  for(unsigned int i = 0 ; i < sizeFront ; i++){
    if( m_EventData->GetFront_Energy(i)>m_StripFront_E_RAW_Threshold && IsValidChannel("Front", m_EventData->GetFront_DetectorNbr(i), m_EventData->GetFront_StripNbr(i)) ){
      double Front_E = fStrip_Front_E(m_EventData , i);
      if( Front_E > m_StripFront_E_Threshold ){
            m_PreTreatedData->SetFront( m_EventData->GetFront_DetectorNbr(i),
            m_EventData->GetFront_StripNbr(i),
            Front_E,
            m_EventData->GetFront_TimeCFD(i),
            m_EventData->GetFront_TimeLED(i));

        StripFront_OriginalE.push_back( m_EventData->GetFront_Energy(i));
      }
    }
  }

  //  Back
  unsigned int sizeBack = m_EventData->GetMultiplicityBack() ;
  for(unsigned int i = 0 ; i < sizeBack ; i++){
    if( m_EventData->GetBack_Energy(i)>m_StripBack_E_RAW_Threshold && IsValidChannel("Back", m_EventData->GetBack_DetectorNbr(i), m_EventData->GetBack_StripNbr(i)) ){
      double Back_E = fStrip_Back_E(m_EventData , i);
      if( Back_E > m_StripBack_E_Threshold ){
        m_PreTreatedData->SetBack( m_EventData->GetBack_DetectorNbr(i),
            m_EventData->GetBack_StripNbr(i),
            Back_E,
            m_EventData->GetBack_TimeCFD(i),
            m_EventData->GetBack_TimeLED(i));

      StripBack_OriginalE.push_back( m_EventData->GetBack_Energy(i));
      }
    }
  }

  
  //m_PreTreatedData->Dump();
  //CHECK
  return;
}


///////////////////////////////////////////////////////////////////////////
int TWAS3ABiPhysics :: CheckEvent(){
  return 1 ; // Regular Event
}

///////////////////////////////////////////////////////////////////////////
vector < TVector2 > TWAS3ABiPhysics :: Match_Front_Back(){
  vector < TVector2 > ArrayOfGoodCouple ;
  // Prevent code from treating very high multiplicity Event
  // Those event are not physical anyway and that improve speed.
  if( m_PreTreatedData->GetMultiplicityFront() > m_MaximumStripMultiplicityAllowed || m_PreTreatedData->GetMultiplicityBack() > m_MaximumStripMultiplicityAllowed )
    return ArrayOfGoodCouple;

  ArrayOfGoodCouple.clear();
  unsigned int mysizeF =  m_PreTreatedData->GetMultiplicityFront();
  unsigned int mysizeB =  m_PreTreatedData->GetMultiplicityBack();
  //cout << " mysizeF/B = " << mysizeF << " " << mysizeB << endl; //CHECK

  for(unsigned int i = 0 ; i < mysizeF ; i++) {
    for(unsigned int j = 0 ; j < mysizeB ; j++){
      //cout << "Det Number F/B = " << m_PreTreatedData->GetFront_DetectorNbr(i) << " =? " << m_PreTreatedData->GetBack_DetectorNbr(j) << endl;
      //   if same detector check energy
      if ( m_PreTreatedData->GetFront_DetectorNbr(i) == m_PreTreatedData->GetBack_DetectorNbr(j) ){
        //cout << "Energy F= " << m_PreTreatedData->GetFront_Energy(i) << endl;
        //cout << "Energy B= " << m_PreTreatedData->GetBack_Energy(j) << endl;
        //cout << "abs(F-B)= " << abs( (m_PreTreatedData->GetFront_Energy(i)-m_PreTreatedData->GetBack_Energy(j))/2. ) << endl;
        //cout << "Allowed difference = " << m_StripEnergyMatchingNumberOfSigma*m_StripEnergyMatchingSigma << endl;

        //   Look if energy match
        if( abs( (m_PreTreatedData->GetFront_Energy(i)-m_PreTreatedData->GetBack_Energy(j))/2. ) < m_StripEnergyMatchingNumberOfSigma*m_StripEnergyMatchingSigma ){
          //cout << " Good Couple!! " << endl;
          ArrayOfGoodCouple.push_back ( TVector2(i,j) ) ;
        }

      }  
    }
  }

  //cout << " Good couples = " << ArrayOfGoodCouple.size()  << "  must.be.less.then  " << m_PreTreatedData->GetMultiplicityFront() << endl;

  //  Prevent to treat event with ambiguous matching beetween Front and Back
  if( ArrayOfGoodCouple.size() > m_PreTreatedData->GetMultiplicityFront() ) 
    ArrayOfGoodCouple.clear() ;
  
  //cout << " final size = " << ArrayOfGoodCouple.size() << endl;

  return ArrayOfGoodCouple;
}


////////////////////////////////////////////////////////////////////////////
bool TWAS3ABiPhysics :: IsValidChannel(const string& DetectorType, const int& telescope , const int& channel){

  if(DetectorType == "Front")
    return *(m_FrontChannelStatus[telescope-1].begin()+channel-1);

  else if(DetectorType == "Back")
    return *(m_BackChannelStatus[telescope-1].begin()+channel-1);


  else return false;
}

///////////////////////////////////////////////////////////////////////////
void TWAS3ABiPhysics::ReadAnalysisConfig(){
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/ConfigWAS3ABi.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigWAS3ABi.dat found: Default parameter loaded for Analayis " << FileName << endl;
    return;
  }
  cout << " Loading user parameter for Analysis from ConfigWAS3ABi.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigWAS3ABi.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    if (LineBuffer.compare(0, 11, "ConfigWAS3ABi") == 0) ReadingStatus = true;

    // loop on tokens and data
    while (ReadingStatus ) {

      whatToDo="";
      AnalysisConfigFile >> whatToDo;

      // Search for comment symbol (%)
      if (whatToDo.compare(0, 1, "%") == 0) {
        AnalysisConfigFile.ignore(numeric_limits<streamsize>::max(), '\n' );
      }

      else if (whatToDo=="MAX_STRIP_MULTIPLICITY") {
        AnalysisConfigFile >> DataBuffer;
        m_MaximumStripMultiplicityAllowed = atoi(DataBuffer.c_str() );
        cout << "MAXIMUN STRIP MULTIPLICITY " << m_MaximumStripMultiplicityAllowed << endl;
      }

      else if (whatToDo=="STRIP_ENERGY_MATCHING_SIGMA") {
        AnalysisConfigFile >> DataBuffer;
        m_StripEnergyMatchingSigma = atof(DataBuffer.c_str() );
        cout << "STRIP ENERGY MATCHING SIGMA " << m_StripEnergyMatchingSigma <<endl;
      }

      else if (whatToDo=="STRIP_ENERGY_MATCHING_NUMBER_OF_SIGMA") {
        AnalysisConfigFile >> DataBuffer;
        m_StripEnergyMatchingNumberOfSigma = atoi(DataBuffer.c_str() );
        cout << "STRIP ENERGY MATCHING NUMBER OF SIGMA " << m_StripEnergyMatchingNumberOfSigma << endl;
      }

      else if (whatToDo== "DISABLE_ALL") {
        AnalysisConfigFile >> DataBuffer;
        cout << whatToDo << "  " << DataBuffer << endl;
        int Detector = atoi(DataBuffer.substr(2,2).c_str());
        vector< bool > ChannelStatus;
        ChannelStatus.resize(24,false);
        m_FrontChannelStatus[Detector-1] = ChannelStatus;
        ChannelStatus.resize(48,false);
        m_BackChannelStatus[Detector-1] = ChannelStatus;
        ChannelStatus.resize(1,false);
      }

      else if (whatToDo == "DISABLE_CHANNEL") {
        AnalysisConfigFile >> DataBuffer;
        int Detector = atoi(DataBuffer.substr(2,2).c_str());
        int channel = -1;
        if (DataBuffer.find("STRF") != string::npos) {
          channel = atoi(DataBuffer.substr(8).c_str());
          *(m_FrontChannelStatus[Detector-1].begin()+channel-1) = false;
          cout << "DISABLE DETECTOR " << Detector << " STRIP FRONT " << channel << endl;
        }

        else if (DataBuffer.find("STRB")!=string::npos) {
          channel = atoi(DataBuffer.substr(8).c_str());
          *(m_BackChannelStatus[Detector-1].begin()+channel-1) = false;
          cout << "DISABLE DETECTOR " << Detector << " STRIP BACK " << channel << endl;

        }


        else cout << "Warning: detector type for WAS3ABi unknown!" << endl;

      }

      else if (whatToDo=="TAKE_E_FRONT") {
        m_Take_E_Front = true;
        cout << whatToDo << endl;
      }

      else if (whatToDo=="TAKE_E_BACK") {
        m_Take_E_Front = false;
        cout << whatToDo << endl;
      }

      else if (whatToDo=="TAKE_T_FRONT") {
        m_Take_T_Back = false;
        cout << whatToDo << endl;
      }

      else if (whatToDo=="TAKE_T_BACK") {
        m_Take_T_Back = true;
        cout << whatToDo << endl;
      }

      else if (whatToDo=="STRIP_FRONT_E_RAW_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_StripFront_E_RAW_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_StripFront_E_RAW_Threshold << endl;
      }

      else if (whatToDo=="STRIP_BACK_E_RAW_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_StripBack_E_RAW_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_StripBack_E_RAW_Threshold << endl;
      }

      else if (whatToDo=="STRIP_FRONT_E_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_StripFront_E_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_StripFront_E_Threshold << endl;
      }

      else if (whatToDo=="STRIP_BACK_E_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_StripBack_E_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_StripBack_E_Threshold << endl;
      }


      else {
        ReadingStatus = false;
      }

    }
  }

}


///////////////////////////////////////////////////////////////////////////
void TWAS3ABiPhysics::Clear(){
  EventMultiplicity = 0;

  //   Provide a Classification of Event
  EventType.clear() ;

  // Detector
  DetectorNumber.clear() ;

  //   DSSD
  Strip_E.clear() ;
  Strip_T.clear() ;
  StripFront_E.clear() ;
  StripFront_T.clear();
  StripBack_E.clear() ;
  StripBack_T.clear() ;
  Strip_Front.clear() ;
  Strip_Back.clear() ;
  StripFront_OriginalE.clear();
  StripBack_OriginalE.clear();
  DeadLayer.clear();
  Strip_Front_RawE.clear(); 
  Strip_Back_RawE.clear();
}
///////////////////////////////////////////////////////////////////////////

////   Innherited from VDetector Class   ////

///////////////////////////////////////////////////////////////////////////
void TWAS3ABiPhysics::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("WAS3ABi");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

 vector<string> tokenBOX = {"Z","ThicknessDetector1","ThicknessDetector2","ThicknessDetector3","ThicknessDetector4","ThicknessPAD1","ThicknessPAD2","ThicknessPAD3","ThicknessPAD4"};
 
  for(unsigned int i = 0 ; i < blocks.size() ; i++){

    if(blocks[i]->GetMainValue()=="BOX" && blocks[i]->HasTokenList(tokenBOX)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  WAS3ABi Box " << i+1 <<  endl;
        double Z = blocks[i]->GetDouble("Z","mm");
        double Thickness1= blocks[i]->GetDouble("ThicknessDetector1","micrometer");
        double Thickness2= blocks[i]->GetDouble("ThicknessDetector2","micrometer");
        double Thickness3= blocks[i]->GetDouble("ThicknessDetector3","micrometer");
        double Thickness4= blocks[i]->GetDouble("ThicknessDetector4","micrometer");
        AddBoxDetector(Z);
    }

    else{
      cout << "Warning: check your input file formatting " << endl;
    }
  }
  InitializeStandardParameter();
  ReadAnalysisConfig();
}
///////////////////////////////////////////////////////////////////////////
void TWAS3ABiPhysics::InitSpectra(){  
  m_Spectra = new TWAS3ABiSpectra(m_NumberOfDetector);
}

///////////////////////////////////////////////////////////////////////////
void TWAS3ABiPhysics::FillSpectra(){  
  m_Spectra -> FillRawSpectra(m_EventData);
  m_Spectra -> FillPreTreatedSpectra(m_PreTreatedData);
  m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}
///////////////////////////////////////////////////////////////////////////
void TWAS3ABiPhysics::CheckSpectra(){  
  m_Spectra->CheckSpectra();  
}
///////////////////////////////////////////////////////////////////////////
void TWAS3ABiPhysics::ClearSpectra(){  
  // To be done
}
///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TWAS3ABiPhysics::GetSpectra() {
  if(m_Spectra)
    return m_Spectra->GetMapHisto();
  else{
    map< string , TH1*> empty;
    return empty;
  }
} 

///////////////////////////////////////////////////////////////////////////
void TWAS3ABiPhysics::WriteSpectra(){
  m_Spectra->WriteSpectra();
}
///////////////////////////////////////////////////////////////////////////
void TWAS3ABiPhysics::AddParameterToCalibrationManager(){
  CalibrationManager* Cal = CalibrationManager::getInstance();

  for(int i = 0 ; i < m_NumberOfDetector ; ++i){

    for( int j = 0 ; j < 24 ; ++j){
      // Front Strip Calibration
      Cal->AddParameter("WAS3ABI", "D"+ NPL::itoa(i+1)+"_STRIP_FRONT"+ NPL::itoa(j+1)+"_E","WAS3ABI_D"+ NPL::itoa(i+1)+"_STRIP_FRONT"+ NPL::itoa(j+1)+"_E")   ;
      Cal->AddParameter("WAS3ABI", "D"+ NPL::itoa(i+1)+"_STRIP_FRONT"+ NPL::itoa(j+1)+"_DEADLAYER",
                        "WAS3ABI_D"+ NPL::itoa(i+1)+"_STRIP_FRONT"+ NPL::itoa(j+1)+"_DEADLAYER")   ;
      Cal->AddParameter("WAS3ABI", "D"+ NPL::itoa(i+1)+"_STRIP_FRONT"+ NPL::itoa(j+1)+"_T","WAS3ABI_D"+ NPL::itoa(i+1)+"_STRIP_FRONT"+ NPL::itoa(j+1)+"_T")   ;

      // Pixel Calibration
      for( int k = 0 ; k < 48 ; ++k){
        Cal->AddParameter("WAS3ABI", "D"+ NPL::itoa(i+1)+"_STRIP_FRONT"+ NPL::itoa(j+1)+"_BACK"+ NPL::itoa(k+1)+"_E","WAS3ABI_D"+ NPL::itoa(i+1)+"_STRIP_FRONT"+ NPL::itoa(j+1)+"_BACK"+ NPL::itoa(k+1)+"_E")   ;
        Cal->AddParameter("WAS3ABI", "D"+ NPL::itoa(i+1)+"_STRIP_FRONT"+ NPL::itoa(j+1)+"_BACK"+ NPL::itoa(k+1)+"_DEADLAYER",
                          "WAS3ABI_D"+ NPL::itoa(i+1)+"_STRIP_FRONT"+ NPL::itoa(j+1)+"_BACK"+ NPL::itoa(k+1)+"_DEADLAYER")   ;

      }
    }

    for( int j = 0 ; j < 48 ; ++j){
      // Back strip Calibration
      Cal->AddParameter("WAS3ABI", "D"+ NPL::itoa(i+1)+"_STRIP_BACK"+ NPL::itoa(j+1)+"_E","WAS3ABI_D"+ NPL::itoa(i+1)+"_STRIP_BACK"+ NPL::itoa(j+1)+"_E")   ;
      Cal->AddParameter("WAS3ABI", "D"+ NPL::itoa(i+1)+"_STRIP_BACK"+ NPL::itoa(j+1)+"_DEADLAYER",
                        "WAS3ABI_D"+ NPL::itoa(i+1)+"_STRIP_BACK"+ NPL::itoa(j+1)+"_DEADLAYER")   ;

      Cal->AddParameter("WAS3ABI", "D"+ NPL::itoa(i+1)+"_STRIP_BACK"+ NPL::itoa(j+1)+"_T","WAS3ABI_D"+ NPL::itoa(i+1)+"_STRIP_BACK"+ NPL::itoa(j+1)+"_T")   ;
    }

  }

  return;

}

///////////////////////////////////////////////////////////////////////////
void TWAS3ABiPhysics::InitializeRootInputRaw(){
  TChain* inputChain = RootInput::getInstance()->GetChain()   ;
  inputChain->SetBranchStatus( "WAS3ABi" , true );
  // The following line is necessary only for system were the tree is splitted
  // (older root version). The found argument silenced the Branches not found
  // warning for non splitted tree.
  if(inputChain->FindBranch("fWAS3ABi_*"))
    inputChain->SetBranchStatus( "fWAS3ABi_*",true);
  inputChain->SetBranchAddress( "WAS3ABi" , &m_EventData );

}

///////////////////////////////////////////////////////////////////////////
void TWAS3ABiPhysics::InitializeRootInputPhysics(){
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus( "EventMultiplicity" , true );
  inputChain->SetBranchStatus( "EventType" , true );
  inputChain->SetBranchStatus( "DetectorNumber" , true );
  inputChain->SetBranchStatus( "Strip_E" , true );
  inputChain->SetBranchStatus( "Strip_T" , true );
  inputChain->SetBranchStatus( "StripFront_E" , true );
  inputChain->SetBranchStatus( "StripFront_T" , true );
  inputChain->SetBranchStatus( "StripBack_E" , true );
  inputChain->SetBranchStatus( "StripBack_T" , true );
  inputChain->SetBranchStatus( "Strip_Front" , true );
  inputChain->SetBranchStatus( "Strip_Back" , true );
  inputChain->SetBranchAddress( "WAS3ABi" , &m_EventPhysics )      ;
}

///////////////////////////////////////////////////////////////////////////
void TWAS3ABiPhysics::InitializeRootOutput(){
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch( "WAS3ABi" , "TWAS3ABiPhysics" , &m_EventPhysics );
}

////////////////////////////////////////////////////////////////////////////////
/////   Specific to WAS3ABiArray   ////
void TWAS3ABiPhysics::AddBoxDetector(double Z){
  // BOX //
  double BOX_PCB_Width  = 61.10;
  double BOX_PCB_Length = 104.00;
  double BOX_PCB_Thickness = 3.4;
  double BOX_PCB_Border_LongSide = 1;
  double BOX_PCB_Border_ShortSide = 2;

  // Single stage box case (DSSD only)
  double BOX_PCB_Slot_Width1 = BOX_PCB_Thickness;
  double BOX_PCB_Slot_Border1 = 4;
  double BOX_PCB_Slot_Deepness1 = BOX_PCB_Border_ShortSide;

  // BOX Wafer
  double BOX_ActiveWafer_Width  = 48;
  double BOX_ActiveWafer_Length = 72;
  double BOX_Wafer_Width  = 52.20;
  double BOX_Wafer_Length = 76.20;  
 
  int    BOX_Wafer_Front_NumberOfStrip = 24 ;
  int    BOX_Wafer_Back_NumberOfStrip = 48 ;

  // Compute
  double BOX_LeftOver1 =  BOX_PCB_Length - BOX_PCB_Border_ShortSide - BOX_Wafer_Length - BOX_PCB_Slot_Border1 - BOX_PCB_Slot_Width1 ;
  double BOX_Exposed_Length1 = BOX_Wafer_Length + BOX_PCB_Slot_Border1 ;

  double BOX_CenterOffset1 = - 0.5 * BOX_PCB_Length+BOX_PCB_Border_ShortSide+0.5*BOX_Exposed_Length1;
  double BOX_DetectorSpacing1 = 0.5*BOX_Exposed_Length1+0.5*BOX_PCB_Slot_Width1;

  double BOX_Wafer_Width_Offset1 = -0.5*BOX_PCB_Width + BOX_PCB_Border_LongSide + 0.5*BOX_Wafer_Width;
  double BOX_Wafer_Length_Offset1 = -0.5*BOX_PCB_Length + BOX_PCB_Border_ShortSide + 0.5*BOX_Wafer_Length;

  double BOX_PCB_Slot_Position1 = 0.5*BOX_PCB_Length-BOX_LeftOver1 - 0.5*BOX_PCB_Slot_Width1;

  double StripPitchFront = BOX_ActiveWafer_Length/BOX_Wafer_Front_NumberOfStrip ; //mm
  double StripPitchBack  = BOX_ActiveWafer_Width/BOX_Wafer_Back_NumberOfStrip ; //mm
  m_BoxPitchBack = StripPitchBack;
  m_BoxPitchFront = StripPitchFront;
  

  
  double A1 = BOX_Exposed_Length1*0.5 -BOX_PCB_Slot_Border1- 0.5*StripPitchFront ; 
  double B1 = BOX_DetectorSpacing1 - 0.5*BOX_PCB_Thickness;
  double Z1 = Z - BOX_Wafer_Width*0.5 + StripPitchBack*0.5 ;


  TVector3 U; TVector3 V;TVector3 Strip_1_1;
  
  for(int i = 0 ; i < 4 ; i++){
    m_NumberOfDetector++;
    if(i==0)      {U=TVector3(1,0,0);V=TVector3(0,0,1);  Strip_1_1=TVector3( -A1 , B1  ,Z1); m_DetectorNormal.push_back(TVector3(0,-1,0));}
    else if(i==1) {U=TVector3(0,1,0);V=TVector3(0,0,1);  Strip_1_1=TVector3( -B1 , -A1 ,Z1); m_DetectorNormal.push_back(TVector3(1,0,0)) ;}
    else if(i==2) {U=TVector3(-1,0,0);V=TVector3(0,0,1); Strip_1_1=TVector3( A1  , -B1 ,Z1); m_DetectorNormal.push_back(TVector3(0,1,0)) ;}
    else if(i==3) {U=TVector3(0,-1,0);V=TVector3(0,0,1); Strip_1_1=TVector3( B1  , A1  ,Z1); m_DetectorNormal.push_back(TVector3(-1,0,0));}
    

    m_U.push_back(U);
    m_V.push_back(V);

    //   Buffer object to fill Position Array
    vector<double> lineX ; vector<double> lineY ; vector<double> lineZ ;

    vector< vector< double > >   OneBoxStripPositionX   ;
    vector< vector< double > >   OneBoxStripPositionY   ;
    vector< vector< double > >   OneBoxStripPositionZ   ;

    TVector3 StripCenter = Strip_1_1;
    for(int f = 0 ; f < BOX_Wafer_Front_NumberOfStrip ; f++){
      lineX.clear()   ;
      lineY.clear()   ;
      lineZ.clear()   ;

      for(int b = 0 ; b < BOX_Wafer_Back_NumberOfStrip ; b++){
        StripCenter = Strip_1_1 + ( StripPitchFront*f*U + StripPitchBack*b*V  );

        lineX.push_back( StripCenter.X() );
        lineY.push_back( StripCenter.Y() );
        lineZ.push_back( StripCenter.Z() );
      }

      OneBoxStripPositionX.push_back(lineX);
      OneBoxStripPositionY.push_back(lineY);
      OneBoxStripPositionZ.push_back(lineZ);
    }
    m_StripPositionX.push_back( OneBoxStripPositionX ) ;
    m_StripPositionY.push_back( OneBoxStripPositionY ) ;
    m_StripPositionZ.push_back( OneBoxStripPositionZ ) ;
  }
}
////////////////////////////////////////////////////////////////////////////////
void TWAS3ABiPhysics::AddQQQDetector( double R,double Phi,double Z){

  if(Z>0)
    m_DetectorNormal.push_back(TVector3(0,0,-1));
  else
    m_DetectorNormal.push_back(TVector3(0,0,1));

  m_U.push_back(TVector3(0,0,0));
  m_V.push_back(TVector3(0,0,0));

  double QQQ_R_Min = 9.+R;
  double QQQ_R_Max = 41.0+R;

  double QQQ_Phi_Min = 2.0*M_PI/180.  ;
  double QQQ_Phi_Max = 83.6*M_PI/180. ;

  int    QQQ_Radial_NumberOfStrip = 16 ;
  int    QQQ_Sector_NumberOfStrip = 24 ;

  double StripPitchSector = (QQQ_Phi_Max-QQQ_Phi_Min)/QQQ_Sector_NumberOfStrip ; //radial strip spacing in rad
  double StripPitchRadial = (QQQ_R_Max-QQQ_R_Min)/QQQ_Radial_NumberOfStrip  ; // ring strip spacing in mm
  m_QQQPitchBack = StripPitchSector;
  m_QQQPitchFront = StripPitchRadial;
  TVector3 Strip_1_1;

  m_NumberOfDetector++;
  Strip_1_1=TVector3(0,0,Z);

  //   Buffer object to fill Position Array
  vector<double> lineX ; vector<double> lineY ; vector<double> lineZ ;

  vector< vector< double > >   OneQQQStripPositionX   ;
  vector< vector< double > >   OneQQQStripPositionY   ;
  vector< vector< double > >   OneQQQStripPositionZ   ;

  TVector3 StripCenter = Strip_1_1;
  for(int f = 0 ; f < QQQ_Radial_NumberOfStrip ; f++){
    lineX.clear()   ;
    lineY.clear()   ;
    lineZ.clear()   ;

    for(int b = 0 ; b < QQQ_Sector_NumberOfStrip ; b++){
      StripCenter = Strip_1_1;
      StripCenter.SetY(QQQ_R_Max-f*StripPitchRadial);
      StripCenter.SetZ(Z);
      StripCenter.RotateZ(Phi+QQQ_Phi_Min+b*StripPitchSector);
      lineX.push_back( StripCenter.X() );
      lineY.push_back( StripCenter.Y() );
      lineZ.push_back( StripCenter.Z() );
    }

    OneQQQStripPositionX.push_back(lineX);
    OneQQQStripPositionY.push_back(lineY);
    OneQQQStripPositionZ.push_back(lineZ);
  }
  m_StripPositionX.push_back( OneQQQStripPositionX ) ;
  m_StripPositionY.push_back( OneQQQStripPositionY ) ;
  m_StripPositionZ.push_back( OneQQQStripPositionZ ) ;

  return;
}
////////////////////////////////////////////////////////////////////////////////
TVector3 TWAS3ABiPhysics::GetDetectorNormal( const int& i) const{
  return (m_DetectorNormal[DetectorNumber[i]-1]);

}
////////////////////////////////////////////////////////////////////////////////
TVector3 TWAS3ABiPhysics::GetPositionOfInteraction(const int& i,bool random) const{
  static TVector3 Position ;

  Position = TVector3 (  GetStripPositionX( DetectorNumber[i], Strip_Front[i], Strip_Back[i] ),
      GetStripPositionY( DetectorNumber[i] , Strip_Front[i], Strip_Back[i] ),
      GetStripPositionZ( DetectorNumber[i] , Strip_Front[i], Strip_Back[i] )) ;
  
  if(random){
    // Box Detector
    if(m_U[ DetectorNumber[i]-1].Mag()!=0){
      Position += m_V[ DetectorNumber[i]-1]*m_Rand->Uniform(-1,1)*m_BoxPitchBack*0.5;
      Position += m_U[ DetectorNumber[i]-1]*m_Rand->Uniform(-1,1)*m_BoxPitchFront*0.5;
    }

    // QQQ Detector
    else{
      Position.SetPerp( Position.Perp() + m_Rand->Uniform(-1,1)*m_QQQPitchFront*0.5);
      Position.RotateZ(m_Rand->Uniform(-1,1)*m_QQQPitchBack*0.5);
    }
  }
  
  return Position ;

}
////////////////////////////////////////////////////////////////////////////////
double TWAS3ABiPhysics::GetDeadLayer(const int& i ) const{
  return DeadLayer[i];
}
////////////////////////////////////////////////////////////////////////////////
void TWAS3ABiPhysics::InitializeStandardParameter()
{
  //   Enable all channel
  vector< bool > ChannelStatus;
  m_FrontChannelStatus.clear()    ;
  m_BackChannelStatus.clear()    ;

  ChannelStatus.resize(24,true);
  for(int i = 0 ; i < m_NumberOfDetector ; i++){
    m_FrontChannelStatus[i] = ChannelStatus;
  }

  ChannelStatus.resize(48,true);
  for(int i = 0 ; i < m_NumberOfDetector ; i++){
    m_BackChannelStatus[i] = ChannelStatus;
  }


  m_MaximumStripMultiplicityAllowed = m_NumberOfDetector   ;

  return;
}


///////////////////////////////////////////////////////////////////////////
namespace WAS3ABi_LOCAL{
  //   DSSD
  //   Front
  double fStrip_Front_E(const TWAS3ABiData* m_EventData , const int& i){
    static CalibrationManager* Cal = CalibrationManager::getInstance();
    static string name ;
    name = "WAS3ABI/D" + NPL::itoa( m_EventData->GetFront_DetectorNbr(i) ) + "_STRIP_FRONT" + NPL::itoa( m_EventData->GetFront_StripNbr(i) ) + "_E";
    return Cal->ApplyCalibration(name,m_EventData->GetFront_Energy(i) );
  }

  double fStrip_Front_T(const TWAS3ABiData* m_EventData , const int& i){
    static CalibrationManager* Cal = CalibrationManager::getInstance();
    static string name ;
    name ="WAS3ABI/D" + NPL::itoa( m_EventData->GetFront_DetectorNbr(i) ) + "_STRIP_FRONT" + NPL::itoa( m_EventData->GetFront_StripNbr(i) ) +"_T"; 

    return Cal->ApplyCalibration(name, m_EventData->GetFront_TimeCFD(i) );
  }

  //   Back
  double fStrip_Back_E(const TWAS3ABiData* m_EventData , const int& i){
    static CalibrationManager* Cal = CalibrationManager::getInstance();
    static string name ;
    name =  "WAS3ABI/D" + NPL::itoa( m_EventData->GetBack_DetectorNbr(i) ) + "_STRIP_BACK" + NPL::itoa( m_EventData->GetBack_StripNbr(i)) +"_E";

    return Cal->ApplyCalibration(name, m_EventData->GetBack_Energy(i) );
  }

  double fStrip_Back_T(const TWAS3ABiData* m_EventData , const int& i){
    static CalibrationManager* Cal = CalibrationManager::getInstance();
    static string name ;
    name = "WAS3ABI/D" + NPL::itoa( m_EventData->GetBack_DetectorNbr(i) ) + "_STRIP_BACK" + NPL::itoa( m_EventData->GetBack_StripNbr(i) ) +"_T"; 

    return Cal->ApplyCalibration(name, m_EventData->GetFront_TimeCFD(i));
  }


}

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TWAS3ABiPhysics::Construct(){
  return (NPL::VDetector*) new TWAS3ABiPhysics();
}

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
class proxy_sharc{
  public:
    proxy_sharc(){
      NPL::DetectorFactory::getInstance()->AddToken("WAS3ABi","WAS3ABi");
      NPL::DetectorFactory::getInstance()->AddDetector("WAS3ABi",TWAS3ABiPhysics::Construct);
    }
};

proxy_sharc p_sharc;
}

