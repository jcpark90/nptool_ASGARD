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
 *  This class describe the WAS3ABi Silicon array                              *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

// C++ headers
#include <sstream>
#include <cmath>
#include <limits>
//G4 Geometry object
#include "G4Box.hh"
#include "G4Tubs.hh"

//G4 sensitive
#include "G4SDManager.hh"

//G4 various object
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4Colour.hh"
#include "G4PVDivision.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
// NPS
#include "WAS3ABi.hh"
#include "DSSDScorers.hh"
#include "InteractionScorers.hh"
#include "MaterialManager.hh"
#include "NPSDetectorFactory.hh"
// NPL
#include "NPOptionManager.h"

#include "RootOutput.h"
using namespace WAS3ABI;

// CLHEP header
#include "CLHEP/Random/RandGauss.h"

using namespace std;
using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// WAS3ABi Specific Method
WAS3ABi::WAS3ABi(){
  InitializeMaterial();
  m_Event = new TWAS3ABiData();
  // Dark Grey
  SiliconVisAtt = new G4VisAttributes(G4Colour(0.3, 0.3, 0.3)) ;
  // Green
  PCBVisAtt = new G4VisAttributes(G4Colour(0.2, 0.5, 0.2)) ;
  // Light Grey
  FrameVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5)) ;

  m_DSSSDScorer = 0 ;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
WAS3ABi::~WAS3ABi(){
  //delete m_DSSSDScorer;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void WAS3ABi::AddDSSSDDetector(G4double Z, unsigned int Nlayers, G4double GapDistance, G4double Thickness){

  m_Type.push_back(true);
  m_ZDSSSD.push_back(Z);
  m_NumLayersDSSSD.push_back(Nlayers);
  m_GapDistanceDSSSD.push_back(GapDistance);
  m_ThicknessDSSSD.push_back(Thickness);


}
void WAS3ABi::AddPassiveStopper(G4double Z, G4String Material, G4double Thickness, G4double Width, G4double Height){

  m_Type.push_back(false);
  m_ZPassive.push_back(Z);  
  m_MaterialPassive.push_back(Material);
  m_ThicknessPassive.push_back(Thickness);
  m_WidthPassive.push_back(Width);
  m_HeightPassive.push_back(Height);

}

void WAS3ABi::AddPassiveCircleStopper(G4double Z, G4String Material, G4double Thickness, G4double R){

  m_Type.push_back(false);
  m_ZPassive.push_back(Z);  
  m_MaterialPassive.push_back(Material);
  m_ThicknessPassive.push_back(Thickness);
  m_RadiPassive.push_back(R);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class
// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void WAS3ABi::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("WAS3ABi");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> tokenDSSSD = {"Z","NumLayers", "GapDistance", "Thickness"};
  vector<string> tokenPassive = {"Z", "Material", "Thickness", "Width", "Height"};
  vector<string> tokenPassiveCircle = {"Z", "Material", "Thickness", "R"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){

    if(blocks[i]->GetMainValue()=="DSSSD" && blocks[i]->HasTokenList(tokenDSSSD)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  WAS3ABi DSSSD " << i+1 <<  endl;
      double Z = blocks[i]->GetDouble("Z","mm");
      unsigned int NumLayers = blocks[i]->GetInt("NumLayers");
      double GapDistance = blocks[i]->GetDouble("GapDistance","mm");
      double Thickness= blocks[i]->GetDouble("Thickness","micrometer");
      AddDSSSDDetector(Z, NumLayers, GapDistance, Thickness);
    }
    else if(blocks[i]->GetMainValue()=="Passive" && blocks[i]->HasTokenList(tokenPassive)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Passive stopper " << i+1 <<  endl;
      double Z = blocks[i]->GetDouble("Z","mm");
      string Material = blocks[i]->GetString("Material");
      double Thickness= blocks[i]->GetDouble("Thickness","micrometer");
      double Width= blocks[i]->GetDouble("Width","mm");
      double Height= blocks[i]->GetDouble("Height","mm");
      
      AddPassiveStopper(Z, Material, Thickness, Width, Height);
      
    }
    else if(blocks[i]->GetMainValue()=="PassiveCircle" && blocks[i]->HasTokenList(tokenPassiveCircle)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  PassiveCircle stopper " << i+1 <<  endl;
      double Z = blocks[i]->GetDouble("Z","mm");
      string Material = blocks[i]->GetString("Material");
      double Thickness= blocks[i]->GetDouble("Thickness","micrometer");
      double Radi = blocks[i]->GetDouble("R","mm");
      
      AddPassiveCircleStopper(Z, Material, Thickness, Radi);
      
    }
    else{
      cout << "Warning: check your input file formatting " << endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Construct detector and inialise sensitive part.
// Called After DetecorConstruction::AddDetector Method
void WAS3ABi::ConstructDetector(G4LogicalVolume* world){
  ConstructDSSSDDetector(world);
  ConstructPassiveStopper(world);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void WAS3ABi::ConstructDSSSDDetector(G4LogicalVolume* world){
  for(unsigned int i = 0 ; i < m_ZDSSSD.size() ; i++){
    for (unsigned int j = 0 ; j < m_NumLayersDSSSD[i]; j++) {

      int DetNbr = j;
      // create the DSSSD 
      // Make the a single detector geometry
      G4Box*  PCBFull = new G4Box("PCBFull"  ,
				  DSSSD_PCB_Length/2.,
				  DSSSD_PCB_Width/2.,
				  DSSSD_PCB_Thickness/2.);

      G4Box*  WaferShape = new G4Box("WaferShape",
				     DSSSD_Wafer_Length/2.,
				     DSSSD_Wafer_Width/2.,
				     DSSSD_PCB_Thickness/2.+0.1*mm);

      G4Box*  Wafer = new G4Box("Wafer",
				DSSSD_Wafer_Length/2.,
				DSSSD_Wafer_Width/2.,
				m_ThicknessDSSSD[i]/2.);

      G4Box*  ActiveWafer = new G4Box("ActiveWafer",
				      DSSSD_ActiveWafer_Length/2.,
				      DSSSD_ActiveWafer_Width/2.,
				      m_ThicknessDSSSD[i]/2.);

      G4double DSSSD_PCB_Slot_Width;
      G4double DSSSD_PCB_Slot_Deepness;
      G4double DSSSD_PCB_Slot_Position;
      G4double DSSSD_DetectorSpacing;



      DSSSD_PCB_Slot_Width = DSSSD_PCB_Slot_Width1;
      DSSSD_PCB_Slot_Deepness = DSSSD_PCB_Slot_Deepness1;
      DSSSD_PCB_Slot_Position = DSSSD_PCB_Slot_Position1;
      DSSSD_DetectorSpacing = DSSSD_DetectorSpacing1 ;


      G4Box*  SlotShape = new G4Box("SlotShape",
				    DSSSD_PCB_Slot_Width/2.,
				    DSSSD_PCB_Width/2.+0.1*mm,
				    DSSSD_PCB_Slot_Deepness);


      G4ThreeVector DSSSD_Wafer_Offset ;
      DSSSD_Wafer_Offset = G4ThreeVector(DSSSD_Wafer_Length_Offset1, DSSSD_Wafer_Width_Offset1,0 );


      G4SubtractionSolid* PCB1 = new G4SubtractionSolid("PCB", PCBFull, SlotShape,new G4RotationMatrix,G4ThreeVector(DSSSD_PCB_Slot_Position, 0,0.5*DSSSD_PCB_Thickness));
      G4SubtractionSolid* PCB = new G4SubtractionSolid("PCB", PCB1, WaferShape,new G4RotationMatrix,DSSSD_Wafer_Offset);

      // Master Volume
      G4LogicalVolume* logicDSSSDDetector =
        new G4LogicalVolume(PCB1,m_MaterialVacuum,"logicDSSSDDetector", 0, 0, 0);
      logicDSSSDDetector->SetVisAttributes(G4VisAttributes::GetInvisible());
      // Sub Volume PCB
      G4LogicalVolume* logicPCB =
        new G4LogicalVolume(PCB,m_MaterialPCB,"logicPCB", 0, 0, 0);
      logicPCB->SetVisAttributes(PCBVisAtt);

      // Sub Volume Wafer
      G4LogicalVolume* logicWafer =
        new G4LogicalVolume(Wafer,m_MaterialSilicon,"logicWafer", 0, 0, 0);
      G4LogicalVolume* logicActiveWafer =
        new G4LogicalVolume(ActiveWafer,m_MaterialSilicon,"logicActiveWafer", 0, 0, 0);

      logicWafer->SetVisAttributes(SiliconVisAtt);
      logicActiveWafer->SetVisAttributes(SiliconVisAtt);

      // Place the sub volume in the master volume
      new G4PVPlacement(new G4RotationMatrix(0,0,0),
			G4ThreeVector(0,0,0),
			logicPCB,"DSSSD_PCB",logicDSSSDDetector,false,DetNbr);

      if(m_ThicknessDSSSD[i]>0){
        new G4PVPlacement(new G4RotationMatrix(0,0,0),
			  DSSSD_Wafer_Offset+G4ThreeVector(0,0,0.5*DSSSD_PCB_Thickness-0.5*m_ThicknessDSSSD[i]),
			  logicWafer,"DSSSD_Wafer",logicDSSSDDetector,false,DetNbr);
        new G4PVPlacement(new G4RotationMatrix(0,0,0),
			  G4ThreeVector(0,0,0),
			  logicActiveWafer,"DSSSD_ActiveWafer",logicWafer,false,DetNbr);

        logicActiveWafer->SetSensitiveDetector(m_DSSSDScorer);
      }

      ///////////////////////////////////////////////////////////////////////////////////
      // Place the detector in the world
      // Position of the center of the PCB

      G4ThreeVector DetectorPosition;

      G4RotationMatrix* DetectorRotation= new G4RotationMatrix;
      // The Rotation Matrix is different for each detector
      double Z= 0;
      //      if (j < m_NumLayers[i]/2)
      Z = m_ZDSSSD[i]+(2.*j+1-(double)m_NumLayersDSSSD[i])/2.*m_GapDistanceDSSSD[i];
      G4cout<<"Z: "<<Z<<G4endl;
      //      DetectorPosition.transform(*DetectorRotation);
      DetectorPosition+=G4ThreeVector(0,0,Z);


      new G4PVPlacement(G4Transform3D(*DetectorRotation,DetectorPosition), logicDSSSDDetector,"DSSSD",world,false,DetNbr);

    }
  }
}
void WAS3ABi::ConstructPassiveStopper(G4LogicalVolume* world){
  for(unsigned int i = 0 ; i < m_ZPassive.size() ; i++){
	
    G4LogicalVolume* logicPassiveStopper;
    G4LogicalVolume* logicBlock;
	if(m_WidthPassive.size()>0){
    G4Box* PassiveStopper = new G4Box("PassiveStopper",
				      m_WidthPassive[i]/2.,
				      m_HeightPassive[i]/2.,
				      m_ThicknessPassive[i]/2.);
	   logicPassiveStopper =
      new G4LogicalVolume(PassiveStopper,m_MaterialVacuum,"logicPassiveStopper", 0, 0, 0);
    logicPassiveStopper->SetVisAttributes(G4VisAttributes::GetInvisible());
    	   logicBlock =
      new G4LogicalVolume(PassiveStopper,
			  MaterialManager::getInstance()->GetMaterialFromLibrary(m_MaterialPassive[i]),
			  "logicBlock", 0, 0, 0);
	}else if(m_RadiPassive.size()>0){
    G4Tubs* PassiveStopper = new G4Tubs("PassvieStopper",
		    			0.0,
					m_RadiPassive[i],
					m_ThicknessPassive[i]/2.,
					0.0*deg,
					360.0*deg);
	   logicPassiveStopper =
      new G4LogicalVolume(PassiveStopper,m_MaterialVacuum,"logicPassiveStopper", 0, 0, 0);
    logicPassiveStopper->SetVisAttributes(G4VisAttributes::GetInvisible());
    	   logicBlock =
      new G4LogicalVolume(PassiveStopper,
			  MaterialManager::getInstance()->GetMaterialFromLibrary(m_MaterialPassive[i]),
			  "logicBlock", 0, 0, 0);
	}
    
    // Sub Volume Block

    logicBlock->SetVisAttributes(PCBVisAtt);

    // Place the sub volume in the master volume
    new G4PVPlacement(new G4RotationMatrix(0,0,0),
		      G4ThreeVector(0,0,0),
		      logicBlock,"DSSSD_Block",logicPassiveStopper,false,0);
      
    ///////////////////////////////////////////////////////////////////////////////////
    // Place the detector in the world
    // Position of the center of the PCB

    G4ThreeVector StopperPosition;
    G4RotationMatrix* StopperRotation= new G4RotationMatrix;
    // The Rotation Matrix is different for each detector
    StopperPosition+=G4ThreeVector(0,0,m_ZPassive[i]);


    new G4PVPlacement(G4Transform3D(*StopperRotation,StopperPosition), logicPassiveStopper,"Passive",world,false,0);

  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void WAS3ABi::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("WAS3ABi")){
    pTree->Branch("WAS3ABi", "TWAS3ABiData", &m_Event) ;
  }
  pTree->SetBranchAddress("WAS3ABi", &m_Event) ;

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void WAS3ABi::ReadSensitive(const G4Event* ){
  m_Event->Clear();

  ///////////
  // DSSSD
  DSSDScorers::PS_Rectangle* DSSSDScorer = (DSSDScorers::PS_Rectangle*) m_DSSSDScorer->GetPrimitive(0);


  // Loop on the DSSSD map
  unsigned int sizeFront= DSSSDScorer->GetLengthMult();

  for (unsigned int i=0 ; i<sizeFront ; i++){

    double Energy = DSSSDScorer->GetEnergyLength(i);

    if(Energy>EnergyThreshold){
      double Time       = DSSSDScorer->GetTimeLength(i);
      int DetNbr        = DSSSDScorer->GetDetectorLength(i);
      int StripFront    = DSSSDScorer->GetStripLength(i);
      //G4cout<<DetNbr<<G4endl;
      m_Event->SetFront(DetNbr,
			StripFront,
			RandGauss::shoot(Energy, ResoEnergy)/keV,
			RandGauss::shoot(Time, ResoTime),
			RandGauss::shoot(Time, ResoTime));
    }
  } 

  unsigned int sizeBack= DSSSDScorer->GetWidthMult();
  for (unsigned int i=0 ; i<sizeBack ; i++){

    double Energy = DSSSDScorer->GetEnergyWidth(i);

    if(Energy>EnergyThreshold){
      double Time       = DSSSDScorer->GetTimeWidth(i);
      int DetNbr        = DSSSDScorer->GetDetectorWidth(i);
      int StripBack    = DSSSDScorer->GetStripWidth(i);

      m_Event->SetBack(DetNbr,
		       DSSSD_Wafer_Back_NumberOfStrip-StripBack+1,
		       RandGauss::shoot(Energy, ResoEnergy)/keV,
		       RandGauss::shoot(Time, ResoTime),
		       RandGauss::shoot(Time, ResoTime));
    }
  }
  // clear map for next event
  DSSSDScorer->clear();


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void WAS3ABi::InitializeScorers(){
  //   Silicon Associate Scorer
  bool already_exist = false;
  m_DSSSDScorer = CheckScorer("WAS3ABi_DSSSDScorer",already_exist);
  // if the scorer were created previously nothing else need to be made
  if(already_exist) return;

  G4VPrimitiveScorer* DSSSDScorer =
    new  DSSDScorers::PS_Rectangle("WAS3ABiDSSSD",0,
				   DSSSD_ActiveWafer_Length,
				   DSSSD_ActiveWafer_Width,
				   DSSSD_Wafer_Front_NumberOfStrip ,
				   DSSSD_Wafer_Back_NumberOfStrip);

  
  G4VPrimitiveScorer* InterScorerDSSSD = 
    new  InteractionScorers::PS_Interactions("WAS3ABiDSSSDInteractionScorer",ms_InterCoord,0);
 
  
  //and register it to the multifunctionnal detector
  m_DSSSDScorer->RegisterPrimitive(DSSSDScorer);
  m_DSSSDScorer->RegisterPrimitive(InterScorerDSSSD);

  G4SDManager::GetSDMpointer()->ListTree();
  //   Add All Scorer to the Global Scorer Manager
  G4SDManager::GetSDMpointer()->AddNewDetector(m_DSSSDScorer) ;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////
/////////////////Material Definition ///////////////////////////
////////////////////////////////////////////////////////////////
void WAS3ABi::InitializeMaterial(){
  m_MaterialSilicon = MaterialManager::getInstance()->GetMaterialFromLibrary("Si");
  m_MaterialPCB = MaterialManager::getInstance()->GetMaterialFromLibrary("PCB");
  m_MaterialVacuum = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");

}


////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* WAS3ABi::Construct(){
  return  (NPS::VDetector*) new WAS3ABi();
}

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_sharc{
  public:
    proxy_nps_sharc(){
      NPS::DetectorFactory::getInstance()->AddToken("WAS3ABi","WAS3ABi");
      NPS::DetectorFactory::getInstance()->AddDetector("WAS3ABi",WAS3ABi::Construct);
    }
  };

  proxy_nps_sharc p_nps_sharc;
}
