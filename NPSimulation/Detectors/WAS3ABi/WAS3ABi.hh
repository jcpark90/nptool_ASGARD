#ifndef WAS3ABi_h
#define WAS3ABi_h 1
/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: a.matta@surrey.ac.uk      *
 *                                                                           *
 * Creation Date  : November 2012                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe the WAS3ABi Silicon detector                           *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/
// C++ header
#include <string>
#include <vector>

// G4 header defining G4 types
#include "globals.hh"

// G4 headers
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4MultiFunctionalDetector.hh"

// NPSimulation header
#include "NPSVDetector.hh"

// NPLib
#include "TWAS3ABiData.h"
#include "NPInputParser.h"
using namespace std;
using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace WAS3ABI{
  // Energy and time Resolution
  const G4double ResoTime    = 0 ;
  const G4double ResoEnergy  = 35*keV ;// = zzkeV of Resolution   //   Unit is MeV/2.35
  const G4double EnergyThreshold = 50*keV;
  // Geometry
  
  // DSSSD //
  // DSSSD PCB
  const G4double DSSSD_PCB_Width  = 51.10*mm;
  const G4double DSSSD_PCB_Length = 84.00*mm;
  const G4double DSSSD_PCB_Thickness = 1.5*mm;
  const G4double DSSSD_PCB_Border_LongSide = 1*mm;
  const G4double DSSSD_PCB_Border_ShortSide = 2*mm;
  
  // Single stage box case (DSSD only)
  const G4double DSSSD_PCB_Slot_Width1 = DSSSD_PCB_Thickness;
  const G4double DSSSD_PCB_Slot_Border1 = 4*mm;
  const G4double DSSSD_PCB_Slot_Deepness1 = DSSSD_PCB_Border_ShortSide;
  
  // DSSSD Wafer
  double DSSSD_ActiveWafer_Width  = 40;
  double DSSSD_ActiveWafer_Length = 60;
  double DSSSD_Wafer_Width  = 44.0;
  double DSSSD_Wafer_Length = 66.0;  
  
  const G4double DSSSD_Wafer_DeadLayer_Thickness = 0.1*um;
  const G4int    DSSSD_Wafer_Front_NumberOfStrip = 60 ;
  const G4int    DSSSD_Wafer_Back_NumberOfStrip = 40 ;
  
  // Compute
  const G4double DSSSD_LeftOver1 =  DSSSD_PCB_Length - DSSSD_PCB_Border_ShortSide - DSSSD_Wafer_Length - DSSSD_PCB_Slot_Border1 - DSSSD_PCB_Slot_Width1 ;
  const G4double DSSSD_Exposed_Length1 = DSSSD_Wafer_Length + DSSSD_PCB_Slot_Border1 ;
  
  const G4double DSSSD_CenterOffset1 = - 0.5 * DSSSD_PCB_Length+DSSSD_PCB_Border_ShortSide+0.5*DSSSD_Exposed_Length1;
  const G4double DSSSD_DetectorSpacing1 = 0.5*DSSSD_Exposed_Length1+0.5*DSSSD_PCB_Slot_Width1;
  const G4double DSSSD_DetectorGap = 5*mm;

  const G4double DSSSD_Wafer_Width_Offset1 = -0.5*DSSSD_PCB_Width + DSSSD_PCB_Border_LongSide + 0.5*DSSSD_Wafer_Width;
  const G4double DSSSD_Wafer_Length_Offset1 = -0.5*DSSSD_PCB_Length + DSSSD_PCB_Border_ShortSide + 0.5*DSSSD_Wafer_Length;
  
  const G4double DSSSD_PCB_Slot_Position1 = 0.5*DSSSD_PCB_Length-DSSSD_LeftOver1 - 0.5*DSSSD_PCB_Slot_Width1;
  
  
  
}

using namespace WAS3ABI ;
class WAS3ABi : public NPS::VDetector
{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
public:
  WAS3ABi() ;
   ~WAS3ABi() ;
  
  ////////////////////////////////////////////////////
  //////// Specific Function of this Class ///////////
  ////////////////////////////////////////////////////
public:
  // To add a WAS3ABi detector
  void AddDSSSDDetector(G4double Z, unsigned int Nlayers, G4double GapDistance, G4double Thickness);
  // To add a Passive stopper
  void AddPassiveStopper(G4double Z, G4String Material, G4double Thickness, G4double Width, G4double Height);
  void AddPassiveCircleStopper(G4double Z, G4String Material, G4double Thickness, G4double R);
  
  // Construct Volume
  void ConstructDSSSDDetector(G4LogicalVolume* world);
  void ConstructPassiveStopper(G4LogicalVolume* world);
  
  ////////////////////////////////////////////////////
  /////////  Inherite from NPS::VDetector class ///////////
  ////////////////////////////////////////////////////
public:
  // Read stream at Configfile to pick-up parameters of detector (Position,...)
  // Called in DetecorConstruction::ReadDetextorConfiguration Method
  void ReadConfiguration(NPL::InputParser) ;
  
  // Construct detector and inialise sensitive part.
  // Called After DetecorConstruction::AddDetector Method
  void ConstructDetector(G4LogicalVolume* world) ;
  
  // Add Detector branch to the EventTree.
  // Called After DetecorConstruction::AddDetector Method
  void InitializeRootOutput() ;
  
  // Read sensitive part and fill the Root tree.
  // Called at in the EventAction::EndOfEventAvtion
  void ReadSensitive(const G4Event* event) ;
  
  ////////////////////////////////////////////////////
  ///////////Event class to store Data////////////////
  ////////////////////////////////////////////////////
private:
  TWAS3ABiData*    m_Event ;
  
  ////////////////////////////////////////////////////
  ///////////////// Scorer Related ///////////////////
  ////////////////////////////////////////////////////
  
private:
  //   Initialize all Scorer
  void InitializeScorers() ;
  
  //   Scorer Associate to the Silicon
  G4MultiFunctionalDetector*   m_DSSSDScorer ;
  
private:
  //    Initialize material used in detector definition
  void InitializeMaterial();
  
  //   List of material
  G4Material* m_MaterialSilicon ;
  G4Material* m_MaterialVacuum  ;
  G4Material* m_MaterialPCB     ;
  
  ////////////////////////////////////////////////////
  ///////////////Private intern Data//////////////////
  ////////////////////////////////////////////////////
private:
  // True if the detector is a DSSSD, false if a passive stopper
  vector<bool>   m_Type  ;

  // Used for Quadrant detectors
  vector<G4ThreeVector>   m_Pos   ; // R , Phi , Z
  
  // Used for DSSSD detectors
  vector<G4double>   m_ZDSSSD;
  vector<unsigned int>   m_NumLayersDSSSD;
  vector<G4double>   m_GapDistanceDSSSD;  
  vector<G4double>   m_ThicknessDSSSD;

  vector<G4double>   m_ZPassive;
  vector<G4String>   m_MaterialPassive;
  vector<G4double>   m_ThicknessPassive;
  vector<G4double>   m_WidthPassive;
  vector<G4double>   m_HeightPassive;
  vector<G4double>   m_RadiPassive;
  
private:/// Visualisation Attribute:
  // Dark Grey
   G4VisAttributes* SiliconVisAtt  ;
  // Green
   G4VisAttributes* PCBVisAtt;
  // Light Grey
   G4VisAttributes* FrameVisAtt ;
  
public:
    static NPS::VDetector* Construct();
};
#endif
