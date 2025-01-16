/*****************************************************************************
 * Copyright (C) 2009-2024   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Jason Park  contact address: jcpark@ibs.re.kr                        *
 *                                                                           *
 * Creation Date  : 9월 2024                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  ASGARD simulation                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

// C++ header

#include <sstream>
#include <cmath>
#include <limits>
//G4 Geometry object
#include "G4Tubs.hh"
#include "G4Box.hh"

//G4 sensitive
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"

//G4 various object
#include "G4Material.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4RunManager.hh"
#include "G4ios.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4ThreeVector.hh"

// NPTool header
#include "ASGARD.hh"
#include "ASGARDScorers.hh"
#include "CalorimeterScorers.hh"
#include "InteractionScorers.hh"
#include "RootOutput.h"
#include "MaterialManager.hh"
#include "NPSDetectorFactory.hh"
#include "NPOptionManager.h"
#include "NPSHitsMap.hh"
// CLHEP header
#include "CLHEP/Random/RandGauss.h"

using namespace std;
using namespace CLHEP;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace {
  // Ge crystal
  // Cylindrical part
  const G4double CrystalOuterRadius   = 30.0*mm; // outer radius for crystal
  const G4double CrystalInnerRadius   =  5.0*mm; // inner radius for hole in crystal
  const G4double CrystalLength        = 90.0*mm; // crystal length
  const G4double CrystalHoleDepth     = 15.0*mm; // depth at which starts the hole

  const G4double CrystalReducedRadius = 27.2*mm;
  const G4double CrystalEdgeRadius = 29.3*mm;

  
  //const G4double CrystaHoleRadius 		= 0*cm;
  const G4double CrystalInterDistance =  0.6*mm; // Distance between two crystal, confirmed by drawing

  // Squared part
  const G4double CrystalWidth         = CrystalReducedRadius + CrystalEdgeRadius;  	// Width of one crystal
  const G4double widthface = 44.0*mm; // 45.5 mm default, JP: reduce to avoid overlaps at 11 cm

  // tapered part
  const G4double SegmentDepth = 31.*mm; // depth, see Ref. [16] in  NIM A 543 (2005), 431
  const G4double TaperDepth = 30.*mm;
  //  const G4double BezelOuterAngle = atan(TaperMagnitude/TaperDepth)*180./M_PI; // 21.8 degrees from NIM A 540 (2005), 348; but 22.5 degrees in NIM A 543 (2005), 431
  const G4double BezelOuterAngle = 22.5; // from Tigress, half of 45 degrees

  // segment part
  const G4double SegmentFullWidth = 27.*mm;
  const G4double SegmentFrontWidth = 20.9*mm;
  const G4double TaperMagnitude = 2.*(SegmentFullWidth-SegmentFrontWidth);  
  const G4double BezelInnerAngle = atan(TaperMagnitude*0.5/TaperDepth)*180./M_PI; // 11.3 degrees

  
  // Exogam Stuff
  const G4double CrystalEdgeOffset1  = CrystalReducedRadius; // distance of the edge from the center of the crystal
  const G4double CrystalEdgeOffset2  = CrystalEdgeRadius; // distance of the edge from the center of the crystal

  const G4double CapsuleWidth        = 1.0*mm;   // capsule width
  const G4double CapsuleLength       = 110.*mm;   // capsule length
  const G4double CapsuleEdgeDepth    = 3.3*cm;   // same as crystal !
  const G4double CrystalToCapsule    = 4.*mm;   // to be adjusted ..

  //const G4double BGOLength           = 120.0*mm;
  //const G4double BGOWidth            = 25.0*mm;

  //const G4double CsILength           = 20.0*mm;


}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// ASGARD Specific Method
ASGARD::ASGARD(){
  InitializeMaterial();

  
  m_ASGARDData = new TASGARDData() ;
  m_ASGARDScorer = 0;

  BlueVisAtt   = new G4VisAttributes(G4Colour(0, 0, 1)) ;
  GreenVisAtt  = new G4VisAttributes(G4Colour(0, 1, 0)) ;
  RedVisAtt    = new G4VisAttributes(G4Colour(1, 0, 0)) ;
  TestVisAtt1   = new G4VisAttributes(G4Colour(0, 1, 1)) ;
  TestVisAtt2   = new G4VisAttributes(G4Colour(1, 0, 1)) ;
  TestVisAtt3   = new G4VisAttributes(G4Colour(1, 1, 0)) ;
  WhiteVisAtt  = new G4VisAttributes(G4Colour(1, 1, 1)) ;
  TrGreyVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 0.5)) ;

  m_LogicClover = 0;

}

ASGARD::~ASGARD(){
}

// Add clover at a free position in space with coordinate
// in spherical coordinate
// Beta are the three angles of rotation in the Clover frame

void ASGARD::AddClover(int CloverId,double R,double Theta,double Phi,double BetaX,double BetaY,double BetaZ){

  m_CloverId.push_back(CloverId);
  m_R.push_back(R);
  m_Theta.push_back(Theta);
  m_Phi.push_back(Phi);
  m_BetaX.push_back(BetaX);
  m_BetaY.push_back(BetaY);
  m_BetaZ.push_back(BetaZ);

}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void ASGARD::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithTokenAndValue("ASGARD","Clover");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " free clovers found " << endl; 

  vector<string> token = {"CloverID","R","Theta","Phi","Beta"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(token)){
      double R = blocks[i]->GetDouble("R","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");
      vector<double> beta = blocks[i]->GetVectorDouble("Beta","deg");
      int     id = blocks[i]->GetInt("CloverID");
      AddClover(id,R,Theta,Phi,beta[0],beta[1],beta[2]);
    }

    else{
      cout << "Warning: check your input file formatting " << endl;
    }
  }

  
}
// Return a G4VSolid modeling the Crystal
// JP: need to divide crystal into 8 segments, need to reduce cutting steps
// Counting from 1:
// 1-4: top left and anticlockwise, front volume (closest to target)
// 5-8: top left and anticlockwise, back volume (furthest from target)
G4LogicalVolume* ASGARD::ConstructSegment(G4int seg){
  // start with cylinder and cutting box
  G4Tubs* Crystal_Cylinder = new G4Tubs("Crystal_Cylinder", 0, CrystalOuterRadius, CrystalLength*0.5, 0, 2*M_PI);
  G4Box* Crystal_Box1 = new G4Box("Crystal_Box1", CrystalWidth*0.6, CrystalWidth*0.6,CrystalLength*0.6); // rectangular prism which can house the cylinder

  G4RotationMatrix* BoxRotation = new G4RotationMatrix(0,0,0);
  G4RotationMatrix* BoxRotationOuterLeftBezel = new G4RotationMatrix(0,0,0);
  BoxRotationOuterLeftBezel->rotate(BezelOuterAngle*deg,G4ThreeVector(1,0,0));
  G4RotationMatrix* BoxRotationOuterTopBezel = new G4RotationMatrix(0,0,0);
  BoxRotationOuterTopBezel->rotate(BezelOuterAngle*deg,G4ThreeVector(0,1,0));

  G4RotationMatrix* BoxRotationInnerRightBezel = new G4RotationMatrix(0,0,0);
  BoxRotationInnerRightBezel->rotate(BezelInnerAngle*deg,G4ThreeVector(1,0,0));
  G4RotationMatrix* BoxRotationInnerBottomBezel = new G4RotationMatrix(0,0,0);
  BoxRotationInnerBottomBezel->rotate(BezelInnerAngle*deg,G4ThreeVector(0,1,0));

  // Central Hole for cold finger (middle of crystal from back)
  G4Tubs* Crystal_Hole = new G4Tubs("Crystal_Hole", 0, CrystalInnerRadius, (CrystalLength-CrystalHoleDepth)*0.5, 0, 2*M_PI);
  G4SubtractionSolid* Crystal_withHole = new G4SubtractionSolid("Crystal_withHole",Crystal_Cylinder,Crystal_Hole,BoxRotation,G4ThreeVector(0,0,CrystalHoleDepth*0.5));

  // Flat surface on the side cut out on 4 sides and z
  // apply variations in these cuts for each of the 8 segments
  // need one more cut for depth (z)
  // orientation descriptors as seen from target
  G4SubtractionSolid* Crystal_RightCut;
  G4SubtractionSolid* Crystal_LeftCut;
  G4SubtractionSolid* Crystal_TopCut;
  G4SubtractionSolid* Crystal_BottomCut;
  G4SubtractionSolid* Crystal_TopBezelCut;
  G4SubtractionSolid* Crystal_LeftBezelCut;
  
  
  
  G4SubtractionSolid* Segment_final;
  
  G4LogicalVolume* logicSegment;
  // CrystalWidth*0.6 = 56.5(0.6) = 33.9
  // CrystalLength*0.6 = 90(0.6) = 54.0
  //45 - 54 = -9 mm, add depth + 9 mm which is 0.1*CrystalLength; i.e. 31 mm depth front volume = 31+9 = 40 mm shift back
  // for back segments, 79 mm shift front
  // Full box dimensions = 67.8 mm
  
  // cut coordinate system: (+down/-up, +left/-right, +back/-front) <-volume to be cut
  // start dis1./costing here into 8 segments (4 front + 4 back)
  if (seg==0){ // top left, front
    // apply tilted cuts for outer corner segment
    Crystal_BottomCut =
      new G4SubtractionSolid("Crystal_BottomCut",Crystal_withHole,Crystal_Box1,
			     BoxRotationInnerBottomBezel,
			     G4ThreeVector(2.*CrystalOuterRadius-SegmentFullWidth-CrystalEdgeRadius+0.6*CrystalWidth*(1./cos(BezelInnerAngle*M_PI/180.)-1.)+CrystalWidth*0.6,0,-0.5*CrystalLength)); // in mm
    Crystal_RightCut =
      new G4SubtractionSolid("Crystal_RightCut",Crystal_BottomCut,Crystal_Box1,
			     BoxRotationInnerRightBezel,
			     G4ThreeVector(0,-(2.*CrystalOuterRadius-SegmentFullWidth-CrystalEdgeRadius+0.6*CrystalWidth*(1./cos(BezelInnerAngle*M_PI/180.)-1.))-CrystalWidth*0.6,-0.5*CrystalLength));
    
    // left bezel (applied on segments 0, 1)
    Crystal_LeftBezelCut =
      new G4SubtractionSolid("Crystal_LeftBezelCut",Crystal_RightCut,Crystal_Box1,
			     BoxRotationOuterLeftBezel,
			     G4ThreeVector(0,(CrystalEdgeRadius-TaperMagnitude+0.6*CrystalWidth*(1./cos(BezelOuterAngle*M_PI/180.)-1.))+CrystalWidth*0.6,-CrystalLength*0.5));

    // need to slice off non-tapered edge at left
    Crystal_LeftCut= new G4SubtractionSolid("Crystal_LeftCut",Crystal_LeftBezelCut,Crystal_Box1,BoxRotation,G4ThreeVector(0, CrystalEdgeRadius+CrystalWidth*0.6,0));

    // top bezel (applied on segments 0, 3)
    Crystal_TopBezelCut =
      new G4SubtractionSolid("Crystal_TopBezelCut",Crystal_LeftCut,Crystal_Box1,
			     BoxRotationOuterTopBezel,
			     G4ThreeVector(-(CrystalEdgeRadius-TaperMagnitude+0.6*CrystalWidth*(1./cos(BezelOuterAngle*M_PI/180.)-1.))-CrystalWidth*0.6,0,-CrystalLength*0.5));

    // need to slice off non-tapered edge at top
    Crystal_TopCut= new G4SubtractionSolid("Crystal_TopCut",Crystal_TopBezelCut,Crystal_Box1,BoxRotation,G4ThreeVector(-CrystalEdgeRadius-CrystalWidth*0.6,0,0));

    
    Segment_final =
      new G4SubtractionSolid("Segment_final",Crystal_TopCut,Crystal_Box1,
			     BoxRotation,
			     G4ThreeVector(0,0,SegmentDepth+0.1*CrystalLength)); // cut back volume
    
  }
  else if (seg==1){ // bottom left, front
    //    Crystal_TopCut = new G4SubtractionSolid("Crystal_TopCut",Crystal_withHole,Crystal_Box1,BoxRotation,G4ThreeVector(-CrystalWidth*0.6,0,0)); // in mm
    //    Crystal_RightCut = new G4SubtractionSolid("Crystal_RightCut",Crystal_TopCut,Crystal_Box1,BoxRotation,G4ThreeVector(0,-CrystalWidth*0.6,0));
    Crystal_TopCut = new G4SubtractionSolid("Crystal_TopCut",Crystal_withHole,Crystal_Box1,BoxRotationInnerBottomBezel,
					    G4ThreeVector((CrystalOuterRadius-SegmentFullWidth+0.6*CrystalWidth*(1./cos(BezelInnerAngle*M_PI/180.)-1.)-CrystalWidth*0.6),0,-0.5*CrystalLength)); // in mm
  
    Crystal_RightCut = new G4SubtractionSolid("Crystal_RightCut",Crystal_TopCut,Crystal_Box1,BoxRotationInnerRightBezel,
					      G4ThreeVector(0,-(2.*CrystalOuterRadius-SegmentFullWidth-CrystalEdgeRadius+0.6*CrystalWidth*(1./cos(BezelInnerAngle*M_PI/180.)-1.))-CrystalWidth*0.6,-0.5*CrystalLength));

    Crystal_BottomCut = new G4SubtractionSolid("Crystal_BottomCut",Crystal_RightCut,Crystal_Box1,BoxRotation,G4ThreeVector(CrystalReducedRadius+CrystalWidth*0.6,0,0));

    // left bezel (applied on segments 0, 1, 4, 5)
    Crystal_LeftBezelCut= new G4SubtractionSolid("Crystal_LeftBezelCut",Crystal_BottomCut,Crystal_Box1,BoxRotationOuterLeftBezel,
						 G4ThreeVector(0, (CrystalEdgeRadius-TaperMagnitude+0.6*CrystalWidth*(1./cos(BezelOuterAngle*M_PI/180.)-1.))+CrystalWidth*0.6,-CrystalLength*0.5));
    // need to slice off non-tapered edge at left
    Crystal_LeftCut= new G4SubtractionSolid("Crystal_LeftCut",Crystal_LeftBezelCut,Crystal_Box1,BoxRotation,G4ThreeVector(0, CrystalEdgeRadius+CrystalWidth*0.6,0));
    
    Segment_final= new G4SubtractionSolid("Segment_final",Crystal_LeftCut,Crystal_Box1,BoxRotation,G4ThreeVector(0,0,SegmentDepth+0.1*CrystalLength)); // cut back volume
  }
  else if (seg==2){ // bottom right, front
    Crystal_TopCut = new G4SubtractionSolid("Crystal_TopCut",Crystal_withHole,Crystal_Box1,BoxRotationInnerBottomBezel,
					    G4ThreeVector((CrystalOuterRadius-SegmentFullWidth+0.6*CrystalWidth*(1./cos(BezelInnerAngle*M_PI/180.)-1.)-CrystalWidth*0.6),0,-0.5*CrystalLength));
    Crystal_LeftCut = new G4SubtractionSolid("Crystal_LeftCut",Crystal_TopCut,Crystal_Box1,BoxRotationInnerRightBezel,
					      G4ThreeVector(0,-(CrystalOuterRadius-SegmentFullWidth+0.6*CrystalWidth*(1./cos(BezelInnerAngle*M_PI/180.)-1.))+CrystalWidth*0.6,-0.5*CrystalLength));

    
    Crystal_RightCut = new G4SubtractionSolid("Crystal_RightCut",Crystal_LeftCut,Crystal_Box1,BoxRotation,G4ThreeVector(0,-CrystalReducedRadius-CrystalWidth*0.6,0));
     Crystal_BottomCut = new G4SubtractionSolid("Crystal_BottomCut",Crystal_RightCut,Crystal_Box1,BoxRotation,G4ThreeVector(CrystalReducedRadius+CrystalWidth*0.6,0,0));
    
    Segment_final= new G4SubtractionSolid("Segment_final",Crystal_BottomCut,Crystal_Box1,BoxRotation,G4ThreeVector(0,0,SegmentDepth+0.1*CrystalLength)); // cut back volume
    
  }
  else if (seg==3){ // top right, front
    Crystal_LeftCut = new G4SubtractionSolid("Crystal_LeftCut",Crystal_withHole,Crystal_Box1,BoxRotationInnerRightBezel,
					     G4ThreeVector(0,-(CrystalOuterRadius-SegmentFullWidth+0.6*CrystalWidth*(1./cos(BezelInnerAngle*M_PI/180.)-1.))+CrystalWidth*0.6,-0.5*CrystalLength));
    Crystal_BottomCut =
      new G4SubtractionSolid("Crystal_BottomCut",Crystal_LeftCut,Crystal_Box1,
			     BoxRotationInnerBottomBezel,
			     G4ThreeVector(2.*CrystalOuterRadius-SegmentFullWidth-CrystalEdgeRadius+0.6*CrystalWidth*(1./cos(BezelInnerAngle*M_PI/180.)-1.)+CrystalWidth*0.6,0,-0.5*CrystalLength)); // in mm
    Crystal_RightCut = new G4SubtractionSolid("Crystal_RightCut",Crystal_BottomCut,Crystal_Box1,BoxRotation,G4ThreeVector(0,-CrystalReducedRadius-CrystalWidth*0.6,0));

    Crystal_TopBezelCut =
      new G4SubtractionSolid("Crystal_TopBezelCut",Crystal_RightCut,Crystal_Box1,
			     BoxRotationOuterTopBezel,
			     G4ThreeVector(-(CrystalEdgeRadius-TaperMagnitude+0.6*CrystalWidth*(1./cos(BezelOuterAngle*M_PI/180.)-1.))-CrystalWidth*0.6,0,-CrystalLength*0.5));

    // need to slice off non-tapered edge at top
    Crystal_TopCut= new G4SubtractionSolid("Crystal_TopCut",Crystal_TopBezelCut,Crystal_Box1,BoxRotation,G4ThreeVector(-CrystalEdgeRadius-CrystalWidth*0.6,0,0));

    
    Segment_final= new G4SubtractionSolid("Segment_final",Crystal_TopCut,Crystal_Box1,BoxRotation,G4ThreeVector(0,0,SegmentDepth+0.1*CrystalLength)); // cut back volume
    
  }
  else if (seg==4){ // top left, back
    Crystal_BottomCut = new G4SubtractionSolid("Crystal_BottomCut",Crystal_withHole,Crystal_Box1,BoxRotation,G4ThreeVector(CrystalWidth*0.6-(SegmentFullWidth-CrystalReducedRadius),0,0)); // in mm
    Crystal_RightCut = new G4SubtractionSolid("Crystal_RightCut",Crystal_BottomCut,Crystal_Box1,BoxRotation,G4ThreeVector(0,-CrystalWidth*0.6+(SegmentFullWidth-CrystalReducedRadius),0));
    Crystal_LeftCut= new G4SubtractionSolid("Crystal_LeftCut",Crystal_RightCut,Crystal_Box1,BoxRotation,G4ThreeVector(0, CrystalEdgeRadius+CrystalWidth*0.6,0));
    Crystal_TopCut= new G4SubtractionSolid("Crystal_TopCut",Crystal_LeftCut,Crystal_Box1,BoxRotation,G4ThreeVector(-CrystalEdgeRadius-CrystalWidth*0.6,0,0));
    
    Segment_final= new G4SubtractionSolid("Segment_final",Crystal_TopCut,Crystal_Box1,BoxRotation,G4ThreeVector(0,0,-(CrystalLength-SegmentDepth)-0.1*CrystalLength)); // cut front volume

  }
  else if (seg==5){ // bottom left, back
    Crystal_BottomCut = new G4SubtractionSolid("Crystal_BottomCut",Crystal_withHole,Crystal_Box1,BoxRotation,G4ThreeVector(CrystalWidth*0.6+CrystalReducedRadius,0,0)); // in mm
    Crystal_RightCut = new G4SubtractionSolid("Crystal_RightCut",Crystal_BottomCut,Crystal_Box1,BoxRotation,G4ThreeVector(0,-CrystalWidth*0.6+(SegmentFullWidth-CrystalReducedRadius),0));
    Crystal_LeftCut= new G4SubtractionSolid("Crystal_LeftCut",Crystal_RightCut,Crystal_Box1,BoxRotation,G4ThreeVector(0, CrystalEdgeRadius+CrystalWidth*0.6,0));
    Crystal_TopCut= new G4SubtractionSolid("Crystal_TopCut",Crystal_LeftCut,Crystal_Box1,BoxRotation,G4ThreeVector(-(SegmentFullWidth-CrystalReducedRadius)-CrystalWidth*0.6,0,0));
    
    Segment_final= new G4SubtractionSolid("Segment_final",Crystal_TopCut,Crystal_Box1,BoxRotation,G4ThreeVector(0,0,-(CrystalLength-SegmentDepth)-0.1*CrystalLength)); // cut front volume
  }
  else if (seg==6){ // bottom right, back
    Crystal_BottomCut = new G4SubtractionSolid("Crystal_BottomCut",Crystal_withHole,Crystal_Box1,BoxRotation,G4ThreeVector(CrystalWidth*0.6+CrystalReducedRadius,0,0)); // in mm
    Crystal_RightCut = new G4SubtractionSolid("Crystal_RightCut",Crystal_BottomCut,Crystal_Box1,BoxRotation,G4ThreeVector(0,-CrystalWidth*0.6-CrystalReducedRadius,0));
    Crystal_LeftCut= new G4SubtractionSolid("Crystal_LeftCut",Crystal_RightCut,Crystal_Box1,BoxRotation,G4ThreeVector(0, (SegmentFullWidth-CrystalReducedRadius)+CrystalWidth*0.6,0));
    Crystal_TopCut= new G4SubtractionSolid("Crystal_TopCut",Crystal_LeftCut,Crystal_Box1,BoxRotation,G4ThreeVector(-(SegmentFullWidth-CrystalReducedRadius)-CrystalWidth*0.6,0,0));
    
    Segment_final= new G4SubtractionSolid("Segment_final",Crystal_TopCut,Crystal_Box1,BoxRotation,G4ThreeVector(0,0,-(CrystalLength-SegmentDepth)-0.1*CrystalLength)); // cut front volume
  }
  else if (seg==7){ // top right, back
    Crystal_BottomCut = new G4SubtractionSolid("Crystal_BottomCut",Crystal_withHole,Crystal_Box1,BoxRotation,G4ThreeVector(CrystalWidth*0.6-(SegmentFullWidth-CrystalReducedRadius),0,0)); // in mm
    Crystal_RightCut = new G4SubtractionSolid("Crystal_RightCut",Crystal_BottomCut,Crystal_Box1,BoxRotation,G4ThreeVector(0,-CrystalWidth*0.6-CrystalReducedRadius,0));
    Crystal_LeftCut= new G4SubtractionSolid("Crystal_LeftCut",Crystal_RightCut,Crystal_Box1,BoxRotation,G4ThreeVector(0, (SegmentFullWidth-CrystalReducedRadius)+CrystalWidth*0.6,0));
    Crystal_TopCut= new G4SubtractionSolid("Crystal_TopCut",Crystal_LeftCut,Crystal_Box1,BoxRotation,G4ThreeVector(-CrystalEdgeRadius-CrystalWidth*0.6,0,0));
    
    Segment_final= new G4SubtractionSolid("Segment_final",Crystal_TopCut,Crystal_Box1,BoxRotation,G4ThreeVector(0,0,-(CrystalLength-SegmentDepth)-0.1*CrystalLength)); // cut front volume
  }
    
  
  else{ 
    // give logical volume to final shape
    logicSegment =
      new G4LogicalVolume(Crystal_LeftBezelCut,m_MaterialGe,"LogicCrystal", 0, 0, 0);
    
  }
  logicSegment =
    new G4LogicalVolume(Segment_final,m_MaterialGe,Form("LogicSegment%d", seg), 0, 0, 0);
  logicSegment->SetSensitiveDetector(m_ASGARDScorer);
  
  return  logicSegment;

}

G4LogicalVolume* ASGARD::ConstructCrystal(){
  G4Tubs* Crystal_Cylinder = new G4Tubs("Crystal_Cylinder", 0, CrystalOuterRadius, CrystalLength*0.5, 0, 2*M_PI);
  G4Box* Crystal_Box1 = new G4Box("Crystal_Box1", CrystalWidth*0.6, CrystalWidth*0.6,CrystalLength*0.6); // rectangular prism which can house the cylinder

  G4RotationMatrix* BoxRotation = new G4RotationMatrix(0,0,0);

  // Central Hole for cold finger (middle of crystal from back), only segments 4-7
  G4Tubs* Crystal_Hole = new G4Tubs("Crystal_Hole", 0, CrystalInnerRadius, (CrystalLength-CrystalHoleDepth)*0.5, 0, 2*M_PI);
  G4SubtractionSolid* Crystal_withHole = new G4SubtractionSolid("Crystal_withHole",Crystal_Cylinder,Crystal_Hole,BoxRotation,G4ThreeVector(-(SegmentFullWidth-CrystalReducedRadius),(SegmentFullWidth-CrystalReducedRadius),CrystalHoleDepth*0.5));

  // Flat surface on the side cut out on 4 sides and z
  // apply variations in these cuts for each of the 8 segments
  // need one more cut for depth (z)
  // orientation descriptors as seen from target
  G4SubtractionSolid* Crystal_BottomCut = new G4SubtractionSolid("Crystal_BottomCut",Crystal_withHole,Crystal_Box1,BoxRotation,G4ThreeVector(CrystalReducedRadius+CrystalWidth*0.6,0,0)); // in mm
  G4SubtractionSolid* Crystal_LeftCut = new G4SubtractionSolid("Crystal_LeftCut",Crystal_BottomCut,Crystal_Box1,BoxRotation,G4ThreeVector(-CrystalEdgeRadius-CrystalWidth*0.6,0,0));
  G4SubtractionSolid* Crystal_TopCut = new G4SubtractionSolid("Crystal_TopCut",Crystal_LeftCut,Crystal_Box1,BoxRotation,G4ThreeVector(0,CrystalEdgeRadius+CrystalWidth*0.6,0));
  G4SubtractionSolid* Crystal_RightCut = new G4SubtractionSolid("Crystal_RightCut",Crystal_TopCut,Crystal_Box1,BoxRotation,G4ThreeVector(0,-CrystalReducedRadius-CrystalWidth*0.6,0));

  // Bezel (tapered cut) on the outer 2 surfaces, affecting 6 segments
  G4RotationMatrix* BoxRotationTopBezel = new G4RotationMatrix(0,0,0);
  BoxRotationTopBezel->rotate(BezelOuterAngle*deg,G4ThreeVector(1,0,0));
  // top bezel (applied on segments 0, 3)
  G4SubtractionSolid* Crystal_TopBezelCut= new G4SubtractionSolid("Crystal_TopBezelCut",Crystal_RightCut,Crystal_Box1,BoxRotationTopBezel,G4ThreeVector(0,(CrystalEdgeRadius-TaperMagnitude+0.6*CrystalWidth*(1./cos(BezelOuterAngle*M_PI/180.)-1.))+CrystalWidth*0.6,-CrystalLength*0.5));
  G4RotationMatrix* BoxRotationLeftBezel = new G4RotationMatrix(0,0,0);
  BoxRotationLeftBezel->rotate(BezelOuterAngle*deg,G4ThreeVector(0,1,0));
  // left bezel (applied on segments 0, 1)
  G4SubtractionSolid* Crystal_LeftBezelCut= new G4SubtractionSolid("Crystal_LeftBezelCut",Crystal_TopBezelCut,Crystal_Box1,BoxRotationLeftBezel,G4ThreeVector(-(CrystalEdgeRadius-TaperMagnitude+0.6*CrystalWidth*(1./cos(BezelOuterAngle*M_PI/180.)-1.))-CrystalWidth*0.6,0,-CrystalLength*0.5));

  G4LogicalVolume* logicCrystal =
    new G4LogicalVolume(Crystal_LeftBezelCut,m_MaterialGe,"LogicCrystal", 0, 0, 0);

  // construct segments here? same? 
  
  return  logicCrystal;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Return a G4VSolid modeling the Capsule
G4LogicalVolume* ASGARD::ConstructCapsule(){

  G4int nbslice = 7;
  G4double zSlice[7] = {  0.0*mm,
			  CapsuleWidth-0.1*mm,
			  CapsuleWidth,
			  CapsuleEdgeDepth,
			  CapsuleLength-CapsuleWidth,
			  CapsuleLength-CapsuleWidth-0.1*mm,
			  CapsuleLength  };
   
  G4double InnNullRad[7] = {0,0,0,0,0,0,0};

  G4double InnRad[7] = {  0.0*mm,
			  0.0*mm,
			  widthface-CapsuleWidth,
			  CrystalEdgeOffset1 + CrystalEdgeOffset2 + CrystalToCapsule - CapsuleWidth,
			  CrystalEdgeOffset1 + CrystalEdgeOffset2 + CrystalToCapsule - CapsuleWidth,
			  0.0*mm,
			  0.0*mm};
 
  G4double OutRad[7] = {  widthface-1.5*mm,
			  widthface,
			  widthface,
			  CrystalEdgeOffset1 + CrystalEdgeOffset2 + CrystalToCapsule,
			  CrystalEdgeOffset1 + CrystalEdgeOffset2 + CrystalToCapsule,
			  CrystalEdgeOffset1 + CrystalEdgeOffset2 + CrystalToCapsule,
			  CrystalEdgeOffset1 + CrystalEdgeOffset2 + CrystalToCapsule};

  // The whole volume of the Capsule, made of N2

  G4Polyhedra* caps = new G4Polyhedra(G4String("Capsule"), 45.*deg, 360.*deg, 4, nbslice, zSlice, InnNullRad, OutRad);
  G4LogicalVolume* LogicCapsule=
    new G4LogicalVolume(caps,m_MaterialN2,"LogicCapsule", 0, 0, 0);
  LogicCapsule->SetVisAttributes(G4VisAttributes::Invisible);

  // The wall of the Capsule made of Al
  G4Polyhedra* capsWall = new G4Polyhedra(G4String("CapsuleWall"), 45.*deg, 360.*deg, 4, nbslice, zSlice, InnRad, OutRad);
  G4LogicalVolume* logicCapsuleWall =
    new G4LogicalVolume(capsWall,m_MaterialAl,"LogicCapsuleWall", 0, 0, 0);

  new G4PVPlacement(G4Transform3D(*(new G4RotationMatrix()), G4ThreeVector(0,0,0)),
		    logicCapsuleWall,"CapsuleWall",LogicCapsule,false,1);

  logicCapsuleWall->SetVisAttributes(TrGreyVisAtt);

  return LogicCapsule;

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* ASGARD::ConstructDewar(){
  G4Tubs* DewarSolid = new G4Tubs("DewarSolid",0,90*mm*0.5,90*mm*0.5,0,2*M_PI);
  G4Tubs* DewarCFSolid = new G4Tubs("DewarCFSolid",0,45*mm*0.5,145*mm*0.5,0,2*M_PI);

  G4UnionSolid* DewarFull =
    new G4UnionSolid("Dewarfull", DewarSolid, DewarCFSolid, new G4RotationMatrix(),G4ThreeVector(0,0,-90*mm-(145-90)*0.5*mm));

  G4LogicalVolume* LogicDewar = new G4LogicalVolume(DewarFull,m_MaterialAl,"LogicDewar",0,0,0);

  G4Tubs* N2Solid = new G4Tubs("N2Solid",0,90*mm*0.5-1*mm,90*mm*0.5-1*mm,0,2*M_PI);
  G4Tubs* N2CFSolid = new G4Tubs("N2CFSolid",0,45*mm*0.5-1*mm,145*mm*0.5-1*mm,0,2*M_PI);

  G4LogicalVolume* LogicN2 = new G4LogicalVolume(N2Solid,m_MaterialN2,"LogicN2",0,0,0);
  G4LogicalVolume* LogicN2CF = new G4LogicalVolume(N2CFSolid,m_MaterialN2,"LogicN2CF",0,0,0);
 
  LogicN2->SetVisAttributes(GreenVisAtt);
  LogicN2CF->SetVisAttributes(GreenVisAtt);
  new G4PVPlacement(G4Transform3D(*(new G4RotationMatrix()), G4ThreeVector(0,0,0)),
		    LogicN2,"N2 Deware",LogicDewar,false,1);
 
  // new G4PVPlacement(G4Transform3D(*(new G4RotationMatrix()), G4ThreeVector(0,0,-90*mm)),
  // LogicN2CF,"N2 Deware",LogicDewar,false,1);
  new G4PVPlacement(G4Transform3D(*(new G4RotationMatrix()), G4ThreeVector(0,0,-120*mm)),
		    LogicN2CF,"N2 Deware",LogicDewar,false,1);

  LogicDewar->SetVisAttributes(TrGreyVisAtt);

  return LogicDewar;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Construct clover detector
// Called After DetecorConstruction::AddDetector Method
void ASGARD::ConstructClover(){
  if(m_LogicClover==0){
    // Construct the clover itself
    m_LogicClover = ConstructCapsule();

    double CrystalOffset = (CrystalReducedRadius+0.5*mm);
    G4LogicalVolume* logicCrystal = ConstructCrystal();
    logicCrystal->SetVisAttributes(TrGreyVisAtt);
  
    G4LogicalVolume* logicSegment[8];
    for(int i = 0; i < 8; i++){
      logicSegment[i] = ConstructSegment(i);
  
    }
    logicSegment[0]->SetVisAttributes(RedVisAtt);
    logicSegment[1]->SetVisAttributes(BlueVisAtt);
    logicSegment[2]->SetVisAttributes(GreenVisAtt);
    logicSegment[3]->SetVisAttributes(WhiteVisAtt);
    logicSegment[4]->SetVisAttributes(TestVisAtt1);
    logicSegment[5]->SetVisAttributes(TestVisAtt2);
    logicSegment[6]->SetVisAttributes(TestVisAtt3);
    logicSegment[7]->SetVisAttributes(WhiteVisAtt);
  
  
    G4RotationMatrix* CrystalRotationB = new G4RotationMatrix(0,0,0); // no rotation
    G4ThreeVector* SegmentPosition = new G4ThreeVector(0., 0., 0.);
    // place segments in crystal
    for(int i = 0; i < 8; i++){
      new G4PVPlacement(G4Transform3D(*CrystalRotationB, *SegmentPosition),
			logicSegment[i],Form("LogicSegment%d", i),logicCrystal,false,i+1); // start from 1 for indexing (0 will be reserved for core)
    }
  
    // place 4 crystals in clover
    G4ThreeVector CrystalPositionB = G4ThreeVector(-CrystalOffset,+CrystalOffset,0.5*CrystalLength+7*mm);
    new G4PVPlacement(G4Transform3D(*CrystalRotationB, CrystalPositionB),
		      logicCrystal,"LogicCrystalB",m_LogicClover,true,0);
  
    G4RotationMatrix* CrystalRotationG = new G4RotationMatrix(0,0,0);
    CrystalRotationG->rotate(-90*deg, G4ThreeVector(0,0,1));
    G4ThreeVector CrystalPositionG = G4ThreeVector(+CrystalOffset,+CrystalOffset,0.5*CrystalLength+7*mm);
    new G4PVPlacement(G4Transform3D(*CrystalRotationG, CrystalPositionG),
		      logicCrystal,"LogicCrystalG",m_LogicClover,true,1);
  
    G4RotationMatrix* CrystalRotationR = new G4RotationMatrix(0,0,0);
    CrystalRotationR->rotate(-180*deg, G4ThreeVector(0,0,1));
    G4ThreeVector CrystalPositionR = G4ThreeVector(+CrystalOffset,-CrystalOffset,0.5*CrystalLength+7*mm);
    new G4PVPlacement(G4Transform3D(*CrystalRotationR, CrystalPositionR),
		      logicCrystal,"LogicCrystalR",m_LogicClover,true,2);
  
    G4RotationMatrix* CrystalRotationW = new G4RotationMatrix(0,0,0);
    CrystalRotationW->rotate(90*deg, G4ThreeVector(0,0,1));
    G4ThreeVector CrystalPositionW = G4ThreeVector(-CrystalOffset,-CrystalOffset,0.5*CrystalLength+7*mm);
    new G4PVPlacement(G4Transform3D(*CrystalRotationW, CrystalPositionW),
		      logicCrystal,"LogicCrystalW",m_LogicClover,true,3);
    
  }
    
}

void ASGARD::ConstructDetector(G4LogicalVolume* world){
  ConstructClover();

  G4RotationMatrix* DetectorRotation = new G4RotationMatrix(0,0,0);
  for (unsigned int i = 0 ;  i < m_CloverId.size(); i++) {

    // Constructing the Detector referential and the transition matrix
    G4ThreeVector U,V,W;
    G4double wX = sin(m_Theta[i]) * cos(m_Phi[i]) ;
    G4double wY = sin(m_Theta[i]) * sin(m_Phi[i]) ;
    G4double wZ = cos(m_Theta[i]);
    W = G4ThreeVector(wX, wY, wZ) ;

    // vector parallel to one axis of the entrance plane
    G4double vX = cos(m_Theta[i]) * cos(m_Phi[i]);
    G4double vY = cos(m_Theta[i]) * sin(m_Phi[i]);
    G4double vZ = -sin(m_Theta[i]);
    V = G4ThreeVector(vX, vY, vZ);

    W = W.unit();
    U = V.cross(W);
    U = U.unit();
    V = W.cross(U);
    V = V.unit();
    // Passage Matrix from Lab Referential to Clover Referential
    delete DetectorRotation;
    DetectorRotation = new G4RotationMatrix(U, V, W);

    DetectorRotation->rotate(m_BetaX[i], U);
    DetectorRotation->rotate(m_BetaY[i], V);
    DetectorRotation->rotate(m_BetaZ[i], W);
    G4ThreeVector DetectorPosition = m_R[i]*W;
  
    new G4PVPlacement(G4Transform3D(*DetectorRotation, DetectorPosition),
		      m_LogicClover,"Clover",world,false,m_CloverId[i]);
  
    // G4LogicalVolume* LogicDewar = ConstructDewar();
  
    // // new G4PVPlacement(G4Transform3D(*DetectorRotation, DetectorPosition+W*((90*mm+(145)*mm)+CapsuleLength*0.5+90*0.5*mm)),
    // //     LogicDewar,"Dewar",world,false,m_CloverId[i]);
    // new G4PVPlacement(G4Transform3D(*DetectorRotation, DetectorPosition+W*((90*mm+(145)*mm)+CapsuleLength*0.5+20*0.5*mm)),
    // 		      LogicDewar,"Dewar",world,false,m_CloverId[i]);

    
    // Determine CoG coordinates for each detector element according to m_R, m_Theta, m_Phi
    
    
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void ASGARD::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("ASGARD")){
    pTree->Branch("ASGARD", "TASGARDData", &m_ASGARDData) ;
  }
  pTree->SetBranchAddress("ASGARD", &m_ASGARDData) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void ASGARD::ReadSensitive(const G4Event* event){
  m_ASGARDData->Clear();

  ///////////
  // HPGE
  NPS::HitsMap<G4double*>*     HPGEHitMap;
  std::map<G4int, G4double**>::iterator    HPGE_itr;
  //  G4cout<<"*************"<<event->GetEventID()<<"*************"<<G4endl;
  G4int HPGECollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("ASGARD_Scorer/ASGARD");
  if(HPGECollectionID == -1) {
    G4cerr << "ERROR: No Collection found for ASGARDScorer: Skipping processing of ASGARD Hit" << G4endl;
    return;
  }
  HPGEHitMap = (NPS::HitsMap<G4double*>*)(event->GetHCofThisEvent()->GetHC(HPGECollectionID));

  // Loop on the HPGE map
  for (HPGE_itr = HPGEHitMap->GetMap()->begin() ; HPGE_itr != HPGEHitMap->GetMap()->end() ; HPGE_itr++){
    G4double* Info = *(HPGE_itr->second);
    G4int m_index = HPGE_itr->first; // identifier for clover/crystal/segment
    G4double energy = Info[0]/keV;
    //    G4double Energy_reso = 0.4*keV + Energy_true*0.0015; // FWHM
    
    G4double Energy_reso = (1.5 + energy*0.0015)/2.355; // sigma: convert FWHM of (1.5 keV + E(keV)*0.15%) by division of 2.355
    
    //    if (energy > 30*keV && energy < 110*keV)

    // see ASGARDScorers.cc for clover/crystal/segment numbering
    G4int clo = m_index/100;
    G4int cry = (m_index%100)/10;
    G4int seg = m_index%10;
    
    //    G4cout<<clo<<" "<<cry<<" "<<seg<<" "<<energy<<G4endl;
    G4double time     =  Info[1]; // global detection time in ns (including beta-decay time)
    
    // need to assess segment hit info for core energy, highest-energy segment, etc.
      
    G4double theta_semitrue = Info[5];
    G4double phi_semitrue = Info[6];


    TVector3 v_pos(Info[2], Info[3], Info[4]);

    
    //
    if (energy > 10*keV){ // JP: increased threshold from 1*keV
      m_ASGARDData->AddGeCloverNbr(clo);
      m_ASGARDData->AddGeCrystalNbr(cry);
      m_ASGARDData->AddGeSegmentNbr(seg);
      m_ASGARDData->AddGeEnergy(energy);
      m_ASGARDData->AddGeTimeLED(RandGauss::shoot(time, 15));
      m_ASGARDData->AddGePosition(v_pos);
      m_ASGARDData->AddGeThetaSemiTrue(theta_semitrue);
      m_ASGARDData->AddGePhiSemiTrue(phi_semitrue);
      // m_ASGARDData->SetGeSegmentEnergy(RandGauss::shoot(Energy, Energy_reso));
      //      m_ASGARDData->AddGeEnergy(RandGauss::shoot(energy, Energy_reso));
    }
  }
  m_ASGARDData->ResetLastCloverCrystal();
  m_ASGARDData->SortSegmentData();
  // clear map for next event
  HPGEHitMap->clear();
  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void ASGARD::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_ASGARDScorer = CheckScorer("ASGARD_Scorer",already_exist) ;

  if(already_exist) 
    return ;

  //   ASGARD Associate Scorer
  G4VPrimitiveScorer* ASGARDScorer = new ASGARDSCORERS::PS_ASGARD("ASGARD",0);

  //and register it to the multifunctionnal detector
  m_ASGARDScorer->RegisterPrimitive(ASGARDScorer);

  // Otherwise the scorer is initialised
  //and register it to the multifunctionnal detector
  G4SDManager::GetSDMpointer()->AddNewDetector(m_ASGARDScorer) ;
}
void ASGARD::InitializeMaterial(){
  // JP: Trying to build Auger transition table for each material, it seg. faults for "vacuum"!
  //  m_MaterialVacuum = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
  m_MaterialGe= MaterialManager::getInstance()->GetMaterialFromLibrary("Ge");
  m_MaterialAl= MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
  m_MaterialN2= MaterialManager::getInstance()->GetMaterialFromLibrary("N2_liquid"); 
}

// void ASGARD::DefineSegmentPositions(){
//   // coordinates of segments of  crystal C
//   v_xyz_segment[0][0] = TVector3(23.4, 27.2, -57);
//   v_xyz_segment[0][1] = TVector3(37.4, 37.2, -29);
//   v_xyz_segment[0][2] = TVector3(13.3, 37.5, -29.2);
//   v_xyz_segment[0][3] = TVector3(13.6, 13.2, -28.9);
//   v_xyz_segment[0][4] = TVector3(36.7, 14.0, -29.3);
//   v_xyz_segment[0][5] = TVector3(35.4, 40.9, -72.6);
//   v_xyz_segment[0][6] = TVector3(9.9, 40.7, -72.6);
//   v_xyz_segment[0][7] = TVector3(9.6, 15.1, -72.6);
//   v_xyz_segment[0][8] = TVector3(35.7, 14.9, -72.6);
// }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* ASGARD::Construct(){
  return  (NPS::VDetector*) new ASGARD();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_ASGARD{
  public:
    proxy_nps_ASGARD(){
      NPS::DetectorFactory::getInstance()->AddToken("ASGARD","ASGARD");
      NPS::DetectorFactory::getInstance()->AddDetector("ASGARD",ASGARD::Construct);
    }
  };

  proxy_nps_ASGARD p_nps_ASGARD;
}
