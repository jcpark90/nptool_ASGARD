/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : January 2009                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This event Generator is used to simulated Isotropic ion Source           *
 *  Very usefull to figure out Geometric Efficacity of experimental Set-Up   *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
// C++
#include<limits>

// G4 headers
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ThreeVector.hh"

// G4 headers including CLHEP headers
// for generating random numbers
#include "Randomize.hh"

// NPS headers
#include "EventGeneratorIsotropic.hh"

// NPL headers
#include "RootOutput.h"
#include "NPNucleus.h"
#include "NPOptionManager.h"
using namespace CLHEP;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
EventGeneratorIsotropic::EventGeneratorIsotropic(){
  m_ParticleStack = ParticleStack::getInstance();
  event_ID=0;
  atomic_CS=-1;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
EventGeneratorIsotropic::~EventGeneratorIsotropic(){
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
EventGeneratorIsotropic::SourceParameters::SourceParameters(){
  m_EnergyLow    =  0  ;
  m_EnergyHigh   =  0  ;
  m_EnergyDistribution = "flat";
  m_x0           =  0  ;
  m_y0           =  0  ;
  m_z0           =  0  ;
  m_source_theta =  0  ;
  m_source_phi   =  0  ;
  m_SigmaX       =  0  ;
  m_SigmaY       =  0  ;
  m_SigmaR       =  0  ;
  m_SigmaZ       =  0  ;
  m_atomic_bgd_file = "";
  m_si_det_thickness = 0  ;
  m_particle     = NULL;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventGeneratorIsotropic::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Isotropic");
  m_Parameters.reserve(blocks.size());
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << endl << "\033[1;35m//// Isotropic reaction found " << endl; 

  vector<string> token = {"EnergyLow","EnergyHigh","HalfOpenAngleMin","HalfOpenAngleMax","x0","y0","z0","Particle"};
  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(token)){
      m_Parameters.push_back(SourceParameters());
      const vector<SourceParameters>::reverse_iterator it = m_Parameters.rbegin();

      it->m_EnergyLow         =blocks[i]->GetDouble("EnergyLow","MeV");
      it->m_EnergyHigh        =blocks[i]->GetDouble("EnergyHigh","MeV");
      if(blocks[i]->HasToken("EnergyDistribution"))
        it->m_EnergyDistribution=blocks[i]->GetString("EnergyDistribution");
      it->m_HalfOpenAngleMin  =blocks[i]->GetDouble("HalfOpenAngleMin","deg");
      it->m_HalfOpenAngleMax  =blocks[i]->GetDouble("HalfOpenAngleMax","deg");
      it->m_x0                =blocks[i]->GetDouble("x0","mm");
      it->m_y0                =blocks[i]->GetDouble("y0","mm");
      it->m_z0                =blocks[i]->GetDouble("z0","mm");
      vector<string> particleName =blocks[i]->GetVectorString("Particle");
      for(unsigned int j = 0 ; j < particleName.size() ; j++){
        if(particleName[j]=="proton"){ it->m_particleName.push_back("1H")  ;} 
        else if(particleName[j]=="deuton"){ it->m_particleName.push_back("2H")  ; }
        else if(particleName[j]=="triton"){ it->m_particleName.push_back("3H")  ; }
        else if(particleName[j]=="3He" || particleName[j]=="He3") { it->m_particleName.push_back("3He") ; }
        else if(particleName[j]=="alpha") { it->m_particleName.push_back("4He") ; }
        else if(particleName[j]=="gamma") { it->m_particleName.push_back("gamma") ;}
        else if(particleName[j]=="mu+") { it->m_particleName.push_back("mu+") ;}
        else if(particleName[j]=="mu-") { it->m_particleName.push_back("mu-") ;}
        else if(particleName[j]=="neutron") {it->m_particleName.push_back("neutron") ;}
        else it->m_particleName.push_back(particleName[j]);
      }

      if(blocks[i]->HasToken("ExcitationEnergy"))
        it->m_ExcitationEnergy =blocks[i]->GetVectorDouble("ExcitationEnergy","MeV");

      if(blocks[i]->HasToken("SigmaX"))
        it->m_SigmaX=blocks[i]->GetDouble("SigmaX","mm");
      if(blocks[i]->HasToken("SigmaY"))
        it->m_SigmaY=blocks[i]->GetDouble("SigmaY","mm");
      if(blocks[i]->HasToken("SigmaZ"))
        it->m_SigmaZ=blocks[i]->GetDouble("SigmaZ","mm");      
      if(blocks[i]->HasToken("SigmaR"))
        it->m_SigmaR=blocks[i]->GetDouble("SigmaR","mm");
      if(blocks[i]->HasToken("source_theta"))
        it->m_source_theta=blocks[i]->GetDouble("source_theta","deg");
      if(blocks[i]->HasToken("source_phi"))
        it->m_source_phi=blocks[i]->GetDouble("source_phi","deg");
      if(blocks[i]->HasToken("atomic_bgd_file"))
        it->m_atomic_bgd_file=blocks[i]->GetString("atomic_bgd_file");
      if(blocks[i]->HasToken("si_det_thickness"))
        it->m_si_det_thickness=blocks[i]->GetDouble("si_det_thickness","mm");
      if(blocks[i]->HasToken("Multiplicity"))
        it->m_Multiplicty=blocks[i]->GetVectorInt("Multiplicity");
    }
    else{
      cout << "ERROR: check your input file formatting \033[0m" << endl; 
      exit(1);
    }
  }
  for(auto& par : m_Parameters) {
    if(par.m_ExcitationEnergy.size()==0)
      par.m_ExcitationEnergy.resize(par.m_particleName.size(),0);
    if(par.m_Multiplicty.size()==0)
      par.m_Multiplicty.resize(par.m_particleName.size(),1);

    if(par.m_EnergyDistribution!="flat"){
      if(par.m_EnergyDistribution=="Watt"){
        fEnergyDist = new TF1("fWatt","0.4865*TMath::SinH(sqrt(2*x))*TMath::Exp(-x)",par.m_EnergyLow,par.m_EnergyHigh);
      }
      else{
        fEnergyDist = new TF1("fDist", par.m_EnergyDistribution, par.m_EnergyLow, par.m_EnergyHigh);
      }

    }
    if (par.m_atomic_bgd_file!=""){
      SetBackgroundFile(par.m_atomic_bgd_file, par.m_si_det_thickness);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventGeneratorIsotropic::GenerateEvent(G4Event*){

  for(auto& par : m_Parameters) {
    for(unsigned int p=0; p<par.m_particleName.size(); p++){
      for(int i=0; i<par.m_Multiplicty[p]; i++){
        par.m_particle=NULL;
        if(par.m_particle==NULL){

          if(par.m_particleName[p]=="gamma" || par.m_particleName[p]=="neutron" ||  par.m_particleName[p]=="opticalphoton"  ||  par.m_particleName[p]=="mu+" ||  par.m_particleName[p]=="mu-"){
            par.m_particle =  G4ParticleTable::GetParticleTable()->FindParticle(par.m_particleName[p].c_str());
          }
          else{
            NPL::Nucleus* N = new NPL::Nucleus(par.m_particleName[p]);
            par.m_particle = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon(N->GetZ(), N->GetA(),par.m_ExcitationEnergy[p]);
            delete N;
          }

        }

        G4double cos_theta_min   = cos(par.m_HalfOpenAngleMin);
        G4double cos_theta_max   = cos(par.m_HalfOpenAngleMax);
        G4double cos_theta       = cos_theta_min + (cos_theta_max - cos_theta_min) * RandFlat::shoot();
        G4double theta           = acos(cos_theta)                                                   ;
        G4double phi             = RandFlat::shoot() * 2 * pi                                        ;
        G4double particle_energy;
        if(par.m_EnergyDistribution=="flat"){
          particle_energy = par.m_EnergyLow + RandFlat::shoot() * (par.m_EnergyHigh - par.m_EnergyLow)    ;
          event_ID++;
        }
        else if(par.m_EnergyDistribution=="Watt"){
          particle_energy = fEnergyDist->GetRandom();
        }
        else{
          particle_energy = fEnergyDist->GetRandom();
        }


        // Direction of particle, energy and laboratory angle
        G4double momentum_x = sin(theta) * cos(phi)  ;
        G4double momentum_y = sin(theta) * sin(phi)  ;
        G4double momentum_z = cos(theta)             ;

        G4double x0;
        G4double y0;
	G4double z0; 
	if (par.m_SigmaX > 0)
	  x0 = RandGauss::shoot(par.m_x0,par.m_SigmaX);
	else
	  x0 = par.m_x0+(RandFlat::shoot()-0.5)*2*par.m_SigmaX; 
	if (par.m_SigmaY > 0)
	  y0 = RandGauss::shoot(par.m_y0,par.m_SigmaY);
	else
	  y0 = par.m_y0+(RandFlat::shoot()-0.5)*2*par.m_SigmaY;
       	if (par.m_SigmaZ > 0)
	  z0 = RandGauss::shoot(par.m_z0,par.m_SigmaZ);
	else
	  z0 = par.m_z0+(RandFlat::shoot()-0.5)*2*par.m_SigmaZ;
	if (par.m_SigmaR > 0){
	  double r = abs(RandGauss::shoot(0, par.m_SigmaR));
	  
	  double th = RandFlat::shoot()* 2 * pi;
	  x0 = par.m_x0 + r*cos(th);
	  y0 = par.m_y0 + r*sin(th);
	}
	else if (par.m_SigmaR!=0){
	  double r = -par.m_SigmaR*sqrt(RandFlat::shoot());
	  double th = RandFlat::shoot()* 2 * pi;
	  x0 = par.m_x0 + r*cos(th);
	  y0 = par.m_y0 + r*sin(th);
	}
	G4ThreeVector sourcePos(x0, y0, z0);
	sourcePos.rotateX(-par.m_source_theta);
	sourcePos.rotateZ(par.m_source_phi-90*deg);

	//        Particle particle(par.m_particle, theta,particle_energy,G4ThreeVector(momentum_x, momentum_y, momentum_z),G4ThreeVector(x0, y0, par.m_z0));
	Particle particle(par.m_particle, theta,particle_energy,G4ThreeVector(momentum_x, momentum_y, momentum_z), sourcePos);

        m_ParticleStack->AddParticleToStack(particle);


	int numberofgammas = RandPoisson::shoot(gammasperevent);
	//	cout << "number of gammas for this event " << numberofgammas << endl;
	par.m_particle = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
	for(int gg = 0; gg<numberofgammas; gg++){
	  phi       = twopi*G4UniformRand();
	  bgdist->GetRandom2(theta,particle_energy);
	  theta *= pi/180;
	  //G4cout <<setprecision(5)<< theta << "\t" << phi << "\t" << energy << endl;
	  G4ThreeVector v(0.0,0.0,1.0);
	  v.setTheta(theta);
	  v.setPhi(phi);
	  // set the particle gun in the direction
	  Particle abgd(par.m_particle, theta,particle_energy*keV,v, sourcePos);
	  m_ParticleStack->AddParticleToStack(abgd);
	}
      }
    }
  }

}
void EventGeneratorIsotropic::SetBackgroundFile(G4String filename, G4double thickness){
  TFile *f = new TFile(filename);
  if(f->IsOpen()){
    bgdist = (TH2F*)f->Get("hSum");
    cout << "Background histogram title: " << bgdist->GetTitle() << endl;
    if(atomic_CS <0){
      atomic_CS =  extractCS(bgdist->GetTitle());
      cout << "atomic cross section " << atomic_CS << " barn" << endl;
      cout << "target thickness " << thickness << " mm " << endl;
      gammasperevent = Avogadro*2.329*thickness/10*atomic_CS*1e-24; // Si density: 2.329 g/cm3
      cout << "gammas per event " << gammasperevent << endl;
    }
  }
  else{
    G4cout << "couldn't open background file" << G4endl;
  }
}
G4double EventGeneratorIsotropic::extractCS(string str){
  vector<string> result;
  istringstream iss(str);
  for(string strr; iss >> strr; )
    result.push_back(strr);
  for(unsigned short i=0;i<result.size();i++){
    //cout << result[i] << "\t";
    double cs = sciToDouble(result[i]);
    //cout << cs << endl;
    if(cs>0)
      return cs;
      
  }
  return 0;
}
G4double EventGeneratorIsotropic::sciToDouble(string& str) {
   stringstream ss(str);
   double d = 0;
   ss >> d;
   //cout << d << endl;
   if (ss.fail()) {
     return -1;
   }
   return (d);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventGeneratorIsotropic::InitializeRootOutput(){

}
