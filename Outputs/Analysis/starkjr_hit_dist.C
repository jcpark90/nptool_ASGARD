TChain* tc;
const double EBEAM = 40.*3.5; // 3.5 MeV/u
const double EX = 1.461; // 1.461-MeV excitation energy for 40Ar
const double AMU = 931.494;
const double MASS[2] = {AMU*40, AMU*197}; // [0]: beam (40Ar), [1]: target (197Au) 
TH1D* hraw_asgard;
TH1D* hdop_asgard_40Ar;
TH1D* hdop_asgard_197Au;
TH1D* hbeta;
TH2D* meg_thpg_40Ar;
TH2D* meg_thpg_197Au;
TH2D* ep_z;
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include "TASGARDData.h"
#include "TSTARKData.h"
TTreeReader fReader;
TH2D* hxy[6][8];
TH1D* hz[6][8];

const double DIST = 110.;
int BeamOrTarget(double ebeam, double ep, double z){ // return 0 for beam, 1 for target based on kinematics
  if (z < 0)
    return 0; // target cannot be scattered backward
  else{
    if (ep > z) // near y = x relationship for 3.7 MeV/u 40Ar beam on 2.36 um 197AU target
      return 0;
    else
      return 1;
  }
}
double GetThetaT( double epb, double thb) {

  /// Returns theta angle of target nucleus using angle and energy of beam
	
  double tau = MASS[0]/MASS[1];
  double Eprime = EBEAM - EX*(1+tau);
  double epsilon = TMath::Sqrt(EBEAM/Eprime);
  double x, y, TTh;
  if( tau > 1 ) { // inverse kinematics: maximum scattering angle may be exceeded...
    y = TMath::ASin(1./(tau*epsilon)); // maximum projectile angle in lab
    if( thb < y ) y = thb;
    y = TMath::Tan(y);
  } else {
    y = TMath::Tan(thb); // y = tan(Theta_projlab)
  }
  // if( tau > 1 && rand.Gaus(epb,30000.)<50000. ) {
  //   x = (y*y*epsilon*tau + TMath::Sqrt(-y*y*epsilon*epsilon*tau*tau + y*y + 1) ) / (1+y*y);
  // }
  // else {
    x = (y*y*epsilon*tau - TMath::Sqrt(-y*y*epsilon*epsilon*tau*tau + y*y + 1) ) / (1+y*y);
    //  }
  //	cout << "Centre of mass angle: " << TMath::ACos(x)*TMath::RadToDeg() << endl;
  TTh = TMath::ATan( TMath::Sqrt( 1 - x*x ) / ( epsilon + x ) ); // choose kinematic branch using energy cut... as I haven't a clue how to do it any other way?!
  if( TTh < 0 ) TTh += TMath::Pi();
  //	cout << "Simulated target angle: " << TTh*TMath::RadToDeg() << endl;
  return TTh;
}
double GetThetaB( double tht ) {

	/// Returns theta angle of B using angle of T
	
	double tau = MASS[0]/MASS[1];
	double Eprime = EBEAM - EX*(1+tau);
	double epsilon = TMath::Sqrt(EBEAM/Eprime);
	double x, y, BTh;
	y = TMath::Tan(tht); // y = tan(Theta_targetlab)
	x = (y*y*epsilon - TMath::Sqrt(-y*y*epsilon*epsilon + y*y + 1) ) / (1+y*y);
	//	cout << "Centre of mass angle: " << TMath::ACos(x)*TMath::RadToDeg() << endl;
	BTh = TMath::ATan( TMath::Sqrt( 1 - x*x ) / ( tau*epsilon + x ) );
	if( BTh < 0 ) BTh += TMath::Pi();
	//	cout << "Simulated beam angle: " << BTh*TMath::RadToDeg() << endl;
	return BTh;

}

double GetBeta(double ekin, double mass){ // in MeV and MeV/c2)
  double gamma = 1.+ekin/mass;
  return sqrt(1.-pow(gamma, -2.));
}
double GetCosThetaPG(double thp, double php, double thg, double phg){
  return sin(thp)*sin(thg)*cos(php-phg)+cos(thp)*cos(thg);
}
double GetDopplerGamma(double eg, double beta, double thp, double php, double thg, double phg){
  double cos_theta_pg = GetCosThetaPG(thp, php, thg, phg);
  return eg*(1.-beta*cos_theta_pg)/sqrt(1.-beta*beta);
}
void DefineHistos(){
  for(int i = 0; i < 6; i++){
    for(int j = 0; j < 8; j++){
      hxy[i][j] = new TH2D(Form("hxy_%d_%d", i+1, j+1), "", 200, -100, 100, 200, -100, 100);
      hz[i][j] = new TH1D(Form("hz_%d_%d", i+1, j+1), "", 200, -100, 100);
    }
  }
}
void starkjr_hit_dist(){
  gStyle->SetOptStat(11111111);
  tc = new TChain("SimulatedTree");
  tc->Add("../Simulation/40Ar_197Au_coulex_real_target_3mm_beam_140mm_asgard.root");
  fReader.SetTree(tc);
  DefineHistos();
  TTreeReaderValue<TSTARKData> STARK = {fReader, "STARK"};
  while(fReader.Next()){           
    // i++;
    // if((i % 100000)==0) cout << i << endl;
    //    cout<<STARK->GetMult()<<endl;
    for(int n = 0; n < STARK->GetMult(); n++){
      TVector3 vec = STARK->GetPos(n);
      //      cout<<vec.X()<<" "<<vec.Y()<<" "<<vec.Z()<<endl;
      int det = STARK->GetDetN(n); // 1 to 6
      int fst = STARK->GetFStN(n); // 1 to 8
      int bst = STARK->GetBStN(n);
      
      if (det==1 && STARK->GetFrE(n) > 5){
	hxy[det-1][fst-1]->Fill(vec.X(), vec.Y());
	hz[det-1][fst-1]->Fill(vec.Z());
      }
      // if (clo==5 && STARKy->GetGeEnergy(n) > 100){ // clo 4: phi==0
      // 	hxy[cry][seg]->Fill(vec.X()-DIST/sqrt(2), vec.Y()-DIST/sqrt(2));
      // 	hz[cry][seg]->Fill(vec.Z());
      // 	hxy[cry][0]->Fill(vec.X()-DIST/sqrt(2), vec.Y()-DIST/sqrt(2));
      // 	hz[cry][0]->Fill(vec.Z());
      // }
    }
  }
  TCanvas* cxy[4];
  TCanvas* cz[4];
  for (int i = 0; i < 1; i++){
    cxy[i] = new TCanvas(Form("cxy%d", i+1), Form("cxy%d", i+1), 1600, 800);
    cxy[i]->Divide(4,2);
    cz[i] = new TCanvas(Form("cz%d", i+1), Form("cz%d", i+1), 1600, 800);
    cz[i]->Divide(4,2);
    for(int j = 0; j < 8; j++){
      cxy[i]->cd(j+1);
      hxy[i][j]->Draw("colz");
      cz[i]->cd(j+1);
      hz[i][j]->Draw("hist");
      cout<<i<<" "<<j<<" "<<hxy[i][j]->GetMean(1)<<" "<<hxy[i][j]->GetMean(2)<<" "<<hz[i][j]->GetMean()<<endl;
    }
  }
}
