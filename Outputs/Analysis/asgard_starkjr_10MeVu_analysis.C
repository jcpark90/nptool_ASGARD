#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include "TASGARDData.h"
#include "TSTARKData.h"
TChain* tc;

// ASGARD + STARK Jr geometry parameters
const bool EXACT_ASGARD_ANGLES = 0; // 0: get from segment CoG, 1: get true theta/phi from hit info
const bool EXACT_STARKJr_ANGLES = 0; // get discrete strip centroids, 1: get true theta/phi from hit info
const double DIST_ASGARD = 110.; // 110 mm from target center
const double X6_STRIP_WIDTH = 40.3/8;
const double X6_RHO = 47.2; // in mm, closest distance to center

// experiment beam/target/energy settings
const char* FILENAME = "asgard_starkjr_40Ar_10MeVu_pencil.root";
const double EBEAM = 40.*10; // in MeV, 40Ar at 10 MeV/u
const double EBEAM_LOSS_MAX = 120.; // maximum energy loss of beam particle for PID, in MeV
const double AMU = 931.494;
const int ZBEAM = 18;
const int ABEAM = 40;
const int ZTARG = 79;
const int ATARG = 197;
const double EG_BEAM_MAIN = 1461; // in keV
const double EG_TARG_MAIN = 548; // in keV

const double EX = (EG_BEAM_MAIN+EG_TARG_MAIN)/1000; // 1461 keV for 40Ar + 548 keV for 197Au, in MeV
const double MASS[2] = {AMU*40-39.040, AMU*197-31.140}; // [0]: beam (40Ar), [1]: target (197Au) in MeV
const double DENSITY_AU = 19.3e3; // in mg/cm3
const double TARGET_THICKNESS = 2.36e-4*DENSITY_AU; // physical thickness 2.36 um = 2.36e-4 cm, should be in in mg/cm2
const int NITER = 5; // number of calculation iterations for energy loss correction in target

// define histograms
TH1D* hraw_asgard;
TH1D* hdop_asgard_40Ar;
TH1D* hdop_asgard_197Au;
TH1D* hbeta;
TH2D* meg_thpg_40Ar;
TH2D* meg_thpg_197Au;
TH2D* ep_z_raw;
TH2D* ep_z_corrected;
TGraph* gSP[2]; // 0 for beam in target, 1 for target in target
TTreeReader fReader;

vector<double> vTHETA_ASGARD;
vector<double> vPHI_ASGARD;
vector<double> vDIST_ASGARD;
vector<double> vRHO_X6;
vector<double> vPHI_X6;
vector<bool> vREV_X6;
vector<TVector3> vXYZ_ASGARD[4][9];
vector<TVector3> vXYZ_X6[8];
// coordinate for Crystal A in clover (the rest will be derived from these coordinates)
TVector3 XYZ_ASGARD[4][9];
TVector3 XYZ_X6[8];

// catkin quantities
double mb; //mass of beam in MeV/c2
double mt; //mass of target in MeV/c2 
double me; //mass of ejectile in MeV/c2 
double mr; //mass of recoil in MeV/c2
double eb; // beam energy in MeV
double qf; // final state reaction Q-value, cell $O$10 in catkin excel file
// values to be determined on initialization
double y; // total mass of system
double betac;
double ecmi;
double ecmf;
double e3cm;
double beta3c;
double y_new;
double c;
// angle-dependent quantities
double cosagl;
double b; 
double a; 
double dsq; // d**2
double b3L1;
double y_loop;
double y_adj;
// calculated values
double cma; // center-of-mass angle in degrees
double el; // lab energy of the ejectile
double rcml; // cm/lab ratio
double ra; // recoil angle in degrees
double er; // lab energy of the recoil
double angl4; // recoil angle
// second solution quantities
double b3L2;
double y_loop2;
double y_adj2;
double angl4_2;
double el2;


void InitializeKinematics(double mb_in, double mt_in, double me_in, double mr_in, double eb_in, double qf_in){
  mb = mb_in;
  mt = mt_in;
  me = me_in;
  mr = mr_in;
  eb = eb_in;
  qf = qf_in;
  y = mb + mt; // total mass of system
  betac = sqrt(eb*(eb+2*mb))/(y + eb);
  ecmi = sqrt(y*y+(2*eb*mt));
  ecmf = ecmi+qf-y+me+mr;
  e3cm = (ecmf*ecmf+(me+mr)*(me-mr))/(2*ecmf);
  beta3c = sqrt(1.-1./pow(e3cm/me, 2.));
  y_new = pow(e3cm/me, 2.)*(1.-betac*betac);
  c = 1.-y_new;
}
void Calculate(double la){ // lab angle in degrees, ra also in degrees!
  cosagl = cos(la*TMath::DegToRad());
  b = -betac*cosagl;
  a = y_new+b*b;
  dsq = b*b-a*c;
  b3L1 = (-b+sqrt(dsq))/a;
  y_loop = (b3L1*cosagl-betac)/((1.-betac*b3L1*cosagl)*beta3c);
  if (fabs(y_loop)<=1.)
    y_adj = y_loop;
  else
    y_adj = 0.;

  cma = TMath::RadToDeg()*acos(y_adj);
  el = me*(1./sqrt(1.-b3L1*b3L1)-1);
  rcml = beta3c*(1+b/b3L1)/(1+betac*beta3c*y_adj)/b3L1;
  er = eb+qf-el;
  angl4 = TMath::RadToDeg()*asin(sqrt((el*(el+2.*me))/(er*(er+2*mr)))*sin(la*TMath::DegToRad()));
  if ((el*(el+2*me)*cosagl*cosagl > eb*(eb+2*mb)) && la < 90.)
    ra = 180.-angl4;
  else
    ra = angl4;
  // not used
  /*
    b3L2 = (-b-sqrt(dsq))/a;
    y_loop2 = (b3L2*cosagl-betac)/((1.-betac*b3L2*cosagl)*beta3c);
    if (fabs(y_loop2)<=1.)
    y_adj2 = y_loop2;
    else
    y_adj2 = 0.;
    if (b3L2 > 1e-7)
    el2 = me*(1./sqrt(1-B3L2*B3L2)-1);
    else
    el2 = 0.;
    if (el2 > 0.)
    angl4_2 = TMath::RadToDeg()*asin(sqrt((el2*(el2+2.*me))/(er*(er+2*mr)))*sin(la*TMath::DegToRad()));
  */
  //  cout<<"cosagl: "<<cosagl<<", b: "<<b<<", a: "<<a<<", dsq: "<<dsq<<", b3L1: "<<b3L1<<", y_loop: "<<y_loop<<"y_adj: "<<y_adj<<", cma: "<<cma<<", el: "<<el<<", rcml: "<<rcml<<", er: "<<er<<", angl4: "<<angl4<<", ra: "<<ra<<endl;
}
void DefineXYZ_X6(){
  // based on X6
  // with reference to phi = 0;
  for(int i = 0 ; i < 8; i++){
    XYZ_X6[i] = TVector3(0, X6_STRIP_WIDTH*(3.5-i), 0); // y from rho, z coordinate from true hit info
  }
  for(int i = 0; i < vRHO_X6.size(); i++){
    for(int j = 0; j < 8; j++){
      TVector3 vec_temp = XYZ_X6[j];
      vec_temp.SetXYZ(vRHO_X6[i]-0.5, vec_temp.Y(), 0.);
      if (vREV_X6[i])
	vec_temp.RotateX(TMath::DegToRad()*180.);
      vec_temp.RotateZ(TMath::DegToRad()*vPHI_X6[i]);
      vXYZ_X6[j].push_back(vec_temp);
      //      cout<<i<<" "<<j<<" "<<vec_temp.X()<<" "<<vec_temp.y()<<" "<<vec_temp.z()<<endl;
    }
  }
}
void DefineXYZ_ASGARD(){
  // from CAD drawing, can be adjusted
  // depth from 197Au and 40Ar gammas in trial Coulex
  // default coordinate system has z-axis as crystal depth
  // add distance R to z before rotating
  XYZ_ASGARD[0][0] = TVector3(-23.4, -27.2, -37);
  XYZ_ASGARD[0][1] = TVector3(-37.4, -37.2, -19);
  XYZ_ASGARD[0][2] = TVector3(-13.3, -37.5, -19.2);
  XYZ_ASGARD[0][3] = TVector3(-13.6, -13.2, -18.9);
  XYZ_ASGARD[0][4] = TVector3(-36.7, -14, -19.3);
  XYZ_ASGARD[0][5] = TVector3(-35.4, -40.9, -57.6);
  XYZ_ASGARD[0][6] = TVector3(-9.9, -40.7, -57.6);
  XYZ_ASGARD[0][7] = TVector3(-9.6, -15.1, -57.6);
  XYZ_ASGARD[0][8] = TVector3(-35.7, -14.9, -57.6);

  XYZ_ASGARD[1][0] = TVector3(-XYZ_ASGARD[0][0].X(), XYZ_ASGARD[0][0].Y(), XYZ_ASGARD[0][0].Z());
  XYZ_ASGARD[1][1] = TVector3(-XYZ_ASGARD[0][1].X(), XYZ_ASGARD[0][1].Y(), XYZ_ASGARD[0][1].Z());
  XYZ_ASGARD[1][2] = TVector3(-XYZ_ASGARD[0][4].X(), XYZ_ASGARD[0][4].Y(), XYZ_ASGARD[0][2].Z());
  XYZ_ASGARD[1][3] = TVector3(-XYZ_ASGARD[0][3].X(), XYZ_ASGARD[0][3].Y(), XYZ_ASGARD[0][3].Z());
  XYZ_ASGARD[1][4] = TVector3(-XYZ_ASGARD[0][2].X(), XYZ_ASGARD[0][2].Y(), XYZ_ASGARD[0][4].Z());
  XYZ_ASGARD[1][5] = TVector3(-XYZ_ASGARD[0][5].X(), XYZ_ASGARD[0][5].Y(), XYZ_ASGARD[0][5].Z());
  XYZ_ASGARD[1][6] = TVector3(-XYZ_ASGARD[0][8].X(), XYZ_ASGARD[0][8].Y(), XYZ_ASGARD[0][6].Z());
  XYZ_ASGARD[1][7] = TVector3(-XYZ_ASGARD[0][7].X(), XYZ_ASGARD[0][7].Y(), XYZ_ASGARD[0][7].Z());
  XYZ_ASGARD[1][8] = TVector3(-XYZ_ASGARD[0][6].X(), XYZ_ASGARD[0][6].Y(), XYZ_ASGARD[0][8].Z());

  XYZ_ASGARD[2][0] = TVector3(-XYZ_ASGARD[0][0].X(), -XYZ_ASGARD[0][0].Y(), XYZ_ASGARD[0][0].Z());
  XYZ_ASGARD[2][1] = TVector3(-XYZ_ASGARD[0][1].X(), -XYZ_ASGARD[0][1].Y(), XYZ_ASGARD[0][1].Z());
  XYZ_ASGARD[2][2] = TVector3(-XYZ_ASGARD[0][2].X(), -XYZ_ASGARD[0][2].Y(), XYZ_ASGARD[0][2].Z());
  XYZ_ASGARD[2][3] = TVector3(-XYZ_ASGARD[0][3].X(), -XYZ_ASGARD[0][3].Y(), XYZ_ASGARD[0][3].Z());
  XYZ_ASGARD[2][4] = TVector3(-XYZ_ASGARD[0][4].X(), -XYZ_ASGARD[0][4].Y(), XYZ_ASGARD[0][4].Z());
  XYZ_ASGARD[2][5] = TVector3(-XYZ_ASGARD[0][5].X(), -XYZ_ASGARD[0][5].Y(), XYZ_ASGARD[0][5].Z());
  XYZ_ASGARD[2][6] = TVector3(-XYZ_ASGARD[0][6].X(), -XYZ_ASGARD[0][6].Y(), XYZ_ASGARD[0][6].Z());
  XYZ_ASGARD[2][7] = TVector3(-XYZ_ASGARD[0][7].X(), -XYZ_ASGARD[0][7].Y(), XYZ_ASGARD[0][7].Z());
  XYZ_ASGARD[2][8] = TVector3(-XYZ_ASGARD[0][8].X(), -XYZ_ASGARD[0][8].Y(), XYZ_ASGARD[0][8].Z());

  XYZ_ASGARD[3][0] = TVector3(XYZ_ASGARD[0][0].X(), -XYZ_ASGARD[0][0].Y(), XYZ_ASGARD[0][0].Z());
  XYZ_ASGARD[3][1] = TVector3(XYZ_ASGARD[0][1].X(), -XYZ_ASGARD[0][1].Y(), XYZ_ASGARD[0][1].Z());
  XYZ_ASGARD[3][2] = TVector3(XYZ_ASGARD[0][4].X(), -XYZ_ASGARD[0][4].Y(), XYZ_ASGARD[0][2].Z());
  XYZ_ASGARD[3][3] = TVector3(XYZ_ASGARD[0][3].X(), -XYZ_ASGARD[0][3].Y(), XYZ_ASGARD[0][3].Z());
  XYZ_ASGARD[3][4] = TVector3(XYZ_ASGARD[0][2].X(), -XYZ_ASGARD[0][2].Y(), XYZ_ASGARD[0][4].Z());
  XYZ_ASGARD[3][5] = TVector3(XYZ_ASGARD[0][5].X(), -XYZ_ASGARD[0][5].Y(), XYZ_ASGARD[0][5].Z());
  XYZ_ASGARD[3][6] = TVector3(XYZ_ASGARD[0][8].X(), -XYZ_ASGARD[0][8].Y(), XYZ_ASGARD[0][6].Z());
  XYZ_ASGARD[3][7] = TVector3(XYZ_ASGARD[0][7].X(), -XYZ_ASGARD[0][7].Y(), XYZ_ASGARD[0][7].Z());
  XYZ_ASGARD[3][8] = TVector3(XYZ_ASGARD[0][6].X(), -XYZ_ASGARD[0][6].Y(), XYZ_ASGARD[0][8].Z());
  
  for(int n = 0; n < vTHETA_ASGARD.size(); n++){
    for(int i = 0; i < 4; i++){
      for(int j = 0; j < 9; j++){
	// rotate each core and segment according to theta, phi (and beta)
	TVector3 vec_temp = XYZ_ASGARD[i][j];
	vec_temp.SetXYZ(vec_temp.X(), vec_temp.Y(), vec_temp.Z()-vDIST_ASGARD[n]);
	vec_temp.RotateY(-TMath::DegToRad()*vTHETA_ASGARD[n]);
	vec_temp.RotateX(TMath::DegToRad()*90.);
	vec_temp.RotateZ(TMath::DegToRad()*vPHI_ASGARD[n]);
	//	cout<<n<<" "<<vTHETA_ASGARD[n]<<" "<<vPHI_ASGARD[n]<<" "<<i<<" "<<j<<": ("<<vec_temp.X()<<", "<<vec_temp.Y()<<", "<<vec_temp.Z()<<")"<<endl;
	vXYZ_ASGARD[i][j].push_back(vec_temp); // later call vXYZ_asgard[cry][seg].at(clo)
      }
    }
  }
  
}
string convertInt( int number ) {
	
  /// Convert an integer into a string
	
  stringstream ss;
  ss << number;
  return ss.str();
	
}


bool stoppingpowers(int pid, int zProj, int aProj, int zTarg, int aTarg) {
	 
  /// Open stopping power files and make TGraphs of data
  /// naming convention of files...
	
  string gElName[110] = {
			 "H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg",
			 "Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr",
			 "Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
			 "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
			 "In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd",
			 "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf",
			 "Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po",
			 "At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm",
			 "Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs",
			 "Mt","Ds" };

  string title = "Stopping powers for ";
  string srimfilename = "";

  // projectile
  srimfilename +=  convertInt(aProj) + gElName[zProj-1];
  title += convertInt(aProj) + gElName[zProj-1];
  // target
  srimfilename += "_in_" + convertInt(aTarg) + gElName[zTarg-1] + ".txt";
  title += " in " + convertInt(aTarg) + gElName[zTarg-1];

  title += ";Ion energy [MeV];Stopping power [MeV/(mg/cm^2)]";
		
	 
  ifstream infile;
  infile.open( srimfilename.c_str(), ios::in );
	 
  if( !infile.is_open() ) {
		  
    cout << "Cannot open " << srimfilename << endl;
    return false;
		  
  }
  cout<<title.c_str()<<endl;
  gSP[pid] = new TGraph();
  gSP[pid]->SetTitle( title.c_str() );

  string line, units, tmp_str;
  stringstream line_ss;
  bool endflag = false;
  double BEn, nucl, elec, total, tmp_dbl;
  int p = 0;
	 
  // Test file format
  getline( infile, line );
  if( line.substr( 0, 5 ) == "=====" ) { // space in front of === (2006) or not (2008)?
		  
    while( line.substr( 0, 5 ) != "-----" ) // Err, what if you use SRIM-2006, which doesn't have the spaces?
      getline( infile, line );
		  
    getline( infile, line ); // read first line of data
		  
  }
	 
  while( !infile.eof() && !endflag ) {
		  
    // Read in data
    line_ss.str("");
    //		cout<<line<<endl;
    line_ss << line;
    line_ss >> BEn >> units >> nucl >> elec >> tmp_dbl >> tmp_str >> tmp_dbl >> tmp_str;
		
    if( units == "eV" ) BEn *= 1E-6;
    else if( units == "keV" ) BEn *= 1E-3;
    else if( units == "MeV" ) BEn *= 1E0;
    else if( units == "GeV" ) BEn *= 1E3;
		
    total = nucl + elec ; // MeV / ( mg / cm^2 )
    //		cout<<pid<<" "<<BEn<<" "<<total<<endl;
    gSP[pid]->SetPoint( p, BEn, total );
		  
    // Get next line
    getline( infile, line );
    p++;
		  
    // If we've reached the end, stop
    if( line.substr( 0, 9 ) == "---------" ) endflag = true;
    if( line.substr( 0, 9 ) == " Multiply" ) endflag = true;
		  
  }
  infile.close();
  cout<<"stopping powers defined"<<endl;
  TCanvas *c = new TCanvas();
  gSP[pid]->Draw("A*");
  //gSP[pid]->GetXaxis()->SetTitleOffset(1.3);
  //gSP[pid]->GetYaxis()->SetTitleOffset(1.3);
  //TGaxis::SetMaxDigits(3);
  string pdfname = srimfilename.substr( 0, srimfilename.find_last_of(".") ) + ".pdf";
  c->SetLogx();
  c->SaveAs( pdfname.c_str() );
	 
  delete c;
	 
  return true;
	 
}
double GetRange( int pid, double Ei, double Ef){ // in mg/cm2
  if (Ef >= Ei)
    return 0.;
  else{
    double dxde = 0.;
    double dedx = 0.;
    int Nmeshpoints = 100; // number of steps to take in integration
    
    double de = (Ei-Ef)/(double)Nmeshpoints; // should be positive 
    double E = Ei;
    double range = 0.;
    
    for( int i = 0; i < Nmeshpoints; i++ ){
      
      if( E < 1. ) break; // when we fall below 1 MeV we assume maximum energy loss
      
      dedx = gSP[pid]->Eval(E);
      dxde = 1./dedx;
      
      range += dxde*de;
     
      E -= de;
	  
    }
    return range; // in mg/cm2
  }
	
  
}
double GetELoss( int pid, double Ei, double dist, int opt ) {

  /// Returns the energy loss at a given initial energy and distance travelled
  /// in the target, the contaminant layer or Si dead layer
  /// Ei is the initial energy in MeV, return value is also in MeV
  /// dist is the distance travelled in the target in mg/cm2
  // pid = 0 beam, 1 target
  /// opt = 0 calculates normal energy loss as particle moves through target (default), Ei = initial energy
  /// opt = 1 calculates energy increase, i.e. tracing particle back to reaction point, Ei = detected energy
  /// folder with the format 62Fe_109Ag.txt, 62Fe_Si.txt, 109Ag_109Ag.txt or 109Ag_Si.txt,
  /// for combo = "BT", "TT", "BS" and "TS", repsectively.
  /// The srim file should be in units of MeV/(mg/cm^2)
	
  double dedx = 0;
  int Nmeshpoints = 100; // number of steps to take in integration
  double dx = dist/(double)Nmeshpoints;
  double E = Ei;
	
  for( int i = 0; i < Nmeshpoints; i++ ){

    if( E < 1. ) break; // when we fall below 1 MeV we assume maximum energy loss

    dedx = gSP[pid]->Eval(E);
		
    if( opt == 1 )
      E += dedx*dx;
    else
      E -= dedx*dx;
		
  }
	
  //if( opt == 0 && combo == "BT" ) cout << "Eloss = " << Ei - E << endl;
	
  if( opt == 0 ) return Ei - E;
  else return E - Ei;

}
int BeamOrTarget(double ebeam, double ep, double theta){ // return 0 for beam, 1 for target based on kinematics
  // if (z < 0)
  //   return 0; // target cannot be scattered backward
  // else{
  //   if (ep > ebeam/z) // near y = x relationship for 3.7 MeV/u 40Ar beam on 2.36 um 197AU target
  //     return 0;
  //   else
  //     return 1;
  // }
  Calculate(theta*180/TMath::Pi()); // get theoretical energy of beam nucleus (assumption) for detected angle theta, "el"
  // neglect energy loss in target
  //  cout<<theta*180/TMath::Pi()<<" "<<ep<<" vs "<<el<<endl;
  if (el-ep < EBEAM_LOSS_MAX)
    return 0;
  else
    return 1;
  
  
}
double GetThetaT( double eb, double thb) {

  /// Returns theta angle of target nucleus using angle and energy of (reacted) beam
  // 
	
  double tau = MASS[0]/MASS[1];
  double Eprime = eb - EX*(1+tau);
  double epsilon = TMath::Sqrt(eb/Eprime);
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
double GetThetaB( double eb, double tht) {	/// Returns theta angle of B using angle of T if T detected
	
  double tau = MASS[0]/MASS[1];
  double Eprime = eb - EX*(1+tau);
  double epsilon = TMath::Sqrt(eb/Eprime);
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
  hraw_asgard = new TH1D("hraw_asgard", "", 1700, 0, 1700);
  hdop_asgard_40Ar = new TH1D("hdop_asgard_40Ar", "", 1700, 0, 1700);
  hdop_asgard_40Ar->GetXaxis()->SetTitle("Energy (keV)");
  hdop_asgard_40Ar->GetYaxis()->SetTitle("Counts / keV");
  hdop_asgard_197Au = new TH1D("hdop_asgard_197Au", "", 1700, 0, 1700);
  hdop_asgard_197Au->GetXaxis()->SetTitle("Energy (keV)");
  hdop_asgard_197Au->GetYaxis()->SetTitle("Counts / keV");
  hbeta = new TH1D("hbeta", "", 1000, 0, 0.1);
  meg_thpg_40Ar = new TH2D("meg_thpg_40Ar", "", 200, -1, 1, 1700, 0, 1700);
  meg_thpg_40Ar->GetXaxis()->SetTitle("cos(#theta_{p#gamma}) of ^{40}Ar");
  meg_thpg_40Ar->GetYaxis()->SetTitle("E_{#gamma} (keV)");
  meg_thpg_197Au = new TH2D("meg_thpg_197Au", "", 200, -1, 1, 1700, 0, 1700);
  meg_thpg_197Au->GetXaxis()->SetTitle("cos(#theta_{p#gamma}) of ^{197}Au");
  meg_thpg_197Au->GetYaxis()->SetTitle("E_{#gamma} (keV)");
  ep_z_raw = new TH2D("ep_z_raw", "", 200, -100, 100, EBEAM, 0, EBEAM);
  ep_z_raw->GetXaxis()->SetTitle("z from STARK Jr. (mm)");
  ep_z_raw->GetYaxis()->SetTitle("Detected particle energy (MeV)");
  ep_z_corrected = new TH2D("ep_z_corrected", "", 200, -100, 100, EBEAM, 0, EBEAM);
  ep_z_corrected->GetXaxis()->SetTitle("z from STARK Jr. (mm)");
  ep_z_corrected->GetYaxis()->SetTitle("Corrected particle energy (MeV)");
}
void AppendASGARDClover(double theta, double phi, double dist){
  vTHETA_ASGARD.push_back(theta);
  vPHI_ASGARD.push_back(phi);
  vDIST_ASGARD.push_back(dist);
}

void AppendX6(double rho, double phi, bool rev){
  vRHO_X6.push_back(rho);
  vPHI_X6.push_back(phi);
  vREV_X6.push_back(rev);
}
void ReadDetectorGeom(){
  // later read from detector geometry file
  // for now add trial
  // theta, phi, dist
  AppendASGARDClover(90, 180, DIST_ASGARD);
  AppendASGARDClover(90, 225, DIST_ASGARD);
  AppendASGARDClover(90, 135, DIST_ASGARD);
  AppendASGARDClover(90, -45, DIST_ASGARD);
  AppendASGARDClover(90, 0, DIST_ASGARD);
  AppendASGARDClover(90, 45, DIST_ASGARD);

  AppendX6(X6_RHO, 0, 1);
  AppendX6(X6_RHO, 60, 1);
  AppendX6(X6_RHO, 120, 1);
  AppendX6(X6_RHO, 180, 1);
  AppendX6(X6_RHO, 240, 1);
  AppendX6(X6_RHO, 300, 1);

  AppendX6(X6_RHO, 0, 0);
  AppendX6(X6_RHO, 60, 0);
  AppendX6(X6_RHO, 120, 0);
  AppendX6(X6_RHO, 180, 0);
  AppendX6(X6_RHO, 240, 0);
  AppendX6(X6_RHO, 300, 0);

}
void asgard_starkjr_10MeVu_analysis(){
  gStyle->SetOptStat(0);
  InitializeKinematics(MASS[0], MASS[1], MASS[0], MASS[1], EBEAM, EX); // beam on target
  tc = new TChain("SimulatedTree");
  //  tc->Add("../Simulation/SimulatedTree.root");
  tc->Add(Form("../Simulation/%s", FILENAME));
  //tc->Add("../Simulation/40Ar_197Au_coulex_thin_target_3mm_beam.root");
  
  fReader.SetTree(tc);
  ReadDetectorGeom();
  DefineXYZ_ASGARD();
  DefineXYZ_X6();
  DefineHistos();
  stoppingpowers(0, ZBEAM, ABEAM, ZTARG, ATARG);
  stoppingpowers(1, ZTARG, ATARG, ZTARG, ATARG);
  //  double eBeam_attenuated = EBEAM-GetELoss(0, EBEAM, 0.5*TARGET_THICKNESS, 0);
  //  cout<<eBeam_attenuated<<endl;
  TTreeReaderValue<TASGARDData> ASGARD = {fReader, "ASGARD"};
  TTreeReaderValue<TSTARKData> STARK = {fReader, "STARK"};
  while(fReader.Next()){           
    // i++;
    // if((i % 100000)==0) cout << i << endl;
    //    cout<<STARK->GetMult()<<endl;
    double ep[2] = {-1, -1};
    double thp[2] = {-1, -1};
    double php[2] = {-1, -1};
    double beta[2] = {-1, -1};
    int np = 0;
    // work out algebra if 2-hit
    // using catkin
    double int_depth = 0.;
    double range = 0.;
    double ebeam_attenuated = 0.;
    int pdet = -1;
    //  cout<<STARK->GetMult()<<endl;
    for(int i = 0; i < STARK->GetMult(); i++){
      double ep_raw = STARK->GetFrE(i); // same as GetBkE(i)
      //cout<<i<<" "<<ep_raw<<endl;
      //      cout<<i<<" "<<ep_raw<<endl;
      if (ep_raw > 5.){ // minimum 5 MeV
	int det = STARK->GetDetN(i)-1;
	int strip = STARK->GetFStN(i)-1;
  	TVector3 vSTARK_true = STARK->GetPos(i);
	TVector3 vSTARK_real = vXYZ_X6[strip].at(det);
	// if (ep_raw > 100)
	// cout<<ep_raw<<" "<<vSTARK_true.Z()<<endl;

	double theta;
	double phi;
	if (EXACT_STARKJr_ANGLES){
	  theta = vSTARK_true.Theta();
	  phi = vSTARK_true.Phi();
	}
	else{
	  vSTARK_real.SetZ(vSTARK_true.Z());
	  theta = vSTARK_real.Theta();
	  phi = vSTARK_real.Phi();	  
	}
	ep_z_raw->Fill(vSTARK_true.Z(), ep_raw);


	InitializeKinematics(MASS[0], MASS[1], MASS[0], MASS[1], EBEAM, EX);
  	int pid = BeamOrTarget(EBEAM, ep_raw, theta);
	// allocate Ep, theta and phi according to pid
  	ep[pid] = ep_raw;
	thp[pid] = theta;
	php[pid] = phi;
	
	if (pid==0){
	  // get energy loss in the target by iterative kinematics/SRIM calculations
	  InitializeKinematics(MASS[0], MASS[1], MASS[0], MASS[1], EBEAM, EX);
	  Calculate(thp[pid]*180/TMath::Pi());
	  //	  cout<<vSTARK_true.Z()<<" "<<thp[pid]*180/TMath::Pi()<<" "<<el<<" "<<ep[pid]<<endl;
	  range = GetRange(pid, el, ep[pid]); // get effective depth, in mg/cm2; el needs to be re-calculated
	  if (thp[pid] < TMath::Pi()/2.){ // forward direction
	    int_depth = TARGET_THICKNESS-cos(thp[pid])*range;
	    if (int_depth < 0)
	      int_depth = 0.;
	  }
	  else{ //backward direction
	    int_depth = cos(TMath::Pi()-thp[pid])*range;
	  }

	  for(int iter = 0; iter < NITER; iter++){
	    ebeam_attenuated = EBEAM-GetELoss(0, EBEAM, int_depth, 0);
	    InitializeKinematics(MASS[0], MASS[1], MASS[0], MASS[1], ebeam_attenuated, EX);
	    Calculate(thp[pid]*180/TMath::Pi());
	    range = GetRange(pid, el, ep[pid]); // get effective depth, in mg/cm2; el needs to be re-calculated
	    if (thp[pid] < TMath::Pi()/2.){ // forward direction
	      int_depth = TARGET_THICKNESS-cos(thp[pid])*range;
	      if (int_depth < 0)
		int_depth = 0.;
	    }
	    else{
	      int_depth = cos(TMath::Pi()-thp[pid])*range;
	    }
	    //	    cout<<iter<<" depth: "<<int_depth<<" mg/cm2, range: "<<range<<" mg/cm2"<<endl;
	  }
	  // add back energy lost while traversing target material
	  ep[pid]+=GetELoss(pid, ep[pid], range, 1);
	  
	  ep_z_corrected->Fill(vSTARK_true.Z(), ep[pid]);
	  beta[pid] = GetBeta(ep[pid], MASS[pid]);
	}
	else{ // 197Au
	  InitializeKinematics(MASS[0], MASS[1], MASS[1], MASS[0], EBEAM, EX); // change to target as ejectile
	  Calculate(thp[pid]*180/TMath::Pi());
	  range = GetRange(pid, el, ep[pid]); // get effective depth, in mg/cm2
	  int_depth = TARGET_THICKNESS-range*cos(thp[pid]); // less than actual
	  if (int_depth < 0)
	    int_depth = 0.;

	  for(int iter = 0; iter < NITER; iter++){
	    ebeam_attenuated = EBEAM-GetELoss(0, EBEAM, int_depth, 0);// higher than actual
	    InitializeKinematics(MASS[0], MASS[1], MASS[1], MASS[0], ebeam_attenuated, EX); // change to target as ejectile
	    Calculate(thp[pid]*180/TMath::Pi());
	    range = GetRange(pid, el, ep[pid]); // get effective depth, in mg/cm2; still longer
	    int_depth = TARGET_THICKNESS-range*cos(thp[pid]); // less than actual
	    if (int_depth < 0)
	      int_depth = 0.;
	  }
	  
	  ep[pid]+=GetELoss(pid, ep[pid], range, 1);
	  ep_z_corrected->Fill(vSTARK_true.Z(), ep[pid]);
	  beta[pid] = GetBeta(ep[pid], MASS[pid]);
	}
      }
    }
    // need to reconstruct 2-body kinematics if only 1p detected
    if (ep[0] < 0 && ep[1] > 0){ // only 197Au detected
      InitializeKinematics(MASS[0], MASS[1], MASS[1], MASS[0], ebeam_attenuated, EX); // set target as ejectile
      Calculate(thp[1]*180/TMath::Pi());
      thp[0] = ra*TMath::Pi()/180.; // GetThetaB(ebeam_attenuated, thp[1]); 
      ep[0] = er; //ebeam_attenuated-ep[1];
      beta[0] = GetBeta(ep[0], MASS[0]);
      if (php[1] > 0)
	php[0] = php[1]-TMath::Pi(); // assume symmetry in phi;
      else
	php[0] = php[1]+TMath::Pi();
	
      np = 1;
      pdet = 1;
    }
    else if (ep[1] < 0 && ep[0] > 0){ // only 40Ar detected
      InitializeKinematics(MASS[0], MASS[1], MASS[0], MASS[1], ebeam_attenuated, EX); // set beam as ejectile
      Calculate(thp[0]*180/TMath::Pi());
      thp[1] = ra*TMath::Pi()/180.;// GetThetaT(ebeam_attenuated, thp[0]);
      ep[1] = er; //ebeam_attenuated-ep[0];
      beta[1] = GetBeta(ep[1], MASS[1]);
      if (php[0] > 0)
	php[1] = php[0]-TMath::Pi(); // assume symmetry;
      else
	php[1] = php[0]+TMath::Pi(); // assume symmetry;
      np = 1;
      pdet = 0;
      //      cout<<thp[1]<<" "<<ep[1]<<" "<<beta[1]<<" "<<php[1]<<endl;
    }
    else if (ep[0] > 0){
      np = 2;
    }
    if (np > 0){
      for(int n = 0; n < ASGARD->GetMultiplicityGe(); n++){
  	double eg_raw = ASGARD->GetGeEnergy(n);
	//	cout<<eg_raw<<endl;
	//	hraw_asgard->Fill(eg_raw);
	int clo = ASGARD->GetGeCloverNbr(n);
	int cry = ASGARD->GetGeCrystalNbr(n);
	int seg = ASGARD->GetGeSegmentNbr(n);
  	double thg;
  	double phg;
	if (EXACT_ASGARD_ANGLES){
	  thg = ASGARD->GetGeThetaSemiTrue(n);
	  phg = ASGARD->GetGePhiSemiTrue(n);
	}
	else{ // get from registered angles for each segment
	  TVector3 pos_asgard = vXYZ_ASGARD[cry][seg].at(clo);
	  thg = pos_asgard.Theta();
	  phg = pos_asgard.Phi();
	}
	if (np > 0){
	// if (np==2){
	//   cout<<php[0]<<" "<<php[1]<<endl;
	  //	  cout<<GetDopplerGamma(eg_raw, beta[0], thp[0], php[0], thg, phg)<<endl;
	  meg_thpg_40Ar->Fill(GetCosThetaPG(thp[0], php[0], thg, phg), GetDopplerGamma(eg_raw, beta[0], thp[0], php[0], thg, phg));
	  meg_thpg_197Au->Fill(GetCosThetaPG(thp[1], php[1], thg, phg), GetDopplerGamma(eg_raw, beta[1], thp[1], php[1], thg, phg));
	  // meg_thpg_40Ar->Fill(GetCosThetaPG(thp[0], php[0], thg, phg), eg_raw);
	  // meg_thpg_197Au->Fill(GetCosThetaPG(thp[1], php[1], thg, phg), eg_raw);
	  // if (GetDopplerGamma(eg_raw, beta[1], thp[1], php[1], thg, phg) > 534 && GetDopplerGamma(eg_raw, beta[1], thp[1], php[1], thg, phg) < 548n)
	  //cout<<"eg = "<<eg_raw<<"->"<<GetDopplerGamma(eg_raw, beta[1], thp[1], php[1], thg, phg)<<", cos(pg) = "<<GetCosThetaPG(thp[1], php[1], thg, phg)<<": "<<beta[1]<<" | "<<thp[1]<<" "<<php[1]<<" | "<<thg<<" "<<phg<<" | "<<thp[0]<<" "<<php[0]<<endl;
	  hdop_asgard_40Ar->Fill(GetDopplerGamma(eg_raw, beta[0], thp[0], php[0], thg, phg));
	  hdop_asgard_197Au->Fill(GetDopplerGamma(eg_raw, beta[1], thp[1], php[1], thg, phg));
	}
      }
    }
  }
  //hbeta->Draw();
  TCanvas* ckin = new TCanvas("ckin", "ckin", 1600, 800);
  ckin->Divide(2,1);
  ckin->cd(1);
  ep_z_raw->Draw("colz");
  ckin->cd(2);
  ep_z_corrected->Draw("colz");
  
  ckin->SaveAs("ckin_40Ar_197Au_pencil.png");


  TCanvas* cdop1D = new TCanvas("cdop1D", "cdop1D", 1200, 800);
  hdop_asgard_197Au->SetLineColor(2);
  hdop_asgard_197Au->Draw();
  hdop_asgard_40Ar->SetLineColor(1);
  hdop_asgard_40Ar->Draw("same");
  hraw_asgard->SetLineColor(4);
  //  hraw_asgard->Draw("same");
  
  TLatex* latex = new TLatex();
  latex->SetTextAlign(12);
  latex->SetTextSize(0.04);
  latex->SetTextAngle(90);
  
  latex->DrawLatex(69, 250, "^{197}Au X-ray");
  latex->DrawLatex(275, 240, "^{197}Au 269/279 keV");
  latex->DrawLatex(EG_TARG_MAIN+35, 480, Form("^{197}Au %.0f keV", EG_TARG_MAIN));
  latex->DrawLatex(EG_BEAM_MAIN, 190, Form("^{40}Ar %.0f keV", EG_BEAM_MAIN));
  
  TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend->AddEntry(hdop_asgard_40Ar, "DC for ^{40}Ar", "l");
  legend->AddEntry(hdop_asgard_197Au, "DC for ^{197}Au", "l");
  //  legend->AddEntry(hraw_asgard, "Uncorrected", "l");
  legend->Draw();
  TF1* f40Ar = new TF1("f40Ar", "gaus(0)+pol1(3)", EG_BEAM_MAIN-80, EG_BEAM_MAIN+80);
  f40Ar->SetParameter(1, EG_BEAM_MAIN);
  f40Ar->SetParameter(2, 10);
  hdop_asgard_40Ar->Fit(f40Ar, "QRN0");
  f40Ar->SetLineColor(1);
  f40Ar->Draw("same");
  latex->SetTextAngle(0);
  double fwhm = fabs(f40Ar->GetParameter(2))*2.355;
  latex->DrawLatex(700, 800, Form("FWHM(%.0f) = %.2f keV (%.2f%%)", EG_BEAM_MAIN, fwhm, fwhm/EG_BEAM_MAIN*100.));
  cdop1D->SaveAs("cdop1D_pencil.png");
  
  TCanvas* ceg_cos2D = new TCanvas("ceg_cos2D", "ceg_cos2D", 1600, 800);
  ceg_cos2D->Divide(2,1);
  ceg_cos2D->cd(1);
  meg_thpg_40Ar->Draw("colz");
  ceg_cos2D->cd(2);
  meg_thpg_197Au->Draw("colz");
  //  ceg_cos2D->SaveAs("ceg_cos2D.pdf");
  ceg_cos2D->SaveAs("ceg_cos2D_pencil.png");
}
