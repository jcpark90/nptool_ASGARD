// masses in MeV
const double UMASS = 931.494;
double mb; //beam
double mt; //target 
double me; //ejectile 
double mr; //recoil

double eb; // beam energy in MeV

double qf; // final state reaction Q-value, cell $O$10

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

vector<pair<double, double> > cd_theta_beam;
vector<pair<double, double> > cd_theta_targ;


void Initialize(double mb, double mt, double me, double mr, double eb, double qf){
  y = mb + mt; // total mass of system
  betac = sqrt(eb*(eb+2*mb))/(y + eb);
  ecmi = sqrt(y*y+(2*eb*mt));
  ecmf = ecmi+qf-y+me+mr;
  e3cm = (ecmf*ecmf+(me+mr)*(me-mr))/(2*ecmf);
  beta3c = sqrt(1.-1./pow(e3cm/me, 2.));
  y_new = pow(e3cm/me, 2.)*(1.-betac*betac);
  c = 1.-y_new;
}
void Calculate(double la){ // lab angle in degrees
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
  //cout<<cosagl<<" "<<b<<" "<<a<<" "<<dsq<<" "<<b3L1<<" "<<y_loop<<" "<<y_adj<<" "<<cma<<" "<<el<<" "<<rcml<<" "<<er<<" "<<angl4<<" "<<ra<<endl;
}
double GetE_ejectile(){ return el; }
double GetE_recoil(){ return er; }
double GetAng_recoil(){ return ra; }

void SetCDAngles(double distance){
  for(int i = 0; i < 16; i++){
    cd_theta_beam.push_back(make_pair(atan((9.0 + 2*i)/distance)*180.0/M_PI, atan((10.9 + 2*i)/distance)*180.0/M_PI));
   
    Calculate(cd_theta_beam[i].first);
    double ra_first = ra;
    Calculate(cd_theta_beam[i].second);
    double ra_second = ra;
    cd_theta_targ.push_back(make_pair(ra_first, ra_second));
  }
  //  angle_beam;
}
 
void catkin_root(){
  // mb = 126598;
  // mt = 184407;
  // me = 126598;
  // mr = 184407;

  // eb = 952.;
  // qf = 0.;
  // Initialize(mb, mt, me, mr, eb, qf);

  // mb = 110.*UMASS;
  // mt = 206.*UMASS;
  // me = 110.*UMASS;
  // mr = 206.*UMASS;
  mb = 102379;
  mt = 191866;
  // me = 102379;
  // mr = 191866;
  mr = 102379;
  me = 191866;

  eb = 450.;
  qf = -2.015;
  Initialize(mb, mt, me, mr, eb, qf);


  Calculate(10); // quantities at lab angle 40 degrees
  double ep = GetE_ejectile();
  double er = GetE_recoil();
  cout<<mb<<" "<<mt<<endl;
  cout<<ep<<" "<<er<<" "<<ra<<endl;
  SetCDAngles(22.94);
  std::cout << std::setprecision(3);
  for(int i = 0; i < 16; i++){
    std::cout<<i<<" "<<cd_theta_beam[i].first<<"-"<<cd_theta_beam[i].second<<" | "<<cd_theta_targ[i].first<<"-"<<cd_theta_targ[i].second<<std::endl;
  }
  
}
