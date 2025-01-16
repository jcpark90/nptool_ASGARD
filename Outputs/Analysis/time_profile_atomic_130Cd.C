TChain* tfatima;
TChain* tkhala;

const int ENERGIES[24] = {100, 200, 300, 400, 500,
			  600, 700, 800, 900, 1000,
			  1100, 1200, 1300, 1400, 1500,
			  1600, 1700, 1800, 1900, 2000,
			  2500, 3000, 3500, 4000};
const int dENERGIES[24] = {20, 20, 25, 30, 35,
			   40, 45, 50, 50, 50,
			   50, 50, 50, 50, 50,
			   50, 50, 50, 50, 50,
			   50, 50, 50, 50};

const int REBIN[24] = {1, 1, 2, 2, 2,
		       2, 2, 2, 2, 5,
		       5, 5, 5, 5, 5,
		       5, 5, 5, 5, 5,
		       10, 10, 10, 10};

TH2D* hfatima;
TH2D* hkhala;
TH2D* hidaten;

void time_profile_atomic_130Cd(){
  //  gStyle->SetOptStat(11111111);
  gStyle->SetOptStat(0);
  TCanvas* cc = new TCanvas("cc", "cc", 1200, 900);
  tfatima = new TChain("SimulatedTree");
  tkhala = new TChain("SimulatedTree");
  
    tfatima->Reset();
    tfatima->Add("idaten_130Cd.root");
    tkhala->Reset();
    tkhala->Add("idaten_130Cd.root");
    
    hfatima = new TH2D("hfatima", "", 5000, 0, 500, 300, 0, 1500);
    hkhala = new TH2D("hkhala", "", 5000, 0, 500, 300, 0, 1500);
    
    tfatima->Project("hfatima", "Fatima.fFATIMA_LaBr3_E_Energy*1000:Fatima.fFATIMA_LaBr3_T_Time");
    tkhala->Project("hkhala", "Khala.fKHALA_LaBr3_E_Energy*1000:Khala.fKHALA_LaBr3_T_Time");
    hidaten = (TH2D*)hfatima->Clone("hidaten");
    hidaten->Add(hkhala);
    hidaten->SetTitle("^{130}Cd isomer decay and atomic background");

    TLatex* latex = new TLatex();
    latex->SetTextSize(0.04);
      hidaten->Draw("colz");
      hidaten->GetXaxis()->SetTitle("Time (arb.)");
      hidaten->GetYaxis()->SetTitle("Energy (keV)");
      hidaten->GetXaxis()->CenterTitle();
      hidaten->GetYaxis()->CenterTitle();
      latex->DrawLatex(0.53, 600, "KHALA");
      latex->DrawLatex(0.59, 600, "FATIMA");
      latex->DrawLatex(0.53, 500, "(early)");
      latex->DrawLatex(0.59, 500, "(late)");
      
  //  cs->SaveAs("idaten_138La_bgd_response.png");
  //    cc->SaveAs("idaten_130Cd_atomic_bgd_time_energy.png");
      
}
