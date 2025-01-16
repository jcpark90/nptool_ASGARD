TChain* tfatima;
TChain* tkhala;
TChain* ttigress;


TChain* teu;
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include "TFatimaData.h"
#include "TTigressData.h"

TH1D* hdt_idaten;
TTreeReader fReader;
void GetGammaProjectionLaBr3WithPlastic(TTree* tr, TH1D* hout, double egate, double ewind){
  fReader.SetTree(tr);
  fReader.Restart();
  TTreeReaderValue<TFatimaData> Fatima = {fReader, "Fatima"};
  TTreeReaderValue<TKhalaData> Khala = {fReader, "Khala"};
  TTreeReaderValue<TTigressData> Tigress = {fReader, "Tigress"};
  TTreeReaderValue<TPlasticData> Plastic = {fReader, "Plastic"};

  long long i = 0;

  while(fReader.Next()){           
    i++;
    bool coinc = 0;
    int index_fatima = -1;
    int index_khala = -1;
    
    if((i % 100000)==0) cout << i << endl;
    if (Plastic->GetMult() > 0){
      for(int n = 0; n < Fatima->GetFatimaLaBr3EMult();n++){
	if (fabs(Fatima->GetFatimaLaBr3EEnergy(n)-egate) < ewind){ // start gamma in index m
	  coinc = 1;
	  index_fatima = n;
	  break;
	}
      }
      if (!coinc){
	for(int n = 0; n < Khala->GetKhalaLaBr3EMult();n++){
	  if (fabs(Khala->GetKhalaLaBr3EEnergy(n)-egate) < ewind){ // start gamma in index n
	    coinc = 1;
	    index_khala = n;
	    break;
	  }
	}
      }
      else{
	for(int n = 0; n < Fatima->GetFatimaLaBr3EMult();n++){
	  if (n==index_fatima) continue;
	  else
	    hout->Fill(Fatima->GetFatimaLaBr3EEnergy(n));
	}
	for(int n = 0; n < Khala->GetKhalaLaBr3EMult();n++){
	  if (n==index_khala) continue;
	  else
	    hout->Fill(Khala->GetKhalaLaBr3EEnergy(n));
	}
	
      }
    }
    
  }
  return hout;
}
void GetGammaProjectionLaBr3(TTree* tr, TH1D* hout, double egate, double ewind){
  fReader.SetTree(tr);
  fReader.Restart();
  TTreeReaderValue<TFatimaData> Fatima = {fReader, "Fatima"};
  TTreeReaderValue<TKhalaData> Khala = {fReader, "Khala"};
  TTreeReaderValue<TTigressData> Tigress = {fReader, "Tigress"};
  long long i = 0;

  while(fReader.Next()){           
    i++;
    bool coinc = 0;
    int index_fatima = -1;
    int index_khala = -1;
    
    if((i % 100000)==0) cout << i << endl;
    for(int n = 0; n < Fatima->GetFatimaLaBr3EMult();n++){
      if (fabs(Fatima->GetFatimaLaBr3EEnergy(n)-egate) < ewind){ // start gamma in index m
	coinc = 1;
	index_fatima = n;
	break;
      }
    }
    if (!coinc){
      for(int n = 0; n < Khala->GetKhalaLaBr3EMult();n++){
	if (fabs(Khala->GetKhalaLaBr3EEnergy(n)-egate) < ewind){ // start gamma in index n
	  coinc = 1;
	  index_khala = n;
	  break;
	}
      }
    }
    else{
      for(int n = 0; n < Fatima->GetFatimaLaBr3EMult();n++){
	if (n==index_fatima) continue;
	else
	  hout->Fill(Fatima->GetFatimaLaBr3EEnergy(n));
      }
      for(int n = 0; n < Khala->GetKhalaLaBr3EMult();n++){
	if (n==index_khala) continue;
	else
	  hout->Fill(Khala->GetKhalaLaBr3EEnergy(n));
      }
      
    }
    
  }
  return hout;
}
void GetGammaProjectionHPGe(TTree* tr, TH1D* hout, double egate, double ewind){
  fReader.SetTree(tr);
  fReader.Restart();
  TTreeReaderValue<TFatimaData> Fatima = {fReader, "Fatima"};
  TTreeReaderValue<TKhalaData> Khala = {fReader, "Khala"};
  TTreeReaderValue<TTigressData> Tigress = {fReader, "Tigress"};
  long long i = 0;

  while(fReader.Next()){           
    i++;
    bool coinc = 0;
    if((i % 100000)==0) cout << i << endl;
    for(int n = 0; n < Tigress->GetMultiplicityGe();n++){
      //      cout<<Tigress->GetGeEnergy(n)<<endl;
      if (fabs(Tigress->GetGeEnergy(n)-egate)<ewind){
	coinc = 1;
	break;
      }
    }
   
    if(coinc){
      for(int n = 0; n < Fatima->GetFatimaLaBr3EMult();n++){
	hout->Fill(Fatima->GetFatimaLaBr3EEnergy(n));
      }
      for(int n = 0; n < Khala->GetKhalaLaBr3EMult();n++){
	hout->Fill(Khala->GetKhalaLaBr3EEnergy(n));
      }
      
    }
    
  }
  return hout;
}

void GetStartStopHist(TTree* tr, TH1D* hout, double estart, double ewind_start, double estop, double ewind_stop){
  fReader.SetTree(tr);
  fReader.Restart();
  TTreeReaderValue<TFatimaData> Fatima = {fReader, "Fatima"};
  TTreeReaderValue<TKhalaData> Khala = {fReader, "Khala"};
  TTreeReaderValue<TTigressData> Tigress = {fReader, "Tigress"};
 
  long long i = 0;
  
  while(fReader.Next()){           
    i++;
    if((i % 100000)==0) cout << i << endl;
    for(int n = 0; n < Fatima->GetFatimaLaBr3EMult();n++){      
      for(int m = n+1; m < Fatima->GetFatimaLaBr3EMult(); m++){  
      	if (fabs(Fatima->GetFatimaLaBr3EEnergy(m)-estart) < ewind_start // start gamma in index m
	    && fabs(Fatima->GetFatimaLaBr3EEnergy(n)-estop) < ewind_stop) // stop gamma in index n
      	  hout->Fill(Fatima->GetFatimaLaBr3TTime(n)-Fatima->GetFatimaLaBr3TTime(m));
      	else if (fabs(Fatima->GetFatimaLaBr3EEnergy(n)-estart) < ewind_start   // start gamma in index n
		 && fabs(Fatima->GetFatimaLaBr3EEnergy(m)-estop) < ewind_stop) // stop gamma in index m 
      	  hout->Fill(Fatima->GetFatimaLaBr3TTime(m)-Fatima->GetFatimaLaBr3TTime(n));
      }
    }
    for(int n = 0; n < Khala->GetKhalaLaBr3EMult();n++){
      for(int m = n+1; m < Khala->GetKhalaLaBr3EMult(); m++){  
       	if (fabs(Khala->GetKhalaLaBr3EEnergy(m)-estart) < ewind_start // start gamma in index m
	    && fabs(Khala->GetKhalaLaBr3EEnergy(n)-estop) < ewind_stop) // stop gamma in index n
      	  hout->Fill(Khala->GetKhalaLaBr3TTime(n)-Khala->GetKhalaLaBr3TTime(m));
      	else if (fabs(Khala->GetKhalaLaBr3EEnergy(n)-estart) < ewind_start   // start gamma in index n
		 && fabs(Khala->GetKhalaLaBr3EEnergy(m)-estop) < ewind_stop) // stop gamma in index m
      	  hout->Fill(Khala->GetKhalaLaBr3TTime(m)-Khala->GetKhalaLaBr3TTime(n));
      	
      }
    }
    //    loop through coincidence data for mixed detector types
    for(int n = 0; n < Fatima->GetFatimaLaBr3EMult();n++){
      for(int m = 0; m < Khala->GetKhalaLaBr3EMult(); m++){
     	if (fabs(Khala->GetKhalaLaBr3EEnergy(m)-estart) < ewind_start // start gamma in index m for khala
	    && fabs(Fatima->GetFatimaLaBr3EEnergy(n)-estop) < ewind_stop) // stop gamma in index n for fatima
      	  hout->Fill(Fatima->GetFatimaLaBr3TTime(n)-Khala->GetKhalaLaBr3TTime(m));
      	else if (fabs(Fatima->GetFatimaLaBr3EEnergy(n)-estart) < ewind_start   // start gamma in index n for fatima
		 && fabs(Khala->GetKhalaLaBr3EEnergy(m)-estop) < ewind_stop) // stop gamma in index m for khala
      	  hout->Fill(Khala->GetKhalaLaBr3TTime(m)-Fatima->GetFatimaLaBr3TTime(n));
      }
    }

  }
  return hout;
}
void GetStartStopHistWithGamma(TTree* tr, TH1D* hout, double estart, double ewind_start, double estop, double ewind_stop, double egam, double ewind_gam){
  fReader.SetTree(tr);
  fReader.Restart();
  TTreeReaderValue<TFatimaData> Fatima = {fReader, "Fatima"};
  TTreeReaderValue<TKhalaData> Khala = {fReader, "Khala"};
  TTreeReaderValue<TTigressData> Tigress = {fReader, "Tigress"};
 
  long long i = 0;
  while(fReader.Next()){           
    i++;
    bool gamma_gated = 0;
    //   cout<<Tigress->GetMultiplicityGe()<<endl;
    for(int n = 0; n < Tigress->GetMultiplicityGe();n++){
      //      cout<<Tigress->GetGeEnergy(n)<<endl;
      if (fabs(Tigress->GetGeEnergy(n)-egam)<ewind_gam){
	gamma_gated = 1;
	break;
      }
    }
    if (gamma_gated){
      for(int n = 0; n < Fatima->GetFatimaLaBr3EMult();n++){      
	for(int m = n+1; m < Fatima->GetFatimaLaBr3EMult(); m++){  
	  if (fabs(Fatima->GetFatimaLaBr3EEnergy(m)-estart) < ewind_start // start gamma in index m
	      && fabs(Fatima->GetFatimaLaBr3EEnergy(n)-estop) < ewind_stop) // stop gamma in index n
	    hout->Fill(Fatima->GetFatimaLaBr3TTime(n)-Fatima->GetFatimaLaBr3TTime(m));
	  else if (fabs(Fatima->GetFatimaLaBr3EEnergy(n)-estart) < ewind_start   // start gamma in index n
		   && fabs(Fatima->GetFatimaLaBr3EEnergy(m)-estop) < ewind_stop) // stop gamma in index m 
	    hout->Fill(Fatima->GetFatimaLaBr3TTime(m)-Fatima->GetFatimaLaBr3TTime(n));
	}
      }
      for(int n = 0; n < Khala->GetKhalaLaBr3EMult();n++){
	for(int m = n+1; m < Khala->GetKhalaLaBr3EMult(); m++){  
	  if (fabs(Khala->GetKhalaLaBr3EEnergy(m)-estart) < ewind_start // start gamma in index m
	      && fabs(Khala->GetKhalaLaBr3EEnergy(n)-estop) < ewind_stop) // stop gamma in index n
	    hout->Fill(Khala->GetKhalaLaBr3TTime(n)-Khala->GetKhalaLaBr3TTime(m));
	  else if (fabs(Khala->GetKhalaLaBr3EEnergy(n)-estart) < ewind_start   // start gamma in index n
		   && fabs(Khala->GetKhalaLaBr3EEnergy(m)-estop) < ewind_stop) // stop gamma in index m
	    hout->Fill(Khala->GetKhalaLaBr3TTime(m)-Khala->GetKhalaLaBr3TTime(n));
	  
	}
      }
      //    loop through coincidence data for mixed detector types
      for(int n = 0; n < Fatima->GetFatimaLaBr3EMult();n++){
	for(int m = 0; m < Khala->GetKhalaLaBr3EMult(); m++){
	  if (fabs(Khala->GetKhalaLaBr3EEnergy(m)-estart) < ewind_start // start gamma in index m for khala
	      && fabs(Fatima->GetFatimaLaBr3EEnergy(n)-estop) < ewind_stop) // stop gamma in index n for fatima
	    hout->Fill(Fatima->GetFatimaLaBr3TTime(n)-Khala->GetKhalaLaBr3TTime(m));
	  else if (fabs(Fatima->GetFatimaLaBr3EEnergy(n)-estart) < ewind_start   // start gamma in index n for fatima
		   && fabs(Khala->GetKhalaLaBr3EEnergy(m)-estop) < ewind_stop) // stop gamma in index m for khala
	    hout->Fill(Khala->GetKhalaLaBr3TTime(m)-Fatima->GetFatimaLaBr3TTime(n));
	}
      }

    }
  }
  return hout;
}
void GetStartLaBr3StopGamma(TTree* tr, TH1D* hout, double estart, double ewind_start, double estop, double ewind_stop){
  fReader.SetTree(tr);
  fReader.Restart();
  TTreeReaderValue<TFatimaData> Fatima = {fReader, "Fatima"};
  TTreeReaderValue<TKhalaData> Khala = {fReader, "Khala"};
  TTreeReaderValue<TTigressData> Tigress = {fReader, "Tigress"};
 
  long long i = 0;
  bool gamma_gated = 0;
  while(fReader.Next()){           
    i++;
    for(int n = 0; n < Fatima->GetFatimaLaBr3EMult(); n++){      
      for(int m = 0; m < Tigress->GetMultiplicityGe(); m++){  
	if (fabs(Fatima->GetFatimaLaBr3EEnergy(n)-estart) < ewind_start // start gamma in index m
	    && fabs(Tigress->GetGeEnergy(m)-estop) < ewind_stop) // stop gamma in index n
	  hout->Fill(Tigress->GetGeTimeLED(m)-Fatima->GetFatimaLaBr3TTime(n));
      }
    }
    for(int n = 0; n < Khala->GetKhalaLaBr3EMult(); n++){      
      for(int m = 0; m < Tigress->GetMultiplicityGe(); m++){  
	if (fabs(Khala->GetKhalaLaBr3EEnergy(n)-estart) < ewind_start // start gamma in index m
	    && fabs(Tigress->GetGeEnergy(m)-estop) < ewind_stop) // stop gamma in index n
	  hout->Fill(Tigress->GetGeTimeLED(m)-Khala->GetKhalaLaBr3TTime(n));
      }
    }
    
  }
  return hout;
}
void GetStartStopHistWithPlastic(TTree* tr, TH1D* hout, double estop, double ewind_stop){
  fReader.SetTree(tr);
  fReader.Restart();
  TTreeReaderValue<TFatimaData> Fatima = {fReader, "Fatima"};
  TTreeReaderValue<TKhalaData> Khala = {fReader, "Khala"};
  TTreeReaderValue<TTigressData> Tigress = {fReader, "Tigress"};
  TTreeReaderValue<TWAS3ABiData> Was3abi = {fReader, "WAS3ABi"};
  TTreeReaderValue<TPlasticData> Plastic = {fReader, "Plastic"};

  int i = 0;
  
  while(fReader.Next()){           
    i++;
    if((i % 100000)==0) cout << i << endl;
    double tpl_temp = 1.e12;
    for(int n = 0; n < Plastic->GetMult();n++){
      if (Plastic->GetTime(n) < tpl_temp)
	tpl_temp = Plastic->GetTime(n); // get the fastest time hit on the plastic detector
    }
    if (tpl_temp < 1.e12){ // check that we have a plastic hit in this event
      for(int n = 0; n < Fatima->GetFatimaLaBr3EMult();n++){
	if (fabs(Fatima->GetFatimaLaBr3EEnergy(n)-estop) < ewind_stop) 
	  hout->Fill(Fatima->GetFatimaLaBr3TTime(n)-tpl_temp);
      }
      for(int n = 0; n < Khala->GetKhalaLaBr3EMult();n++){
      	if (fabs(Khala->GetKhalaLaBr3EEnergy(n)-estop) < ewind_stop) 
	  hout->Fill(Khala->GetKhalaLaBr3TTime(n)-tpl_temp);
      }
    }    
  }
  //  return hout;
}
void GetLaBr3EnergySpectrum(TTree* tr, TH1D* hlabr3){
  TH1D* hfatima = (TH1D*)hlabr3->Clone("hfatima");
  TH1D* hkhala = (TH1D*)hlabr3->Clone("hkhala");
  tr->Project("hfatima", "Fatima.fFATIMA_LaBr3_E_Energy");
  tr->Project("hkhala", "Khala.fKHALA_LaBr3_E_Energy");
  hlabr3->Add(hfatima);
  hlabr3->Add(hkhala);
}
void GetLaBr3EnergySpectrumWithPlastic(TTree* tr, TH1D* hlabr3){
  TH1D* hfatima = (TH1D*)hlabr3->Clone("hfatima");
  TH1D* hkhala = (TH1D*)hlabr3->Clone("hkhala");
  tr->Project("hfatima", "Fatima.fFATIMA_LaBr3_E_Energy", "Plastic.fPlastic_Time > 0");
  tr->Project("hkhala", "Khala.fKHALA_LaBr3_E_Energy", "Plastic.fPlastic_Time > 0");
  hlabr3->Add(hfatima);
  hlabr3->Add(hkhala);
}

void GetHPGeEnergySpectrum(TTree* tr, TH1D* hHPGe){
  tr->Project(hHPGe->GetName(), "Tigress.fTIG_Ge_Energy");
}
void GetHPGeEnergySpectrumWithPlastic(TTree* tr, TH1D* hHPGe){
  tr->Project(hHPGe->GetName(), "Tigress.fTIG_Ge_Energy", "Plastic.fPlastic_Time > 0");
}
void start_stop_analysis_watanabe(){
  // gStyle->SetOptStat(0);

  teu = new TChain("SimulatedTree");

  // 130Cd
  // teu->Add("130Cd_isomer.root");
  // TH1D* he_idaten = new TH1D("he_idaten", "", 2000, 0, 2000);
  // GetLaBr3EnergySpectrum(teu, he_idaten);
  // TH1D* he_hpge = new TH1D("he_hpge", "", 2000, 0, 2000);
  // GetHPGeEnergySpectrum(teu, he_hpge);  
  // TH1D* he_idaten_hpge_gate = new TH1D("he_idaten_hpge_gate", "", 2000, 0, 2000);
  // GetGammaProjectionHPGe(teu, he_idaten_hpge_gate, 128, 2);
  // TH1D* ht_idaten_128_hpge_138 = new TH1D("ht_idaten_138_hpge_138", "", 1000, -50, 950);  
  // GetStartLaBr3StopGamma(teu, ht_idaten_128_hpge_138, 128, 5, 138, 2);
  // TH1D* ht_idaten_138_539_hpge_128 = new TH1D("ht_idaten_138_539_hpge_128", "", 2000, -10, 10);
  // GetStartStopHistWithGamma(teu, ht_idaten_138_539_hpge_128, 138, 10, 539, 25, 128, 2);
  // TH1D* ht_idaten_138_1325_hpge_128 = new TH1D("ht_idaten_138_1325_hpge_128", "", 2000, -10, 10);
  // GetStartStopHistWithGamma(teu, ht_idaten_138_1325_hpge_128, 138, 10, 1325, 40, 128, 2);
  // TFile* f130Cd = new TFile("outputs_130Cd.root", "recreate");
  // he_idaten->Write();
  // he_hpge->Write();
  // he_idaten_hpge_gate->Write();
  // ht_idaten_128_hpge_138->Write();
  // ht_idaten_138_539_hpge_128->Write();
  // ht_idaten_138_1325_hpge_128->Write();
  // f130Cd->Close();

  // 126Pd
  // teu->Add("126Pd_isomer.root");
  // TH1D* he_idaten = new TH1D("he_idaten", "", 2000, 0, 2000);
  // GetLaBr3EnergySpectrum(teu, he_idaten);
  // TH1D* he_hpge = new TH1D("he_hpge", "", 2000, 0, 2000);
  // GetHPGeEnergySpectrum(teu, he_hpge);  
  // TH1D* he_idaten_labr3_gate = new TH1D("he_idaten_labr3_gate", "", 2000, 0, 2000);
  // GetGammaProjectionLaBr3(teu, he_idaten_labr3_gate, 788, 30);
  // TH1D* ht_idaten_788_693 = new TH1D("ht_idaten_788_693", "", 2000, -10, 10);  
  // GetStartStopHist(teu, ht_idaten_788_693, 788, 30, 693, 30);
  // TH1D* ht_idaten_1330_693 = new TH1D("ht_idaten_1330_693", "", 2000, -10, 10);
  // GetStartStopHist(teu, ht_idaten_1330_693, 1330, 40, 693, 30);
  // TH1D* ht_idaten_542_788 = new TH1D("ht_idaten_542_788", "", 2000, -10, 10);
  // GetStartStopHist(teu, ht_idaten_542_788, 542, 25, 788, 30);
  // TFile* f126Pd = new TFile("outputs_126Pd.root", "recreate");
  // he_idaten->Write();
  // he_hpge->Write();
  // he_idaten_labr3_gate->Write();
  // ht_idaten_788_693->Write();
  // ht_idaten_1330_693->Write();
  // ht_idaten_542_788->Write();
  // f126Pd->Close();
  
  // 128Pd
  // teu->Add("128Pd_isomer.root");
  // TH1D* he_idaten = new TH1D("he_idaten", "", 2000, 0, 2000);
  // GetLaBr3EnergySpectrum(teu, he_idaten);
  // TH1D* he_hpge = new TH1D("he_hpge", "", 2000, 0, 2000);
  // GetHPGeEnergySpectrum(teu, he_hpge);  
  // TH1D* he_idaten_labr3_gate = new TH1D("he_idaten_labr3_gate", "", 2000, 0, 2000);
  // GetGammaProjectionLaBr3(teu, he_idaten_labr3_gate, 504, 25);
  // TH1D* ht_idaten_75_260 = new TH1D("ht_idaten_75_260", "", 200, -10, 190);  
  // GetStartStopHist(teu, ht_idaten_75_260, 75, 10, 260, 20);
  // TH1D* ht_idaten_75_504 = new TH1D("ht_idaten_75_504", "", 200, -10, 190);  
  // GetStartStopHist(teu, ht_idaten_75_504, 75, 10, 504, 25);
  // TH1D* ht_idaten_75_1311 = new TH1D("ht_idaten_75_1311", "", 200, -10, 190);  
  // GetStartStopHist(teu, ht_idaten_75_1311, 75, 10, 1311, 40);
  // TH1D* ht_idaten_260_504 = new TH1D("ht_idaten_260_504", "", 200, -10, 10);  
  // GetStartStopHist(teu, ht_idaten_260_504, 260, 20, 504, 25);
  // TH1D* ht_idaten_260_1311 = new TH1D("ht_idaten_260_1311", "", 200, -10, 10);  
  // GetStartStopHist(teu, ht_idaten_260_1311, 260, 20, 1311, 40);
  // TFile* f128Pd = new TFile("outputs_128Pd.root", "recreate");
  // he_idaten->Write();
  // he_hpge->Write();
  // he_idaten_labr3_gate->Write();
  // ht_idaten_75_260->Write();
  // ht_idaten_75_504->Write();
  // ht_idaten_75_1311->Write();
  // ht_idaten_260_504->Write();
  // ht_idaten_260_1311->Write();
  // f128Pd->Close();
  
  // 128Cd
  // teu->Add("128Cd_isomer_0.root");
  // teu->Add("128Cd_isomer_1.root");
  // teu->Add("128Cd_isomer_2.root");
  // teu->Add("128Cd_isomer_3.root");
  // teu->Add("128Cd_isomer_4.root");
  // TH1D* he_idaten = new TH1D("he_idaten", "", 2000, 0, 2000);
  // GetLaBr3EnergySpectrum(teu, he_idaten);
  // TH1D* he_hpge = new TH1D("he_hpge", "", 2000, 0, 2000);
  // GetHPGeEnergySpectrum(teu, he_hpge);  
  // TH1D* he_idaten_labr3_gate = new TH1D("he_idaten_labr3_gate", "", 2000, 0, 2000);
  // GetGammaProjectionLaBr3(teu, he_idaten_labr3_gate, 785, 30);
  // TH1D* ht_idaten_785_646 = new TH1D("ht_idaten_785_646", "", 2000, -10, 10);  
  // GetStartStopHist(teu, ht_idaten_785_646, 785, 30, 646, 30);
  // TH1D* ht_idaten_440_785 = new TH1D("ht_idaten_440_785", "", 2000, -10, 10);  
  // GetStartStopHist(teu, ht_idaten_440_785, 440, 25, 785, 30);
  // TFile* f128Cd = new TFile("outputs_128Cd.root", "recreate");
  // he_idaten->Write();
  // he_hpge->Write();
  // he_idaten_labr3_gate->Write();
  // ht_idaten_785_646->Write();
  // ht_idaten_440_785->Write();
  // f128Cd->Close();
  
  // 150Cs
  // teu->Add("150Cs_decay_reduced.root");
  // TH1D* he_idaten = new TH1D("he_idaten", "", 2000, 0, 2000);
  // GetLaBr3EnergySpectrumWithPlastic(teu, he_idaten);
  // TH1D* he_hpge = new TH1D("he_hpge", "", 2000, 0, 2000);
  // GetHPGeEnergySpectrum(teu, he_hpge);  
  // TH1D* he_idaten_labr3_gate = new TH1D("he_idaten_labr3_gate", "", 2000, 0, 2000);
  // GetGammaProjectionLaBr3WithPlastic(teu, he_idaten_labr3_gate, 101, 12);
  // TH1D* ht_idaten_plastic_101 = new TH1D("ht_idaten_plastic_101", "", 500, -10, 40);
  // GetStartStopHistWithPlastic(teu, ht_idaten_plastic_101, 101, 12);
  // TH1D* ht_idaten_plastic_512 = new TH1D("ht_idaten_plastic_512", "", 500, -10, 40);
  // GetStartStopHistWithPlastic(teu, ht_idaten_plastic_512, 512, 25);
  // TH1D* ht_idaten_plastic_380 = new TH1D("ht_idaten_plastic_380", "", 500, -10, 40);
  // GetStartStopHistWithPlastic(teu, ht_idaten_plastic_380, 380, 25);
  // TH1D* ht_idaten_plastic_613 = new TH1D("ht_idaten_plastic_613", "", 500, -10, 40);
  // GetStartStopHistWithPlastic(teu, ht_idaten_plastic_613, 613, 10);
  // TH1D* ht_idaten_plastic_597 = new TH1D("ht_idaten_plastic_597", "", 500, -10, 40);
  // GetStartStopHistWithPlastic(teu, ht_idaten_plastic_597, 597, 10);
  // TFile* f150Cs = new TFile("outputs_150Cs_reduced.root", "recreate");
  // he_idaten->Write();
  // he_hpge->Write();
  // he_idaten_labr3_gate->Write();
  // ht_idaten_plastic_101->Write();
  // ht_idaten_plastic_512->Write();
  // ht_idaten_plastic_380->Write();
  // ht_idaten_plastic_613->Write();
  // ht_idaten_plastic_597->Write();
  // f150Cs->Close();

  // 148Cs
  // teu->Add("148Cs_decay_0.root");
  // teu->Add("148Cs_decay_1.root");
  // teu->Add("148Cs_decay_2.root");
  // teu->Add("148Cs_decay_3.root");
  // teu->Add("148Cs_decay_4.root");
  teu->Add("148Cs_decay_reduced.root");
  TH1D* he_idaten = new TH1D("he_idaten", "", 2000, 0, 2000);
  GetLaBr3EnergySpectrumWithPlastic(teu, he_idaten);
  TH1D* he_hpge = new TH1D("he_hpge", "", 2000, 0, 2000);
  GetHPGeEnergySpectrum(teu, he_hpge);  
  TH1D* he_idaten_labr3_gate = new TH1D("he_idaten_labr3_gate", "", 2000, 0, 2000);
  GetGammaProjectionLaBr3WithPlastic(teu, he_idaten_labr3_gate, 142, 10);
  TH1D* ht_idaten_plastic_142 = new TH1D("ht_idaten_plastic_142", "", 2500, -5, 20);
  GetStartStopHistWithPlastic(teu, ht_idaten_plastic_142, 142, 10);
  TH1D* ht_idaten_plastic_687 = new TH1D("ht_idaten_plastic_687", "", 2500, -5, 20);
  GetStartStopHistWithPlastic(teu, ht_idaten_plastic_687, 687, 30);
  TH1D* ht_idaten_plastic_546 = new TH1D("ht_idaten_plastic_546", "", 2500, -5, 20);
  GetStartStopHistWithPlastic(teu, ht_idaten_plastic_546, 546, 25);
  TH1D* ht_idaten_plastic_633 = new TH1D("ht_idaten_plastic_633", "", 2500, -5, 20);
  GetStartStopHistWithPlastic(teu, ht_idaten_plastic_633, 633, 30);

  //  TFile* f148Cs = new TFile("outputs_148Cs.root", "recreate");
  TFile* f148Cs = new TFile("outputs_148Cs_reduced.root", "recreate");
  he_idaten->Write();
  he_hpge->Write();
  he_idaten_labr3_gate->Write();
  ht_idaten_plastic_142->Write();
  ht_idaten_plastic_687->Write();
  ht_idaten_plastic_546->Write();
  ht_idaten_plastic_633->Write();
  f148Cs->Close();
}
