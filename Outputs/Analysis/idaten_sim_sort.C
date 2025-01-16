TChain* tfatima;
TChain* tkhala;
TChain* ttigress;


TChain* teu = new TChain("SimulatedTree");
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include "TFatimaData.h"
#include "TTigressData.h"

TH1D* hdt_idaten;
TTreeReader fReader;
#define HAS_FATIMA 0
#define HAS_KHALA 1
// void SetFatimaOn(bool onoff){
//   HAS_FATIMA = onoff;
// }
// void SetKhalaOn(bool onoff){
//   HAS_KHALA = onoff;
// }
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
      if (HAS_FATIMA){
	for(int n = 0; n < Fatima->GetFatimaLaBr3EMult();n++){
	  if (fabs(Fatima->GetFatimaLaBr3EEnergy(n)-egate) < ewind){ // start gamma in index m
	    coinc = 1;
	    index_fatima = n;
	    break;
	  }
	}
      }
      if (!coinc){
	if (HAS_KHALA){
	  for(int n = 0; n < Khala->GetKhalaLaBr3EMult();n++){
	    if (fabs(Khala->GetKhalaLaBr3EEnergy(n)-egate) < ewind){ // start gamma in index n
	      coinc = 1;
	      index_khala = n;
	      break;
	    }
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
void GetGammaGammaMatrix(TTree* tr, TH2D* mout){ // hpge and IDATEN coincidence matrix
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
      for(int m = 0; m < Fatima->GetFatimaLaBr3EMult();m++){
	mout->Fill(Fatima->GetFatimaLaBr3EEnergy(m), Tigress->GetGeEnergy(n));
      }
      for(int m = 0; m < Khala->GetKhalaLaBr3EMult();m++){
	mout->Fill(Khala->GetKhalaLaBr3EEnergy(m), Tigress->GetGeEnergy(n));
      }
      
    }
    
  }
  return mout;
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
  if (HAS_FATIMA)
  TTreeReaderValue<TFatimaData> Fatima = {fReader, "Fatima"};
  TTreeReaderValue<TKhalaData> Khala = {fReader, "Khala"};
  TTreeReaderValue<TTigressData> Tigress = {fReader, "Tigress"};
 
  long long i = 0;
  //  TH1D* hbgd = hout->
  while(fReader.Next()){           
    i++;
    if((i % 100000)==0) cout << i << endl;
    if (HAS_FATIMA){
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
    }
    if (HAS_KHALA){
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
    }
    //    loop through coincidence data for mixed detector types
    if (HAS_FATIMA && HAS_KHALA){
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
    //    double tpl_temp = 0;
    //    double ene_temp = 30;
    for(int n = 0; n < Plastic->GetMult();n++){
      if (Plastic->GetTime(n) < tpl_temp)
	//      if (Plastic->GetTime(n) > tpl_temp)
	tpl_temp = Plastic->GetTime(n); // get the fastest time hit on the plastic detector
    }
    if (tpl_temp < 1.e12){ // check that we have a plastic hit in this event
      //if (tpl_temp > 0){ // check that we have a plastic hit in this event
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
void GetLaBr3EnergySpectrumWAS3ABiGate(TTree* tr, TH1D* hlabr3, double emax_WAS3ABi){
  TH1D* hfatima = (TH1D*)hlabr3->Clone("hfatima");
  TH1D* hkhala = (TH1D*)hlabr3->Clone("hkhala");
  tr->Project("hfatima", "Fatima.fFATIMA_LaBr3_E_Energy", Form("0.5*(WAS3ABi.fWAS3ABi_StripFront_Energy+WAS3ABi.fWAS3ABi_StripFront_Energy) < %.0lf", emax_WAS3ABi));
  tr->Project("hkhala", "Khala.fKHALA_LaBr3_E_Energy", Form("0.5*(WAS3ABi.fWAS3ABi_StripFront_Energy+WAS3ABi.fWAS3ABi_StripFront_Energy) < %.0lf", emax_WAS3ABi));
  hlabr3->Add(hfatima);
  hlabr3->Add(hkhala);
}
void GetLaBr3EnergySpectrumTimeGate(TTree* tr, TH1D* hlabr3, double tmin, double tmax){
  TH1D* hfatima = (TH1D*)hlabr3->Clone("hfatima");
  TH1D* hkhala = (TH1D*)hlabr3->Clone("hkhala");
  tr->Project("hfatima", "Fatima.fFATIMA_LaBr3_E_Energy", Form("Fatima.fFATIMA_LaBr3_T_Time > %e && Fatima.fFATIMA_LaBr3_T_Time < %e", tmin, tmax));
  tr->Project("hkhala", "Khala.fKHALA_LaBr3_E_Energy", Form("Khala.fKHALA_LaBr3_T_Time > %e && Khala.fKHALA_LaBr3_T_Time < %e", tmin, tmax));
  hlabr3->Add(hfatima);
  hlabr3->Add(hkhala);
}
void GetLaBr3EnergySpectrumWithPlastic(TTree* tr, TH1D* hlabr3){
  // TH1D* hfatima = (TH1D*)hlabr3->Clone("hfatima");
  // TH1D* hkhala = (TH1D*)hlabr3->Clone("hkhala");
  // tr->Project("hfatima", "Fatima.fFATIMA_LaBr3_E_Energy", "Plastic.fPlastic_Time > 0");
  // tr->Project("hkhala", "Khala.fKHALA_LaBr3_E_Energy", "Plastic.fPlastic_Time > 0");
  // hlabr3->Add(hfatima);
  // hlabr3->Add(hkhala);

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
    long long tpl_temp = 1.e12;
    for(int n = 0; n < Plastic->GetMult();n++){
      if (Plastic->GetTime(n) < tpl_temp)
	tpl_temp = Plastic->GetTime(n); // get the fastest time hit on the plastic detector
    }
    if (tpl_temp < 1.e12){ // check that we have a plastic hit in this event
      for(int n = 0; n < Fatima->GetFatimaLaBr3EMult();n++){
	if (fabs(Fatima->GetFatimaLaBr3TTime(n)-tpl_temp)<10)
	hlabr3->Fill(Fatima->GetFatimaLaBr3EEnergy(n));
      }
      for(int n = 0; n < Khala->GetKhalaLaBr3EMult();n++){
	if(fabs(Khala->GetKhalaLaBr3TTime(n)-tpl_temp)<10)
	hlabr3->Fill(Khala->GetKhalaLaBr3EEnergy(n)); 
      }
    }    
  }
}

void GetHPGeEnergySpectrum(TTree* tr, TH1D* hHPGe){
  tr->Project(hHPGe->GetName(), "Tigress.fTIG_Ge_Energy");
}
void GetHPGeEnergySpectrumWAS3ABiGate(TTree* tr, TH1D* hHPGe, double emax_WAS3ABi){
  tr->Project(hHPGe->GetName(), "Tigress.fTIG_Ge_Energy", Form("0.5*(WAS3ABi.fWAS3ABi_StripFront_Energy+WAS3ABi.fWAS3ABi_StripFront_Energy) < %.0lf", emax_WAS3ABi));
}

void GetHPGeEnergySpectrumTimeGate(TTree* tr, TH1D* hHPGe, double tmin, double tmax){
  tr->Project(hHPGe->GetName(), "Tigress.fTIG_Ge_Energy", Form("Tigress.fTIG_Ge_TimeLED > %e && Tigress.fTIG_Ge_TimeLED < %e", tmin, tmax));
}
void GetHPGeEnergySpectrumWithPlastic(TTree* tr, TH1D* hHPGe){
  tr->Project(hHPGe->GetName(), "Tigress.fTIG_Ge_Energy", "Plastic.fPlastic_Time > 0");
}
void start_stop_analysis_v2(){
  gStyle->SetOptStat(0);

  //  teu = new TChain("SimulatedTree");

  // 100Sn
  //  teu->Add("100Sn_decay_40.root");
  // // teu->Add("100Sn_decay_40_doubleTau.root");
  // TH1D* he_idaten = new TH1D("he_idaten", "", 1050, 0, 2100);
  // GetLaBr3EnergySpectrumWithPlastic(teu, he_idaten);
  // TH1D* he_hpge = new TH1D("he_hpge", "", 2100, 0, 2100);
  // GetHPGeEnergySpectrum(teu, he_hpge);  
  // TH1D* he_idaten_labr3_gate = new TH1D("he_idaten_labr3_gate", "", 1050, 0, 2100);
  // GetGammaProjectionLaBr3WithPlastic(teu, he_idaten_labr3_gate, 141,7);
  // TH1D* ht_idaten_plastic_1297 = new TH1D("ht_idaten_plastic_1297", "", 80, -2, 2);
  // GetStartStopHistWithPlastic(teu, ht_idaten_plastic_1297, 1297, 40);
  // TH1D* ht_idaten_plastic_436 = new TH1D("ht_idaten_plastic_436", "", 80, -2, 2);
  // GetStartStopHistWithPlastic(teu, ht_idaten_plastic_436, 436, 20);
  // TH1D* ht_idaten_plastic_141 = new TH1D("ht_idaten_plastic_141", "", 80, -2, 2);
  // GetStartStopHistWithPlastic(teu, ht_idaten_plastic_141, 141,7);
  
  // TH1D* ht_idaten_plastic_95 = new TH1D("ht_idaten_plastic_95", "", 80, -2, 2);
  // GetStartStopHistWithPlastic(teu, ht_idaten_plastic_95, 95, 8);
  // TH1D* ht_idaten_plastic_40 = new TH1D("ht_idaten_plastic_40", "", 80, -2, 2);
  // GetStartStopHistWithPlastic(teu, ht_idaten_plastic_40, 40, 5);
  // TH1D* ht_idaten_1297_141 = new TH1D("ht_idaten_1297_141", "", 80, -2, 2);
  // GetStartStopHist(teu, ht_idaten_1297_141, 1297, 40, 141,7);
  // TH1D* ht_idaten_436_141 = new TH1D("ht_idaten_436_141", "", 80, -2, 2);
  // GetStartStopHist(teu, ht_idaten_436_141, 436, 20, 141,7);
  // TH1D* ht_idaten_1297_95 = new TH1D("ht_idaten_1297_95", "", 80, -2, 2);
  // GetStartStopHist(teu, ht_idaten_1297_95, 1297, 40, 95, 7);
  // TH1D* ht_idaten_436_95 = new TH1D("ht_idaten_436_95", "", 80, -2, 2);
  // GetStartStopHist(teu, ht_idaten_436_95, 436, 20, 95, 7);
  // TH1D* ht_idaten_141_95 = new TH1D("ht_idaten_141_95", "", 80, -2, 2);
  // GetStartStopHist(teu, ht_idaten_141_95, 141,7, 95, 7);
  // TH1D* ht_idaten_1297_40 = new TH1D("ht_idaten_1297_40", "", 80, -2, 2);
  // GetStartStopHist(teu, ht_idaten_1297_40, 1297, 40, 40, 5);
  // TH1D* ht_idaten_436_40 = new TH1D("ht_idaten_436_40", "", 80, -2, 2);
  // GetStartStopHist(teu, ht_idaten_436_40, 436, 20, 40, 5);
  // TH1D* ht_idaten_141_40 = new TH1D("ht_idaten_141_40", "", 80, -2, 2);
  // GetStartStopHist(teu, ht_idaten_141_40, 141,7, 40, 5);
  // TH1D* ht_idaten_95_40 = new TH1D("ht_idaten_95_40", "", 80, -2, 2);
  // GetStartStopHist(teu, ht_idaten_95_40, 95, 7, 40, 5);

  // TH1D* ht_idaten_141_sum = (TH1D*)ht_idaten_1297_141->Clone("ht_iaten_141_sum");
  // ht_idaten_141_sum->Add(ht_idaten_436_141);

  // TH1D* ht_idaten_95_sum = (TH1D*)ht_idaten_1297_95->Clone("ht_iaten_95_sum");
  // ht_idaten_95_sum->Add(ht_idaten_436_95);
  // ht_idaten_95_sum->Add(ht_idaten_141_95);

  // TH1D* ht_idaten_40_sum = (TH1D*)ht_idaten_1297_40->Clone("ht_iaten_40_sum");
  // ht_idaten_40_sum->Add(ht_idaten_436_40);
  // ht_idaten_40_sum->Add(ht_idaten_141_40);
  // ht_idaten_40_sum->Add(ht_idaten_95_40);

  
  // //  TFile* f100Sn = new TFile("outputs_100Sn.root", "recreate");
  // TFile* f100Sn = new TFile("outputs_100Sn_doubleTau.root", "recreate");
  // he_idaten->Write();
  // he_hpge->Write();
  // he_idaten_labr3_gate->Write();
  // ht_idaten_plastic_1297->Write();
  // ht_idaten_plastic_436->Write();
  // ht_idaten_plastic_141->Write();
  // ht_idaten_plastic_95->Write();
  // ht_idaten_plastic_40->Write();
  
  // ht_idaten_1297_141->Write();
  // ht_idaten_436_141->Write();
  // ht_idaten_1297_95->Write();
  // ht_idaten_436_95->Write();
  // ht_idaten_141_95->Write();
  // ht_idaten_1297_40->Write();
  // ht_idaten_436_40->Write();
  // ht_idaten_141_40->Write();
  // ht_idaten_95_40->Write();

  // ht_idaten_141_sum->Write();
  // ht_idaten_95_sum->Write();
  // ht_idaten_40_sum->Write();
  
  // f100Sn->Close();

  // 96Cd
  // teu->Add("96Cd_isomer.root");
  // TH1D* he_idaten = new TH1D("he_idaten", "", 2000, 0, 2000);
  // GetLaBr3EnergySpectrum(teu, he_idaten);
  // TH1D* he_hpge = new TH1D("he_hpge", "", 2000, 0, 2000);
  // GetHPGeEnergySpectrum(teu, he_hpge);  
  // // TH1D* he_idaten_labr3_gate = new TH1D("he_idaten_labr3_gate", "", 2000, 0, 2000);
  // // GetGammaProjectionLaBr3(teu, he_idaten_labr3_gate, 788, 30);
  // TH1D* ht_idaten_307_419 = new TH1D("ht_idaten_307_419", "", 200, -3, 7);  
  // GetStartStopHist(teu, ht_idaten_307_419, 307, 20, 415, 15);
  // TH1D* ht_idaten_419_441 = new TH1D("ht_idaten_419_441", "", 200, -3, 7);
  // GetStartStopHist(teu, ht_idaten_419_441, 410, 15, 441, 5);
  // TH1D* ht_idaten_419_1026_1040 = new TH1D("ht_idaten_419_1026_1040", "", 200, -3, 7);
  // GetStartStopHist(teu, ht_idaten_419_1026_1040, 410, 15, 1035, 35);
  // TH1D* ht_idaten_419_811 = new TH1D("ht_idaten_419_811", "", 200, -3, 7);
  // GetStartStopHist(teu, ht_idaten_419_811, 410, 15, 811, 20);
  
  // TH1D* ht_idaten_1561_441 = new TH1D("ht_idaten_1561_441", "", 200, -3, 7);
  // GetStartStopHist(teu, ht_idaten_1561_441, 1561, 40, 441, 5);
  // TH1D* ht_idaten_1561_1026_1040 = new TH1D("ht_idaten_1561_1026_1040", "", 200, -3, 7);
  // GetStartStopHist(teu, ht_idaten_1561_1026_1040, 1561, 40, 1035, 35);
  // TH1D* ht_idaten_1561_811 = new TH1D("ht_idaten_1561_811", "", 200, -3, 7);
  // GetStartStopHist(teu, ht_idaten_1561_811, 1561, 40, 811, 20);

  // TH1D* ht_idaten_1104_441 = new TH1D("ht_idaten_1104_441", "", 200, -3, 7);
  // GetStartStopHist(teu, ht_idaten_1104_441, 1104, 30, 441, 5);
  // TH1D* ht_idaten_1104_1026_1040 = new TH1D("ht_idaten_1104_1026_1040", "", 200, -3, 7);
  // GetStartStopHist(teu, ht_idaten_1104_1026_1040, 1104, 30, 1035, 35);
  // TH1D* ht_idaten_1104_811 = new TH1D("ht_idaten_1104_811", "", 200, -3, 7);
  // GetStartStopHist(teu, ht_idaten_1104_811, 1104, 30, 811, 20);

  // TH1D* ht_idaten_457_441 = new TH1D("ht_idaten_457_441", "", 200, -3, 7);
  // GetStartStopHist(teu, ht_idaten_457_441, 467, 10, 441, 5);
  // TH1D* ht_idaten_457_1026_1040 = new TH1D("ht_idaten_457_1026_1040", "", 200, -3, 7);
  // GetStartStopHist(teu, ht_idaten_457_1026_1040, 467, 10, 1035, 35);
  // TH1D* ht_idaten_457_811 = new TH1D("ht_idaten_457_811", "", 200, -3, 7);
  // GetStartStopHist(teu, ht_idaten_457_811, 467, 10, 811, 20);

  
  // TFile* f96Cd = new TFile("outputs_96Cd.root", "recreate");
  // he_idaten->Write();
  // he_hpge->Write();
  // // he_idaten_labr3_gate->Write();
  // ht_idaten_307_419->Write();
  // ht_idaten_419_441->Write();
  // ht_idaten_419_1026_1040->Write();
  // ht_idaten_419_811->Write();
  
  // ht_idaten_1561_441->Write();
  // ht_idaten_1561_1026_1040->Write();
  // ht_idaten_1561_811->Write();
  
  // ht_idaten_1104_441->Write();
  // ht_idaten_1104_1026_1040->Write();
  // ht_idaten_1104_811->Write();
  
  // ht_idaten_457_441->Write();
  // ht_idaten_457_1026_1040->Write();
  // ht_idaten_457_811->Write();
  
  
  // f96Cd->Close();
  
  // 95Cd
  // teu->Add("95Cd_332_isomer.root");
  // TH1D* he_idaten = new TH1D("he_idaten", "", 1400, 0, 1400);
  // GetLaBr3EnergySpectrumTimeGate(teu, he_idaten, 0, 2e4);
  // TH1D* he_hpge = new TH1D("he_hpge", "", 1400, 0, 1400);
  // GetHPGeEnergySpectrumTimeGate(teu, he_hpge, 0, 2e4); 
  // TFile* f95Cd = new TFile("outputs_95Cd_332_realistic.root", "recreate");
  // he_idaten->Write();
  // he_hpge->Write();
  
  // f95Cd->Close();

  // teu->Add("95Cd_232_isomer.root");
  // TH1D* he_idaten = new TH1D("he_idaten", "", 1400, 0, 1400);
  // GetLaBr3EnergySpectrumWAS3ABiGate(teu, he_idaten, 500);
  // TH1D* he_hpge = new TH1D("he_hpge", "", 1400, 0, 1400);
  // GetHPGeEnergySpectrumWAS3ABiGate(teu, he_hpge, 500);
  // TH2D* me_idaten_hpge = new TH2D("me_idaten_hpge", "", 1400, 0, 1400, 1400, 0, 1400);
  // GetGammaGammaMatrix(teu, me_idaten_hpge);
  // // TH1D* he_idaten_labr3_gate = new TH1D("he_idaten_labr3_gate", "", 2000, 0, 2000);
  // // GetGammaProjectionLaBr3(teu, he_idaten_labr3_gate, 788, 30);
  
  // TFile* f95Cd = new TFile("outputs_95Cd_232.root", "recreate");
  // he_idaten->Write();
  // he_hpge->Write();
  // me_idaten_hpge->Write();
  // f95Cd->Close();

  // teu->Add("95Cd_12_isomer.root");
  // TH1D* he_idaten = new TH1D("he_idaten", "", 1400, 0, 1400);
  // GetLaBr3EnergySpectrumWAS3ABiGate(teu, he_idaten, 300);
  // TH1D* he_hpge = new TH1D("he_hpge", "", 1400, 0, 1400);
  // GetHPGeEnergySpectrumWAS3ABiGate(teu, he_hpge, 300); 
  // TFile* f95Cd = new TFile("outputs_95Cd_12.root", "recreate");
  // he_idaten->Write();
  // he_hpge->Write();
  
  // f95Cd->Close();

  //99In for 99Cd
  /*
  bool GAMMA_GAMMA = 0;
   teu->Add("99In_decay.root");
  TH1D* he_idaten = new TH1D("he_idaten", "", 400, 0, 2000);
  GetLaBr3EnergySpectrumWithPlastic(teu, he_idaten);
  TH1D* he_hpge = new TH1D("he_hpge", "", 2000, 0, 2000);
  GetHPGeEnergySpectrum(teu, he_hpge);  
  TH1D* ht_idaten_1234_441;
  TF1* fhalf;
  if (GAMMA_GAMMA){
    ht_idaten_1234_441 = new TH1D("ht_idaten_1234_441", "", 35, -0.5, 3);  
    GetStartStopHist(teu, ht_idaten_1234_441, 1234, 30, 441, 20);
    TH1D* ht_idaten_607_441 = new TH1D("ht_idaten_607_441", "", 35, -0.5, 3);
    GetStartStopHist(teu, ht_idaten_607_441, 607, 17, 441, 20);
    ht_idaten_1234_441->Add(ht_idaten_607_441);
    fhalf = new TF1("fhalf", "[1]*(exp(0.5*(log(2)/[0])*(2*[3]+(log(2)/[0])*[2]**2-2*x))*TMath::Erfc(([3]+(log(2)/[0])*[2]**2-x)/(sqrt(2)*[2])))+[4]", -0.4, 3);
  }
  else{
    ht_idaten_1234_441 = new TH1D("ht_idaten_1234_441", "", 55, -0.5, 5);  
    GetStartStopHistWithPlastic(teu, ht_idaten_1234_441, 441, 20);
    fhalf = new TF1("fhalf", "[1]*0.5**(x/[0])+[2]", 0.4, 5);
  }

  fhalf->SetParameter(0, 0.33);
  fhalf->SetParameter(1, 25);
  fhalf->SetParameter(2, 1);
  fhalf->SetParameter(3, 0);
  fhalf->FixParameter(4, 0);
  //  fhalf->SetParLimits(4, 0, 100);
  //  ht_idaten_1234_441->Rebin(5);  
  ht_idaten_1234_441->Fit(fhalf, "RN0");
  
  TFile* f99In;
  if (GAMMA_GAMMA){
    f99In = new TFile("outputs_99In_decay_gg.root", "recreate");
  }
  else{
    f99In = new TFile("outputs_99In_decay_single.root", "recreate");
  }
  he_idaten->Write();
  he_hpge->Write();
  ht_idaten_1234_441->Write();
  fhalf->Write();
  
  f99In->Close();
  TCanvas* cc = new TCanvas("cc", "cc", 1200, 800);
  cc->cd()->SetLogy(1);
  ht_idaten_1234_441->Draw();
  ht_idaten_1234_441->GetXaxis()->SetTitle("Time (ns)");
  ht_idaten_1234_441->GetYaxis()->SetTitle("Counts / 0.1 ns");
  fhalf->Draw("same");

  TLatex* latex = new TLatex();
  
 
  if (GAMMA_GAMMA){
    latex->DrawLatex(0.8, ht_idaten_1234_441->GetMaximum()*0.9, "^{99}Cd, 7/2^{+}; 1234(start)-441(stop)");
    latex->DrawLatex(1.1, ht_idaten_1234_441->GetMaximum()*0.45, "T_{1/2} (input)  = 330 ps");
    //  latex->DrawLatex(2, 150, Form("T_{1/2} = %.0lf #pm %.0lf ps", fhalf->GetParameter(0)*1000, fhalf->GetParError(0)*1000));
    latex->DrawLatex(1.1, ht_idaten_1234_441->GetMaximum()*0.3, Form("T_{1/2} (measured) = %.0lf(%.0lf) ps", fhalf->GetParameter(0)*1000, fhalf->GetParError(0)*1000));
    cc->SaveAs("99In_decay_gg.png");
    cc->SaveAs("99In_decay_gg.pdf");
    cc->SaveAs("99In_decay_gg.eps");
  }
  else{
    latex->DrawLatex(0.8, ht_idaten_1234_441->GetMaximum()*0.9, "^{99}Cd, 7/2^{+}; plastic(start)-441(stop)");
    latex->DrawLatex(1.1, ht_idaten_1234_441->GetMaximum()*0.45, "T_{1/2} (input)  = 330 ps");
    //  latex->DrawLatex(2, 150, Form("T_{1/2} = %.0lf #pm %.0lf ps", fhalf->GetParameter(0)*1000, fhalf->GetParError(0)*1000));
    latex->DrawLatex(1.1, ht_idaten_1234_441->GetMaximum()*0.25, Form("T_{1/2} (measured) = %.0lf(%.0lf) ps", fhalf->GetParameter(0)*1000, fhalf->GetParError(0)*1000));
    cc->SaveAs("99In_decay_single.png");
    cc->SaveAs("99In_decay_single.pdf");
    cc->SaveAs("99In_decay_single.eps");
  }
  */
 //97Ag for 97Pd (gates on 1257, 1448 and 1490
  bool GAMMA_GAMMA = 0;
   teu->Add("97Ag_decay.root");
  TH1D* he_idaten = new TH1D("he_idaten", "", 400, 0, 2000);
  GetLaBr3EnergySpectrumWithPlastic(teu, he_idaten);
  TH1D* he_hpge = new TH1D("he_hpge", "", 2000, 0, 2000);
  GetHPGeEnergySpectrum(teu, he_hpge);  
  TH1D* ht_idaten_1257_686;
  TF1* fhalf;
  if (GAMMA_GAMMA){
    ht_idaten_1257_686 = new TH1D("ht_idaten_1257_686", "", 50, -0.5, 2);
    GetStartStopHist(teu, ht_idaten_1257_686, 1257, 30, 686, 25);
    TH1D* ht_idaten_1312_686 = new TH1D("ht_idaten_1312_686", "", 50, -0.5, 2);
    GetStartStopHist(teu, ht_idaten_1312_686, 1312, 25, 686, 25);
    TH1D* ht_idaten_1448_1490_686 = new TH1D("ht_idaten_1448_686", "", 50, -0.5, 2);
    GetStartStopHist(teu, ht_idaten_1448_1490_686, 1465, 55, 686, 25);
    ht_idaten_1257_686->Add(ht_idaten_1312_686);
    ht_idaten_1257_686->Add(ht_idaten_1448_1490_686);
    fhalf = new TF1("fhalf", "[1]*(exp(0.5*(log(2)/[0])*(2*[3]+(log(2)/[0])*[2]**2-2*x))*TMath::Erfc(([3]+(log(2)/[0])*[2]**2-x)/(sqrt(2)*[2])))+[4]", -0.5, 3);
  }
  else{
    ht_idaten_1257_686 = new TH1D("ht_idaten_1257_686", "", 70, -0.5, 3);  
    GetStartStopHistWithPlastic(teu, ht_idaten_1257_686, 686, 25);
    fhalf = new TF1("fhalf", "[1]*0.5**(x/[0])+[2]", 0.4, 3);
  }

  fhalf->SetParameter(0, 0.09);
  fhalf->SetParameter(1, 25);
  fhalf->SetParameter(2, 1);
  fhalf->SetParameter(3, 0);
  fhalf->FixParameter(4, 0);
  //  fhalf->SetParLimits(4, 0, 100);
  //  ht_idaten_1257_686->Rebin(5);  
  ht_idaten_1257_686->Fit(fhalf, "RN0");
  
  TFile* f97Ag;
  if (GAMMA_GAMMA){
    f97Ag = new TFile("outputs_97Ag_decay_gg.root", "recreate");
  }
  else{
    f97Ag = new TFile("outputs_97Ag_decay_single.root", "recreate");
  }
  he_idaten->Write();
  he_hpge->Write();
  ht_idaten_1257_686->Write();
  fhalf->Write();
  
  f97Ag->Close();
  TCanvas* cc = new TCanvas("cc", "cc", 1200, 800);
  cc->cd()->SetLogy(1);
  ht_idaten_1257_686->Draw();
  ht_idaten_1257_686->GetXaxis()->SetTitle("Time (ns)");
  ht_idaten_1257_686->GetYaxis()->SetTitle("Counts / 0.1 ns");
  fhalf->Draw("same");

  TLatex* latex = new TLatex();
  
 
  if (GAMMA_GAMMA){
    latex->DrawLatex(0.6, ht_idaten_1257_686->GetMaximum()*0.9, "^{97}Pd, 7/2^{+}; 1257(start)-686(stop)");
    latex->DrawLatex(0.8, ht_idaten_1257_686->GetMaximum()*0.45, "T_{1/2} (input)  = 90 ps");
    //  latex->DrawLatex(2, 150, Form("T_{1/2} = %.0lf #pm %.0lf ps", fhalf->GetParameter(0)*1000, fhalf->GetParError(0)*1000));
    latex->DrawLatex(0.8, ht_idaten_1257_686->GetMaximum()*0.3, Form("T_{1/2} (measured) = %.0lf(%.0lf) ps", fhalf->GetParameter(0)*1000, fhalf->GetParError(0)*1000));
    cc->SaveAs("97Ag_decay_gg.png");
    cc->SaveAs("97Ag_decay_gg.pdf");
    cc->SaveAs("97Ag_decay_gg.eps");
  }
  else{
    latex->DrawLatex(0.8, ht_idaten_1257_686->GetMaximum()*0.9, "^{97}Pd, 7/2^{+}; plastic(start)-686(stop)");
    latex->DrawLatex(1.1, ht_idaten_1257_686->GetMaximum()*0.45, "T_{1/2} (input)  = 90 ps");
    //  latex->DrawLatex(2, 150, Form("T_{1/2} = %.0lf #pm %.0lf ps", fhalf->GetParameter(0)*1000, fhalf->GetParError(0)*1000));
    latex->DrawLatex(1.1, ht_idaten_1257_686->GetMaximum()*0.25, Form("T_{1/2} (measured) = %.0lf(%.0lf) ps", fhalf->GetParameter(0)*1000, fhalf->GetParError(0)*1000));
    cc->SaveAs("97Ag_decay_single.png");
    cc->SaveAs("97Ag_decay_single.pdf");
    cc->SaveAs("97Ag_decay_single.eps");
  }
  
}
