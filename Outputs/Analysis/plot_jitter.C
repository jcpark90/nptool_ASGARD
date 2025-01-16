void plot_jitter(){
  TGraph* gjitter = new TGraph("fatima_jitter.txt", "%lg %lg");
  gjitter->SetMarkerStyle(20);
  gjitter->Draw("alp");
  TF1* f1 = new  TF1("f1", "[0]/sqrt(x+[1])+[2]", 30, 3000);
  f1->SetParameter(0, 3000);
  f1->SetParameter(1, 50);
  //f1->SetParLimits(1, 0, 1000);
  f1->SetParameter(2, 230);
  
  gjitter->Fit(f1, "QRN0");
  f1->Draw("same");

  cout<<f1->GetParameter(0)<<endl;
  cout<<f1->GetParameter(1)<<endl;
  cout<<f1->GetParameter(2)<<endl;
  
}
