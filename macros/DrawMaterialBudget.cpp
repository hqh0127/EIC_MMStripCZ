void DrawMaterialBudget(const char* name){
	TGeoManager::Import(name);
	string filename= name;
	const size_t period_idx = filename.rfind('.');
  if (std::string::npos != period_idx)
  {
    filename.erase(period_idx);
  }
	string output_plot_name=filename+"_MB.pdf";
	TCanvas *cradlen = new TCanvas("radlen", "radiation length");
	cradlen->SetRightMargin(0.18);
  TH2F *X0 = gGeoManager->GetTopVolume()->LegoPlot();
  X0->SetXTitle( "#varphi (deg)");
  X0->SetYTitle( "#theta (deg)");
  X0->SetZTitle( "X/X0");
  X0->GetZaxis()->SetTitleOffset(1.4);
  
	gStyle->SetPalette(kInvertedDarkBodyRadiator);
  X0->SetStats(kFALSE);
  X0->Draw("colz");
  cradlen->Update();
  TPaletteAxis *palette = (TPaletteAxis*)X0->GetListOfFunctions()->FindObject("palette");
  double y1 = palette->GetY1NDC();
  double y2 = palette->GetY2NDC();
  double x1 = palette->GetX1NDC();
  double x2 = palette->GetX2NDC();
  palette->SetX1NDC(x1+0.02);
  palette->SetX2NDC(x2+0.02);
  cradlen->Print(output_plot_name.data());
  output_plot_name = filename+"_MB.png";
  cradlen->Print(output_plot_name.data());
  new TCanvas;
  gGeoManager->GetTopNode()->Draw("ogl");
}
