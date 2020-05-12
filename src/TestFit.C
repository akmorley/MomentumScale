#include "TROOT.h"



#include "TH3.h"
#include "TH2.h"
#include "TH1.h"
#include "TProfile2D.h"
#include "TFile.h"
#include "TCanvas.h"

#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooMomentMorph.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"

#include "RooAbsPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"



#include "HelperFunctions.h"
#include "HelperClasses.h"

using namespace RooFit;
void TestFit();
void TestFitSimple();

int main()
{
  TestFitSimple();
  return 0 ;
}


void TestFit()
{
  bool isZ = true ;
  bool fitBackground =true;
  bool drawPlots = true;
  bool fitData = true;
    // How the data is cut up
  int nBinsEta = 3;
  int nBinsPhi = 5;
  int nAbins = 10; // dz bins
  int nMbins = 5;  // dS bins
  //TFile infile("JPsiMass2.root");
  //TFile infileData("JPsiMass_Data.root");

  TFile infile("ZMass-JN.root");
  TFile infileData("ZMass-JN-400.root");



  SimFitConfig  simFitConfig( fitBackground, fitData, drawPlots, false );
  
  auto w = new RooWorkspace("w");
  RooRealVar* x;
  if(isZ){
    x = (RooRealVar*) w->factory("x[4000,14400]");
  }  else{
    x = (RooRealVar*) w->factory("x[6,18]");
  }
  x->SetTitle("m^{2}_{#mu#mu} [GeV^{2}]");

  for( int i(1); i <= nBinsPhi; ++i){
    for( int j(3); j <= nBinsEta; ++j){
      w->factory(Form("deltaS_%i_%i[0,-0.1,0.1]",j,i));
      w->factory(Form("deltadZ_%i_%i[0,-0.1,0.1]",j,i));      
    }
  }


  for(int eta1(3); eta1 <= nBinsEta; ++eta1){ 
    for(int eta2(eta1); eta2 <= nBinsEta; ++eta2){ 
      for(int phi1(1); phi1 <= nBinsPhi; ++phi1){ 
        for(int phi2(phi1); phi2 <= nBinsPhi; ++phi2){ 
          TString nameSub = Form( "%i_%i_%i_%i", eta1, eta2, phi1, phi2 );
          auto Mass2AA_MC = (TH3*)infile.Get("Mass2AA_"+nameSub);
          Mass2AA_MC->SetName("Reference_"+nameSub);
          auto Mass2AA_Data = (TH3*) infileData.Get("Mass2AA_"+nameSub);

          std::cout << nameSub << " : " <<  Mass2AA_MC->GetEntries() << " " << Mass2AA_Data->GetEntries() << std::endl;

          // Histgrams containing the average value of each quantity in each bin
          auto Az1_Hist = (TProfile2D*) infile.Get("Az1_"+nameSub);
          auto Az2_Hist = (TProfile2D*) infile.Get("Az2_"+nameSub);
          auto As1_Hist = (TProfile2D*) infile.Get("As1_"+nameSub);
          auto As2_Hist = (TProfile2D*) infile.Get("As2_"+nameSub);




          // Variables we are going to fit
          auto deltaS1  = (RooRealVar*) w->var(Form("deltaS_%i_%i",eta1,phi1)); //w->factory("deltaS1[0,-0.1,0.1]");
          auto deltadZ1 = (RooRealVar*)  w->var(Form("deltadZ_%i_%i",eta1,phi1));//w->factory("deltadZ1[0,-0.1,0.1]");
          auto deltaS2  = (RooRealVar*) w->var(Form("deltaS_%i_%i",eta2,phi2)); //w->factory("deltaS1[0,-0.1,0.1]");
          auto deltadZ2 = (RooRealVar*)  w->var(Form("deltadZ_%i_%i",eta2,phi2));//w->factory("deltadZ1[0,-0.1,0.1]");
          //auto deltaS2  = (RooRealVar*) w->factory("deltaS2[0,-0.1,0.1]");
          //auto deltadZ2 = (RooRealVar*) w->factory("deltadZ2[0,-0.1,0.1]");



          RooArgList argList( *deltaS1, *deltaS2, *deltadZ1, *deltadZ2);
          //RooArgList argList( *deltaS1, *deltaS1, *deltadZ1, *deltadZ1);

          w->defineSet("obs","x"); //observables
          w->defineSet("poi",argList); //parameters of interest

          ScaleFitData  sfData( Mass2AA_Data, Mass2AA_MC, As1_Hist, As2_Hist, Az1_Hist, Az2_Hist, nAbins, nMbins, isZ);
          sfData.entryThreshold     = 1000;
          sfData.considerThreshold  = 100;

          ScaleFit scaleFit( w, nameSub, sfData);
          scaleFit.constructModel( x, argList, simFitConfig);      
        
        
        }
      }
    }
  }


  
  



  // Associate model with the physics state and model_ctl with the control state
  RooSimultaneous simPdf("simPdf", "simultaneous pdf", simFitConfig.pdfMap, simFitConfig.categories );
  RooDataHist combData("combData", "combData", *x, simFitConfig.categories, simFitConfig.dhistMap);
  simPdf.fitTo( combData );
  
  //auto nps = w->set("nuisParams");
  //for(auto np : *nps)
  //  w->var( np->GetName() )->setConstant(false);
  //simPdf.fitTo( combData );

  return;

  
}



void TestFitSimple()
{
  bool isZ = true ;
  bool fitBackground =true;
  bool drawPlots = true;
  bool fitData = true;
    // How the data is cut up
  int nAxbins = 25; // dz bins
  int nAybins = 200; // dz bins
  int nMbins = 1;  // dS bins
  //TFile infile("JPsiMass2.root");
  //TFile infileData("JPsiMass_Data.root");

  TFile infile("ZMass-JN-25-200.root");
  TFile infileData("ZMass-JN-400-25-200.root");

  TH3D* Mass2AA_MC, *Mass2AA_Data;
  TProfile2D* Az1_Hist, *Az2_Hist, *As1_Hist, *As2_Hist;

  int intputHits(0);
  for(int phi1(1); phi1 <= 5; ++phi1){ 
    for(int phi2(phi1); phi2 <= 5; ++phi2){ 
      TString nameSub = Form( "%i_%i_%i_%i", 3, 3, phi1, phi2 );
      if(intputHits==0){
        Mass2AA_MC = (TH3D*)infile.Get("Mass2AA_"+nameSub);
        Mass2AA_MC->SetName("Reference_"+nameSub);
        Mass2AA_Data = (TH3D*) infileData.Get("Mass2AA_"+nameSub);

        std::cout << nameSub << " : " <<  Mass2AA_MC->GetEntries() << " " << Mass2AA_Data->GetEntries() << std::endl;

        // Histgrams containing the average value of each quantity in each bin
        Az1_Hist = (TProfile2D*) infile.Get("Az1_"+nameSub);
        Az2_Hist = (TProfile2D*) infile.Get("Az2_"+nameSub);
        As1_Hist = (TProfile2D*) infile.Get("As1_"+nameSub);
        As2_Hist = (TProfile2D*) infile.Get("As2_"+nameSub);
      } else {
        Mass2AA_MC->Add( (TH3*)infile.Get("Mass2AA_"+nameSub) );
        Mass2AA_Data->Add( (TH3*) infileData.Get("Mass2AA_"+nameSub) );

        std::cout << nameSub << " : " <<  Mass2AA_MC->GetEntries() << " " << Mass2AA_Data->GetEntries() << std::endl;

        // Histgrams containing the average value of each quantity in each bin
        Az1_Hist->Add( (TProfile2D*) infile.Get("Az1_"+nameSub) );
        Az2_Hist->Add( (TProfile2D*) infile.Get("Az2_"+nameSub) );
        As1_Hist->Add( (TProfile2D*) infile.Get("As1_"+nameSub) );
        As2_Hist->Add( (TProfile2D*) infile.Get("As2_"+nameSub) );
      }
      ++intputHits;
    }
  }

  SimFitConfig  simFitConfig( fitBackground, fitData, drawPlots, false );
  
  auto w = new RooWorkspace("w");
  RooRealVar* x;
  if(isZ){
    x = (RooRealVar*) w->factory("x[4000,14400]");
  }  else{
    x = (RooRealVar*) w->factory("x[6,18]");
  }
  x->SetTitle("m^{2}_{#mu#mu} [GeV^{2}]");

  w->factory("deltaS[0,-0.1,0.1]");
  w->factory("deltadZ[0,-0.1,0.1]");      

  int  eta1, eta2, phi1, phi2; 
  eta1 = eta2 = phi1 = phi2 =0;
  
  TString nameSub = Form( "%i_%i_%i_%i", eta1, eta2, phi1, phi2 );
 /*
  auto Mass2AA_MC = (TH3*)infile.Get("Mass2AA");
  Mass2AA_MC->SetName("Reference_"+nameSub);
  auto Mass2AA_Data = (TH3*)infile.Get("Mass2AA_Pt");
//(TH3*) infileData.Get("Mass2AA");

  std::cout << nameSub << " : " <<  Mass2AA_MC->GetEntries() << " " << Mass2AA_Data->GetEntries() << std::endl;

  // Histgrams containing the average value of each quantity in each bin
  auto Az1_Hist = (TProfile2D*) infile.Get("Az1");
  auto Az2_Hist = (TProfile2D*) infile.Get("Az2");
  auto As1_Hist = (TProfile2D*) infile.Get("As1");
  auto As2_Hist = (TProfile2D*) infile.Get("As2");
*/

  // Variables we are going to fit
  auto deltaS  = (RooRealVar*) w->var("deltaS"); 
  auto deltadZ = (RooRealVar*)  w->var("deltadZ");
  RooArgList argList( *deltaS, *deltaS, *deltadZ, *deltadZ);

  w->defineSet("obs","x"); //observables
  w->defineSet("poi",argList); //parameters of interest

  ScaleFitData  sfData( Mass2AA_Data, Mass2AA_MC, As1_Hist, As2_Hist, Az1_Hist, Az2_Hist, nAxbins, nAybins, nMbins, isZ);
  sfData.entryThreshold     = 5000;
  sfData.considerThreshold  = 100;

  ScaleFit scaleFit( w, nameSub, sfData);
  scaleFit.constructModel( x, argList, simFitConfig);      
      


  // Associate model with the physics state and model_ctl with the control state
  RooSimultaneous simPdf("simPdf", "simultaneous pdf", simFitConfig.pdfMap, simFitConfig.categories );
  RooDataHist combData("combData", "combData", *x, simFitConfig.categories, simFitConfig.dhistMap);
  simPdf.fitTo( combData );
  
  //auto nps = w->set("nuisParams");
  //for(auto np : *nps)
  //  w->var( np->GetName() )->setConstant(false);
  //simPdf.fitTo( combData );

  return;

  
}
