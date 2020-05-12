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
void TextMorph();

int main()
{
  TextMorph();
  return 0 ;
}


void TextMorph()
{
  //gROOT->ProcessLine(".L ./HelperClasses.cxx+");
  //gROOT->ProcessLine(".L ./RooTwoSidedCBShape.cxx+");
  //gROOT->ProcessLine(".L ./RooMyPoly.cxx+");

  bool isZ = true ;
  bool fitBackground =true;
  bool drawPlots = true;
  bool fitData = true;

  //TFile infile("JPsiMass2.root");
  //TFile infileData("JPsiMass_Data.root");

  TFile infile("ZMass_MC.root");
  TFile infileData("ZMass_Data.root");

  int aNbins = 10;
  int mNbins = 5;
  int mmNbins = 5;
  
  VariableAxis axisA( -5, 5, aNbins);
  VariableAxis axisMA( 0.95125, 1.00125, mNbins);
  VariableAxis axisDA( 0., 1.0, mmNbins);


  auto Mass2AA = (TH3*)infile.Get("Mass2AA");
  Mass2AA->SetName("Reference");
  auto Mass2AA_Data = (TH3*) infileData.Get("Mass2AA");

  auto Az1_Hist = (TProfile2D*) infile.Get("Az1");
  auto Az2_Hist = (TProfile2D*) infile.Get("Az2");
  auto As1_Hist = (TProfile2D*) infile.Get("As1");
  auto As2_Hist = (TProfile2D*) infile.Get("As2");


  TH2* rebinTest =  (TH2*)Mass2AA->Project3D("yx");
  double entryThreshold = 100000;
  double considerThreshold = 30000;
  

  BinMergeHelper binMerger( aNbins, mNbins);

  int binXmax = rebinTest->GetNbinsX();
  int binYmax = rebinTest->GetNbinsY();
  for( int binX(0); binX < binXmax; ++binX ){
    for( int binY(0); binY < binYmax; ++binY ){
      auto val = rebinTest->GetBinContent( binX+1, binY+1 );
      if(val < considerThreshold){
        rebinTest->SetBinContent( binX+1, binY+1 , 0);
        continue;
      }
      binMerger.bins.insert( { {binX,binY}, BinContent( binX, binY, val) }  );
    }
  }

  binMerger.mergeBins(entryThreshold);
  binMerger.printBins();
  binMerger.mergeHistBins( rebinTest );
  


  auto data3D = fitData ? Mass2AA_Data : (TH3*) infile.Get("Mass2AA_P");
  auto mc3D = Mass2AA;

  std::vector<double> biases = {-0.04,-0.004,+0.004,+0.04};
  std::vector<TH1*> hists;

  auto w = new RooWorkspace("w");
  RooRealVar* x;
  if(isZ){
    x = (RooRealVar*) w->factory("x[4000,14400]");
  }  else{
    x = (RooRealVar*) w->factory("x[6,18]");
  }
  x->SetTitle("m^{2}_{#mu#mu} [GeV^{2}]");

  auto deltaS1 = (RooRealVar*) w->factory("deltaS1[0,-0.1,0.1]");
  auto deltadZ1 = (RooRealVar*) w->factory("deltadZ1[0,-0.1,0.1]");
  auto deltaS2 = (RooRealVar*) w->factory("deltaS2[0,-0.1,0.1]");
  auto deltadZ2 = (RooRealVar*) w->factory("deltadZ2[0,-0.1,0.1]");

  RooCategory categories("AllCategories","AllCategories");
  std::map< std::string, RooDataHist *>   dhistMap;
  std::map< std::string, RooAbsPdf *>     pdfMap;
  RooArgSet nps;

    std::vector<std::vector<double>> constantsUsed;
  
  //3.686097^2 - 3.096900^2  = 3.9965214834
  w->factory("delta2S[3.9965214834]");
   

  TCanvas* cValid  =  new TCanvas("cValid","cValid",1600,800);
  cValid->Divide(2,1);
  TString PrintToFile("TestFits.pdf");
  cValid->Print(PrintToFile+"[");

  
  int counter(0);
  for( auto&  binPair : binMerger.bins ){
      auto& bin = binPair.second;
      int binX = bin.binX;
      int binY = bin.binY;
      
      std::cout << binX << " "<< binY << " Bin "<< counter << " weight : " << bin.totalWeight <<  std::endl;
  

      auto temp1D = bin.getProjection( mc3D ); 
      temp1D->SetName("MCProj");
      if(temp1D->GetEntries()<50000){
        delete temp1D;
        continue;
      }

      std::cout << binX << " "<< binY << " Bin "<< counter << " has entries : " <<  temp1D->GetEntries()  << std::endl;
      delete temp1D;
      ABins A1 = binXtoArAz( axisA, axisMA, binX );
      ABins A2 = binYtoArAz( axisDA, axisMA, binX, A1.Ar );
      
      //  deltaM/M = ar_1 * deltaR + ar_2 * deltaR + az_1 * deltaZ + az_2 * deltaZ

      

      auto as_1 = (RooRealVar*) w->factory(Form("as_1_%i[%f]",counter,2.*bin.getValue(As1_Hist))); //As1_Hist->GetBinContent(binX+1,binY+1)));
      auto as_2 = (RooRealVar*) w->factory(Form("as_2_%i[%f]",counter,2.*bin.getValue(As2_Hist))); //As2_Hist->GetBinContent(binX+1,binY+1)));
     
      auto adz_1 = (RooRealVar*) w->factory(Form("adz_1_%i[%f]",counter,2.*bin.getValue(Az1_Hist))); //Az1_Hist->GetBinContent(binX+1,binY+1)));
      auto adz_2 = (RooRealVar*) w->factory(Form("adz_2_%i[%f]",counter,2.*bin.getValue(Az2_Hist))); //Az2_Hist->GetBinContent(binX+1,binY+1)));

      
      
      auto deltax = (RooAbsReal*) w->factory( Form("expr::deltax_%i('@0*@4+@1*@6+@2*@5+@3*@7',{as_1_%i,as_2_%i,adz_1_%i,adz_2_%i,deltaS1,deltadZ1,deltaS1,deltadZ1} )", counter,counter,counter,counter,counter));
  
      //auto deltax = (RooAbsReal*) w->factory( Form("deltax_%i[0,-0.05,0.05]", counter));
      
      if(isZ)
        w->factory(Form("RooTwoSidedCBShape::pdfDCB_ctl_%i(%s, muCB_ctl_%i[8281,7500,9000], sigmaCB_ctl_%i[1000,100,5000], alphaCBLo_%i[1,0.1,10], nCBLo_%i[2,0.01,100], alphaCBHi_%i[1,0.1,10], nCBHi_%i[2,0.01,100])", counter, x->GetName(),counter, counter,counter,counter,counter, counter ) );
      else
        w->factory(Form("RooTwoSidedCBShape::pdfDCB_ctl_%i(%s, muCB_ctl_%i[9.6,8,10], sigmaCB_ctl_%i[1,0.01,4], alphaCBLo_%i[1,0.1,10], nCBLo_%i[2,0.01,100], alphaCBHi_%i[1,0.1,10], nCBHi_%i[2,0.01,100])", counter, x->GetName(),counter, counter,counter,counter,counter, counter ) );
      
      w->factory(Form("RooTwoSidedCBShape::pdfDCB_%i(%s, expr::muCB_%i('@0*(1.+@1)',  {muCB_ctl_%i,deltax_%i}), prod::sigmaCB_%i(sigmaCB_ctl_%i,sigmaCB_scale_%i[1.,0.5,5]), alphaCBLo_%i, nCBLo_%i, alphaCBHi_%i, nCBHi_%i)", counter, x->GetName(), counter, counter,counter,counter,counter, counter, counter,counter,counter,counter  ) );
      
      w->factory(Form("RooTwoSidedCBShape::pdfDCB_2S_%i(%s, expr::muCB_2S_%i('@0+@1', {muCB_ctl_%i,delta2S}), prod::sigmaCB_2S_%i(sigmaCB_ctl_%i,sigmaCB_2S_scale_%i[1.,0.5,5]), alphaCBLo_%i, nCBLo_%i, alphaCBHi_%i, nCBHi_%i)", counter, x->GetName(), counter, counter, counter, counter, counter, counter, counter, counter, counter ) );
 //RooPolynomial
      w->factory(Form("RooPolynomial::pdf_bkg_%i(%s, {const_a_%i[0,-1000,1000],const_b_%i[0,-1000,1000]} )", counter, x->GetName(), counter, counter ) );
      w->factory(Form("SUM::pdf_data_%i(coeff_sgnl_%i[1,0.1,1]*pdfDCB_%i, coeff_2S_%i[0.0,0,0.2]*pdfDCB_2S_%i, pdf_bkg_%i)",counter,counter,counter,counter,counter,counter));

      auto pdfDCB_ctl = w->pdf(Form("pdfDCB_ctl_%i", counter));
      auto pdfData = w->pdf(Form("pdf_data_%i", counter));
      
      TH1* mc1D = bin.getProjection( mc3D );
      mc1D->SetName("MCProj");
      RooDataHist* dataDataHist_ctl = new RooDataHist(Form("DataHist_ctl_%i",counter), Form("DataHist_%i",counter), RooArgList(*x), mc1D );
      delete mc1D;
      
       
      auto fitOk = pdfDCB_ctl->fitTo( *dataDataHist_ctl, PrintLevel(-1), Save() );
      if( !fitOk || fitOk->status() > 1 )
      {  
        delete dataDataHist_ctl;
        continue;
      }

      pdfMap.insert( {Form("Bin_ctl_%i",counter),  pdfDCB_ctl} );
      pdfMap.insert( {Form("Bin_%i",counter),  pdfData} );
      dhistMap.insert( {Form("Bin_ctl_%i",counter),  dataDataHist_ctl} );

      w->var(Form("alphaCBLo_%i",counter))->setConstant();
      w->var(Form("nCBLo_%i",counter))->setConstant();
      w->var(Form("alphaCBHi_%i",counter))->setConstant();
      w->var(Form("nCBHi_%i",counter))->setConstant();
      w->var(Form("sigmaCB_ctl_%i",counter))->setConstant();
      w->var(Form("muCB_ctl_%i",counter))->setConstant();
      w->var(Form("sigmaCB_scale_%i",counter))->setConstant();

      nps.add(  *w->var(Form("muCB_ctl_%i",counter)) );
      nps.add(  *w->var(Form("sigmaCB_ctl_%i",counter)) );
   
      constantsUsed.push_back( {as_1->getVal(),as_2->getVal(),adz_1->getVal(),adz_2->getVal(), as_1->getVal()+as_2->getVal(),adz_1->getVal()+adz_2->getVal()});
      
      //constantsUsed.push_back( {as_1->getVal(),as_2->getVal(),adz_1->getVal(),adz_2->getVal(), as_1->getVal()+as_2->getVal(),adz_1->getVal()+adz_2->getVal(),w->var(Form("deltax_%i",counter))->getVal()});
      
      
      auto frameFit_ctl  = x->frame();
      if(drawPlots){
        dataDataHist_ctl->plotOn(frameFit_ctl);
        pdfDCB_ctl->plotOn(frameFit_ctl, RooFit::LineColor(kBlue));
      }
  

      //w->import(deltax);
      categories.defineType(Form("Bin_%i",counter));
      categories.defineType(Form("Bin_ctl_%i",counter));
      

        
      TH1* data1D = bin.getProjection( data3D );
      data1D->SetName("Data");
      RooDataHist* dataDataHist = new RooDataHist(Form("DataHist_%i",counter), Form("DataHist_%i",counter), RooArgList(*x), data1D );
      delete data1D;
      dhistMap.insert( {Form("Bin_%i",counter),  dataDataHist} );
      if(!isZ)
        w->var(Form("const_a_%i",counter))->setConstant();
      w->var(Form("const_b_%i",counter))->setConstant();
      
      if(!fitData){
        w->var(Form("coeff_sgnl_%i",counter))->setConstant();
        w->var(Form("coeff_2S_%i",counter))->setConstant();
      }
      if(isZ)
        w->var(Form("coeff_2S_%i",counter))->setConstant();

      if(fitBackground)pdfData->fitTo( *dataDataHist, PrintLevel(-1), PrintEvalErrors(-1) );
      
      w->var(Form("const_a_%i",counter))->setConstant();
      w->var(Form("const_b_%i",counter))->setConstant();
      w->var(Form("coeff_sgnl_%i",counter))->setConstant();
      w->var(Form("coeff_2S_%i",counter))->setConstant();
      w->var(Form("sigmaCB_2S_scale_%i",counter))->setConstant();
      
      if(drawPlots){
        auto frameFit = x->frame();
        dataDataHist->plotOn(frameFit);
        pdfData->plotOn(frameFit, RooFit::LineColor(kBlue));
        pdfData->plotOn(frameFit, RooFit::Components(Form("pdf_bkg_%i",counter)) , RooFit::LineColor(kBlue), RooFit::LineStyle(kDashed));
        pdfData->plotOn(frameFit, RooFit::Components(Form("pdfDCB_2S_%i",counter)) , RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));
        cValid->cd(1);
        frameFit_ctl->Draw();
        cValid->cd(2);
        frameFit->Draw();
        cValid->Print(PrintToFile);
      
        delete frameFit;
      }
      delete frameFit_ctl;
      
      ++counter;   
    
  }
  cValid->Print(PrintToFile+"]");

  
  //define RooArgSets for convenience
  w->defineSet("obs","x"); //observables
  w->defineSet("poi","deltaS1,deltadZ1"); //parameters of interest
  w->defineSet("np",nps); //nuisance parameters
  w->Print();
  std::cout << "Number of bins to fit: " << counter << " " << data3D->GetEntries() << std::endl;
  
  
  
  for( auto constants:  constantsUsed){
    for( auto val:  constants){
      std::cout << val << ", \t" ;
    }

    std::cout << std::endl;
  }

  // Associate model with the physics state and model_ctl with the control state
  RooSimultaneous simPdf("simPdf", "simultaneous pdf", pdfMap, categories );
  RooDataHist combData("combData", "combData", *x, categories, dhistMap);
 
  simPdf.fitTo( combData );
  
  //RooAbsReal* nll = simPdf.createNLL(combData);
  
  //TCanvas* c2  =  new TCanvas("c2","c2",800,800);

  
  //RooPlot* frame = w->var("deltaS1")->frame();
  //nll->plotOn(frame,ShiftToZero()); //the ShiftToZero option puts the minimum at 0 on the y-axis
  //frame->Draw();
  //c2->Print("nll.pdf");
  
  return;

  
}


      /*

      
      
      // Build PDF


       
//      hists.push_back( Mass2AADn2->ProjectionZ("dn2", binX+1,binX+1, binY+1,binY+1) );
      hists.push_back( Mass2AADn1->ProjectionZ("dn1", binX+1,binX+1, binY+1,binY+1) );
       //hists.push_back( Mass2AA->ProjectionZ("nb") );
      hists.push_back( Mass2AAUp1->ProjectionZ("up1",  binX+1,binX+1, binY+1,binY+1) );
      hists.push_back( Mass2AAUp2->ProjectionZ("up2",  binX+1,binX+1, binY+1,binY+1) );
      
      buildRooMomentMorph(  hists, biases, x, &deltax, Form("Morph_%i", counter), w) ;
      auto morph =  w->pdf(Form("Morph_Morph_%i", counter));
      pdfMap.insert( {Form("Bin_%i",counter),  morph} );
      for(auto hist: hists)
        delete hist;
      hists.clear();
       
       */
