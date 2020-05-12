#ifndef HELPER_FUNCTIONS_H
#define HELPER_FUNCTIONS_H

#include "TH3.h"
#include "TH2.h"
#include "TH1.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "TF1.h"
#include "TSpectrum.h"

#include "HelperClasses.h"

using namespace RooFit;

TH1* smoothPloyND(TH1* data, double sigma, int order)
{
  if(order < 0){
    return 0;
  }

  TF1 gaus("mygaus","gaus",-200,200);
  gaus.SetNpx(10000);
  gaus.SetParameters(1.,0.,sigma);

  TH1D* clone =  (TH1D*) data->Clone(Form("%s_%s_%dD",data->GetName(),"Smoothed",order));
  const int ndata =  clone->GetNbinsX();


  if( ndata < 2 )
    return clone;

  std::vector<double> x(ndata,0);
  std::vector<double> y(ndata,0);
  std::vector<double> sig(ndata,0);
  std::vector<double> wt2(ndata,0);

  for( int j(0); j < ndata; ++j){
    int j1 = j+1;
    x[j] = clone->GetBinCenter(j1);
    y[j] = clone->GetBinContent(j1);
    if(y[j]!=0){
      sig[j] = clone->GetBinError(j1);
      wt2[j] = pow(sig[j],-2);
    } else {
      sig[j] = 1.;
      wt2[j] = 1.;
    }
  }

  for( int j(0); j < ndata; ++j){
    TMatrixD A(ndata,order+1); // matrix Npar by X
    TMatrixDSym W(ndata); // weight
    TVectorD m(ndata); //  measuements

    for (int i(0);i < ndata; ++i) { // Accumulate sums ...
      //Weight normal weight multiplied by a gaussian weight
      double dx = x[i] - x[j];
      double gausval = gaus.Eval( dx );
      A(i,0) =  1;
      for(int k(1);  k < order+1; ++k){
        A(i,k) = A(i,k-1)*dx;
      }

      W(i,i) = gausval * wt2[i]; //...with weights
      //std::cout << "i " << i << wt2[i]<< " " << W(i,i) <<  std::endl;
      m(i)   = y[i];
    }

    TMatrixD  ATW = A.T() * W;
    TMatrixD C =  ATW * A.T();
    C.Invert();

    //Then the result;
    auto results =  C * ATW * m;

    if(results(0)<0){
      //cout << "Warning Smoothing issue a<0: "<< results(0) << endl;
      results(0)=1.e-6;
    }

    clone->SetBinContent(j+1,results(0));
    clone->SetBinError(j+1,sqrt(C(0,0)));
  }

  return clone;
}


RooHistPdf makeRooHistPdf( const TH1* hist, const TString name, RooRealVar* x )
{
  RooDataHist roohist("dhist_"+name, "dhist_"+name, RooArgList(*x), hist );
  return RooHistPdf("pdfhist_"+name, "pdfhist_"+name, RooArgSet(*x), roohist );
}

void checkHist( TH1* hist){
  for( int i(0); i <= hist->GetNbinsX(); ++i){
    if(hist->GetBinContent(i)<0){
      std::cout << "Bin is invalid:  " << i << " " << hist->GetBinContent(i) << std::endl;
    }
  }
}

TH1* smoothHist( TH1* hist ){
  TSpectrum s;

  TH1* tempData = (TH1*)hist->Clone("temp1D");
  int nbins = hist->GetNbinsX();
  std::vector<double> source(nbins,0);
  for (int j = 0; j < nbins; j++)
    source[j]= hist->GetBinContent(j + 1);
  s.SmoothMarkov(&source[0],nbins,3);
  for (int j = 0; j < nbins; j++){
    tempData->SetBinContent(j + 1,source[j]);
  }
  return tempData;
}


void buildRooMomentMorph( const std::vector<TH1*> hists, const std::vector<double> values, RooAbsReal* x, RooAbsReal* deltax, TString name, RooWorkspace* w )
{
  if( hists.size() != values.size() )
    return;
  RooArgList datahists;
  RooArgList pdfs;
  TVectorD params(hists.size());
  
  for( unsigned int i(0); i < hists.size(); ++i){
    checkHist( hists[i] );
    auto tempData = smoothPloyND( hists[i], 0.1, 2 );//smoothHist( hists[i] );
    checkHist( tempData );

    TString compName = Form( "%s_%i", name.Data(), i) ;
    RooDataHist* roohist = new RooDataHist("dhist_"+compName, "dhist_"+compName, RooArgList(*x), tempData );
    RooHistPdf* roopdf   = new RooHistPdf("pdfhist_"+compName, "pdfhist_"+compName, RooArgSet(*x), *roohist );
    pdfs.add(*roopdf);
    params[i] = values[i];
    delete tempData;
  }

  auto morph = new RooMomentMorph("Morph_"+name,"Morph_"+name,*deltax, RooArgList(*x), pdfs, params, RooMomentMorph::Linear);
  w->import(*morph);   
}



//binX = calcBin( axisA, Ar1 )*axisMA.nBins + calcBin( axisMA, Mrz1 )
//binY = calcBin( axisDA, Ar1+Ar2 )*axisMA.nBins  + calcBin( axisMA, Mrz2 )

ABins binXtoArAz( VariableAxis& axisR, VariableAxis& axisZ, int binX ){
  int binR = binX; // axisZ.nBins
  int binZ = binX % axisZ.nBins;
  double Ar = axisR.binCenter(binR);
  double Mz = axisR.binCenter(binZ);
  double Az = 0.5*Mz-Ar;
  return ABins(Ar,Az);
}

ABins binYtoArAz( VariableAxis& axisR, VariableAxis& axisZ, int binX, double Ar1 ){
  int binR = binX; // axisZ.nBins
  int binZ = binX % axisZ.nBins;
  double dAr = axisR.binCenter(binR);
  double Mz = axisZ.binCenter(binZ);
  double Ar = dAr-Ar1;
  double Az = 0.5*Mz-Ar;
  return ABins(Ar,Az);
}


#endif