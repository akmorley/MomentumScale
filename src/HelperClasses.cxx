#include "HelperClasses.h"
#include <iostream>
#include "TCanvas.h"

#include "RooAbsPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooFormulaVar.h"

ScaleFit::ScaleFit(  RooWorkspace* w, TString name, ScaleFitData& sfd ):
 m_name(name),
 m_w(w),
 m_sfd(sfd)
{

  m_bhm = BinMergeHelper( m_sfd.nAx, m_sfd.nAy, m_sfd.nM);

  m_bhm.entrySummary = (TH2*) m_sfd.mc3D->Project3D("yx");
  m_bhm.entrySummary->SetName( "entrySummary"+name );

  int binXmax = m_bhm.entrySummary->GetNbinsX();
  int binYmax = m_bhm.entrySummary->GetNbinsY();
  for( int binX(0); binX < binXmax; ++binX ){
    for( int binY(0); binY < binYmax; ++binY ){
      double  val = m_bhm.entrySummary->GetBinContent( binX+1, binY+1 );
      //std::cout << "Checking bin: " <<  binX << ", " <<  binY << " : " << val  << std::endl;
      if(val < sfd.considerThreshold){
        continue;
      }
      //std::cout << "Inserting bin: " <<  binX << ", " <<  binY << " : " << val  << std::endl;
      m_bhm.bins.insert( { {binX,binY}, BinContent( binX, binY, val) }  );
    }
  }

  m_bhm.mergeBins(sfd.entryThreshold);
  m_bhm.printBins();


}

void ScaleFit::constructModel( RooRealVar* x, RooArgList& list, SimFitConfig& sfc )
{

  std::vector<std::vector<double>> constantsUsed;
  RooArgSet nps;

  if(!m_sfd.isZ)
    m_w->factory("delta2S[3.9965214834]");      //3.686097^2 - 3.096900^2  = 3.9965214834
 

  TCanvas* cValid(0); 
  TString PrintToFile("TestFits_"+m_name+".pdf");
  if(sfc.doPrint){
    cValid =  new TCanvas("cValid","cValid",1600,800);
    cValid->Divide(2,1);
    cValid->Print(PrintToFile+"[");
  }
  
  int counter(0);
  for( auto&  binPair : m_bhm.bins ){
    const TString binName = Form( "%s_%i", m_name.Data(), counter);
    auto& bin = binPair.second;
    int binX = bin.binX;
    int binY = bin.binY;
    
    std::cout << binX << " "<< binY << " Bin "<< counter << " weight : " << bin.totalWeight <<  std::endl;

    
    // Setup linear equations relating the change in mass to the change in scale
    // deltaM/M = as_1 * deltaS_1 + as_2 * deltaS_2+ az_1 * deltadZ_1 + az_2 * deltadZ_2
    auto as_1   = (RooRealVar*) m_w->factory(Form("as_1_%s[%f]",binName.Data(),2.*bin.getValue(m_sfd.As1_2D))); 
    auto as_2   = (RooRealVar*) m_w->factory(Form("as_2_%s[%f]",binName.Data(),2.*bin.getValue(m_sfd.As2_2D)));
    auto adz_1  = (RooRealVar*) m_w->factory(Form("adz_1_%s[%f]",binName.Data(),2.*bin.getValue(m_sfd.Az1_2D))); 
    auto adz_2  = (RooRealVar*) m_w->factory(Form("adz_2_%s[%f]",binName.Data(),2.*bin.getValue(m_sfd.Az2_2D))); 
    auto deltax = (RooFormulaVar*) m_w->factory( Form("expr::deltax_%s('@0*@4+@1*@5+@2*@6+@3*@7',{as_1_%s,as_2_%s,adz_1_%s,adz_2_%s,%s,%s,%s,%s} )", binName.Data(),binName.Data(),binName.Data(),binName.Data(),binName.Data(),list[0].GetName(),list[1].GetName(),list[2].GetName(),list[3].GetName()));
    
    std::cout << "Built Linear Equation: " << binName.Data() << std::endl;
    // Dummy test fit 
    //auto deltax = (RooAbsReal*) m_w->factory( Form("deltax_%i[0,-0.05,0.05]", counter));
    
    //Create PDFs for the signal shape and the background
    if(m_sfd.isZ){
      std::cout << Form("RooTwoSidedCBShape::pdfDCB_ctl_%s(%s, muCB_ctl_%s[8281,7500,9000], sigmaCB_ctl_%s[1000,100,5000], alphaCBLo_%s[1,0.1,10], nCBLo_%s[2,0.01,100], alphaCBHi_%s[1,0.1,10], nCBHi_%s[2,0.01,100])", binName.Data(), x->GetName(),binName.Data(), binName.Data(), binName.Data(), binName.Data(), binName.Data(), binName.Data() ) << std::endl;
      m_w->factory(Form("RooTwoSidedCBShape::pdfDCB_ctl_%s(%s, muCB_ctl_%s[8281,7500,9000], sigmaCB_ctl_%s[1000,100,5000], alphaCBLo_%s[1,0.1,10], nCBLo_%s[2,0.01,100], alphaCBHi_%s[1,0.1,10], nCBHi_%s[2,0.01,100])", binName.Data(), x->GetName(),binName.Data(), binName.Data(), binName.Data(), binName.Data(), binName.Data(), binName.Data() ) );
      std::cout << "Built pdfDCB_ctl_ " << binName.Data() << std::endl;
      std::cout << Form("RooTwoSidedCBShape::pdfDCB_%s(%s, expr::muCB_%s('@0*(1.+@1)', {muCB_ctl_%s,deltax_%s} ), prod::sigmaCB_%s(sigmaCB_ctl_%s,sigmaCB_scale_%s[1.,0.5,5]), alphaCBLo_%s, nCBLo_%s, alphaCBHi_%s, nCBHi_%s)", binName.Data(), x->GetName(), binName.Data(), binName.Data(),binName.Data(),binName.Data(),binName.Data(), binName.Data(), binName.Data(),binName.Data(),binName.Data(),binName.Data()  )  << std::endl;
      m_w->factory(Form("RooTwoSidedCBShape::pdfDCB_%s(%s, expr::muCB_%s('@0*(1.+@1)', {muCB_ctl_%s,deltax_%s} ), prod::sigmaCB_%s(sigmaCB_ctl_%s,sigmaCB_scale_%s[1.,0.5,5]), alphaCBLo_%s, nCBLo_%s, alphaCBHi_%s, nCBHi_%s)", binName.Data(), x->GetName(), binName.Data(), binName.Data(),binName.Data(),binName.Data(),binName.Data(), binName.Data(), binName.Data(),binName.Data(),binName.Data(),binName.Data()  ) );    
      std::cout << "Built pdfDCB_" << std::endl;
      m_w->factory(Form("RooPolynomial::pdf_bkg_%s(%s, {const_a_%s[0,-1000,1000],const_b_%s[0,-1000,1000]} )", binName.Data(), x->GetName(), binName.Data(), binName.Data() ) );
      std::cout << "Built pdf_bkg_" << std::endl;
      m_w->factory(Form("SUM::pdf_data_%s(coeff_sgnl_%s[1,0.1,1]*pdfDCB_%s, pdf_bkg_%s)",binName.Data(),binName.Data(),binName.Data(),binName.Data()));
      std::cout << "Built pdf_data_" << std::endl;
    }else{
      m_w->factory(Form("RooTwoSidedCBShape::pdfDCB_ctl_%s(%s, muCB_ctl_%s[9.6,8,10], sigmaCB_ctl_%s[1,0.01,4], alphaCBLo_%s[1,0.1,10], nCBLo_%s[2,0.01,100], alphaCBHi_%s[1,0.1,10], nCBHi_%s[2,0.01,100])", binName.Data(), x->GetName(),binName.Data(), binName.Data(),binName.Data(),binName.Data(),binName.Data(), binName.Data() ) );
      m_w->factory(Form("RooTwoSidedCBShape::pdfDCB_%s(%s, expr::muCB_%s('@0*(1.+@1)',  {muCB_ctl_%s,deltax_%s}), prod::sigmaCB_%s(sigmaCB_ctl_%s,sigmaCB_scale_%s[1.,0.5,5]), alphaCBLo_%s, nCBLo_%s, alphaCBHi_%s, nCBHi_%s)", binName.Data(), x->GetName(), binName.Data(), binName.Data(),binName.Data(),binName.Data(),binName.Data(), binName.Data(), binName.Data(),binName.Data(),binName.Data(),binName.Data()  ) );    
      m_w->factory(Form("RooTwoSidedCBShape::pdfDCB_2S_%s(%s, expr::muCB_2S_%s('@0+@1', {muCB_ctl_%s,delta2S}), prod::sigmaCB_2S_%s(sigmaCB_ctl_%s,sigmaCB_2S_scale_%s[1.,0.5,5]), alphaCBLo_%s, nCBLo_%s, alphaCBHi_%s, nCBHi_%s)", binName.Data(), x->GetName(), binName.Data(), binName.Data(), binName.Data(), binName.Data(), binName.Data(), binName.Data(), binName.Data(), binName.Data(), binName.Data() ) );
      m_w->factory(Form("RooPolynomial::pdf_bkg_%s(%s, {const_a_%s[0,-1000,1000],const_b_%s[0,-1000,1000]} )", binName.Data(), x->GetName(), binName.Data(), binName.Data() ) );
      m_w->factory(Form("SUM::pdf_data_%s(coeff_sgnl_%s[1,0.1,1]*pdfDCB_%s, coeff_2S_%s[0.0,0,0.2]*pdfDCB_2S_%s, pdf_bkg_%s)",binName.Data(),binName.Data(),binName.Data(),binName.Data(),binName.Data(),binName.Data()));
    }

    std::cout << "Built PDFs" << std::endl;


    auto pdfDCB_ctl = m_w->pdf(Form("pdfDCB_ctl_%s", binName.Data()));
    auto pdfData = m_w->pdf(Form("pdf_data_%s", binName.Data()));
    
    // Get Histogram for MC -- this will be the control sample
    TH1* mc1D = bin.getProjection( m_sfd.mc3D );
    mc1D->SetName("MCProj");
    RooDataHist* dataDataHist_ctl = new RooDataHist(Form("DataHist_ctl_%s",binName.Data()), Form("DataHist_%s",binName.Data()), RooArgList(*x), mc1D );
    delete mc1D;
    
    // Fit the control sample if the fit fail's dont use the bin
    auto fitOk = pdfDCB_ctl->fitTo( *dataDataHist_ctl, RooFit::PrintLevel(-1), RooFit::Save() );
    if( !fitOk || fitOk->status() > 1 )
    {  
      delete dataDataHist_ctl;
      continue;
    }
    m_w->import(*dataDataHist_ctl);

    // Define category name to be use in the simulatenous fit
    sfc.categories.defineType(Form("Bin_%s",binName.Data()));
    sfc.categories.defineType(Form("Bin_ctl_%s",binName.Data()));
   

    // Store the pdf and data hist in maps for simulatenous fit
    sfc.pdfMap.insert( {Form("Bin_ctl_%s",binName.Data()),  pdfDCB_ctl} );
    sfc.pdfMap.insert( {Form("Bin_%s",binName.Data()),  pdfData} );
    sfc.dhistMap.insert( {Form("Bin_ctl_%s",binName.Data()),  dataDataHist_ctl} );

    // Fix shape parameters to reduce the DOF
    m_w->var(Form("alphaCBLo_%s",binName.Data()))->setConstant();
    m_w->var(Form("nCBLo_%s",binName.Data()))->setConstant();
    m_w->var(Form("alphaCBHi_%s",binName.Data()))->setConstant();
    m_w->var(Form("nCBHi_%s",binName.Data()))->setConstant();
    m_w->var(Form("sigmaCB_ctl_%s",binName.Data()))->setConstant();
    
    // Ideally these two paramters would be left to float
    m_w->var(Form("muCB_ctl_%s",binName.Data()))->setConstant();
    m_w->var(Form("sigmaCB_scale_%s",binName.Data()))->setConstant();
    nps.add(  *m_w->var(Form("muCB_ctl_%s",binName.Data())) );
    nps.add(  *m_w->var(Form("sigmaCB_ctl_%s",binName.Data())) );
  
    if( sfc.doFloatNP )
    {
      m_w->var(Form("muCB_ctl_%s",binName.Data()))->setConstant(false);
      m_w->var(Form("sigmaCB_scale_%s",binName.Data()))->setConstant(false);
    }
    
    
    // Do some plotting
    auto frameFit_ctl  = x->frame();
    if(sfc.doPrint){
      dataDataHist_ctl->plotOn(frameFit_ctl);
      pdfDCB_ctl->plotOn(frameFit_ctl, RooFit::LineColor(kBlue));
    }


    // Get Histogram for Data -- this will be used to measure the bias
    TH1* data1D = bin.getProjection( m_sfd.data3D );
    data1D->SetName("Data");
    RooDataHist* dataDataHist = new RooDataHist(Form("DataHist_%s",binName.Data()), Form("DataHist_%s",binName.Data()), RooArgList(*x), data1D );
    delete data1D;
    sfc.dhistMap.insert( {Form("Bin_%s",binName.Data()),  dataDataHist} );
    m_w->import(*dataDataHist);

    // Reduce the NDOF for the background shape
    if(!m_sfd.isZ)
      m_w->var(Form("const_a_%s",binName.Data()))->setConstant();
    m_w->var(Form("const_b_%s",binName.Data()))->setConstant();
    
    // If an MC only test with not background do not fit background fractions
    if(!sfc.fitBackground){
      m_w->var(Form("coeff_sgnl_%s",binName.Data()))->setConstant();
      m_w->var(Form("coeff_2S_%s",binName.Data()))->setConstant();
    }
    
    if(sfc.fitData)
      pdfData->fitTo( *dataDataHist, RooFit::PrintLevel(-1), RooFit::PrintEvalErrors(-1) );
    
    //Store values of the linear equations for quick reference
    constantsUsed.push_back( {as_1->getVal(),as_2->getVal(),adz_1->getVal(),adz_2->getVal(), deltax->evaluate() });
    
    //Set all background shapes == constant again to reduce the NDOF
    m_w->var(Form("const_a_%s",binName.Data()))->setConstant();
    m_w->var(Form("const_b_%s",binName.Data()))->setConstant();
    m_w->var(Form("coeff_sgnl_%s",binName.Data()))->setConstant();
    if(!m_sfd.isZ){
      m_w->var(Form("coeff_2S_%s",binName.Data()))->setConstant();
      m_w->var(Form("sigmaCB_2S_scale_%s",binName.Data()))->setConstant();
    }

    if(sfc.doPrint){
      auto frameFit = x->frame();
      dataDataHist->plotOn(frameFit);
      pdfData->plotOn(frameFit, RooFit::LineColor(kBlue));
      pdfData->plotOn(frameFit, RooFit::Components(Form("pdf_bkg_%s",binName.Data())) , RooFit::LineColor(kBlue), RooFit::LineStyle(kDashed));
      if(!m_sfd.isZ)
        pdfData->plotOn(frameFit, RooFit::Components(Form("pdfDCB_2S_%s",binName.Data())) , RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));
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
  if(sfc.doPrint)
    cValid->Print(PrintToFile+"]");

  delete cValid;

  //define RooArgSets for convenience
  m_w->defineSet("nuisParams",nps); //nuisance parameters
  m_w->Print();
  
  std::cout << "Number of bins to fit: " << counter << std::endl;
  for( auto constants:  constantsUsed){
    for( auto val:  constants){
      std::cout << val << ", \t" ;
    }

    std::cout << std::endl;
  }



}

BinContent::BinContent():
  binX( -1 ),
  binY( -1 ),
  weight( 0 ),
  totalWeight( 0 )
  {}

BinContent::BinContent(int x, int y, double w ):
  binX( x ),
  binY( y ),
  weight( w ),
  totalWeight( w )
  {}

BinContent::BinContent( const BinContent& bin ):
binX( bin.binX ),
binY( bin.binY ),
weight( bin.weight ),
totalWeight( bin.totalWeight )
{
  binsToMerge = bin.binsToMerge;
}


void BinContent::addBinContent( BinContent& binContent)
{
  //while (!binContent.binsToMerge.empty())
  //  binsToMerge.merge(binContent.binsToMerge.extract(binContent.binsToMerge.begin()));
  binsToMerge.merge(binContent.binsToMerge);
  binsToMerge.insert( binContent );
  totalWeight += binContent.totalWeight;
}

bool BinContent::operator<(BinContent other) const
{
  if ( binX > other.binX )
    return true;
  else if ( binX < other.binX )
    return false;

  if ( binY >= other.binY )
    return true;
  else if ( binY < other.binY )
    return false;  
  return true;  
}

void BinContent::printBin() const{
  std::cout << Form("T**** %3i, %3i  %9.0f ", binX, binY, totalWeight ) << std::endl; 
}
  
void BinContent::printBinsToMerge() const{
  auto theBin = binsToMerge.begin();
  while(  theBin != binsToMerge.end() ){
    std::cout << Form(" |*** %3i, %3i  %9.0f ", theBin->binX, theBin->binY, theBin->weight ) << std::endl; 
    ++theBin;
  }
}

TH1* BinContent::getProjection(const TH3* hist3D){
  auto hist = hist3D->ProjectionZ("AProjection", binX+1, binX+1, binY+1, binY+1 );
  //std::cout << " projection " << hist->GetEntries() << std::endl;
  for( auto&  bin : binsToMerge ){
    auto hist2 = hist3D->ProjectionZ("Projection2", bin.binX+1, bin.binX+1, bin.binY+1, bin.binY+1 );
    hist->Add(hist2);
    delete hist2;
    //std::cout << " ------------- " << hist->GetEntries() << std::endl;
  }
  return hist;
}

double BinContent::getValue(const TH2* hist2D){
  double value = hist2D->GetBinContent( binX+1, binY+1 ) * weight;
  for( auto&  bin : binsToMerge ){
    value += hist2D->GetBinContent( bin.binX+1, bin.binY+1 ) * bin.weight;
  }
  return value/totalWeight;
}



  
BinMergeHelper::BinMergeHelper( int ax, int ay, int m):
axNbins(ax),
ayNbins(ay),
mNbins(m),
entrySummary(nullptr){};


void BinMergeHelper::mergeHistBins( TH2* hist, bool add ){
  for (auto it = bins.cbegin(); it != bins.cend(); ++it)
  {
    double weight = add ? 1 : it->second.weight;
    double sumTotal  = hist->GetBinContent( it->second.binX+1,it->second.binY+1 ) * weight;
    double sumWeight = weight;
    if(verbose) std::cout << Form("M**** %3i, %3i  %9.0f %9.4f %9.4f ", it->second.binX, it->second.binY, it->second.weight,sumTotal/sumWeight, hist->GetBinContent( it->second.binX+1,it->second.binY+1 ) ) << std::endl;
    for (auto it2 = it->second.binsToMerge.cbegin(); it2 !=  it->second.binsToMerge.cend(); ++it2)
    {
      weight = add ? 1 : it2->weight;
      sumTotal  += hist->GetBinContent( it2->binX+1,it2->binY+1 ) * weight;
      sumWeight += add ? 0. : it2->weight;
      if(verbose) std::cout << Form(" |*** %3i, %3i  %9.0f  %9.4f %9.4f ", it2->binX, it2->binY, it2->weight,sumTotal/sumWeight,  hist->GetBinContent( it2->binX+1,it2->binY+1 ) ) << std::endl;
      hist->SetBinContent( it2->binX+1,it2->binY+1,0 ); 
      hist->SetBinError( it2->binX+1,it2->binY+1,0 ); 
      
    }
    hist->SetBinContent( it->second.binX+1,it->second.binY+1, sumTotal/sumWeight );
  }
}

void BinMergeHelper::printBins() const{
  auto theBin = bins.begin();
  while(  theBin != bins.end() ){
    theBin->second.printBin();
    theBin->second.printBinsToMerge();
    ++theBin;
  }
}

void BinMergeHelper::mergeBins(double threshold){
  //Search procedes in this order merge 
  // Bins with the same Az1,Az2 as As's have little impact on the final result
  for( int mM(0); mM < mNbins; ++mM ){  
    // For each az1 az2 pair
    for( int binXa(0); binXa < axNbins; ++binXa ){
      for( int binYa(0); binYa < ayNbins; ++binYa ){
        // Merge X,Ym
        int x(0);
        while( x<mM ){
          auto theBin = bins.find( { x*axNbins + binXa, mM*ayNbins + binYa} );
          if(theBin!=bins.end()){
            mergeBin(theBin,threshold);
            printBins();
          }              
          ++x;
        }
        // Merge Xm,Y
        int y(0);
        while( y<mM ){
          auto theBin = bins.find( { mM*axNbins + binXa, y*ayNbins + binYa} );
          if(theBin!=bins.end())
            mergeBin(theBin,threshold);              
          ++y;
        }
        // Merge  Xm=Ym 
        auto theBin = bins.find( { mM*axNbins + binXa, mM*ayNbins + binYa} );
        if(theBin!=bins.end())
          mergeBin(theBin,threshold);              
      }
    }
  }
  if( entrySummary ){
   
    int binXmax = entrySummary->GetNbinsX();
    int binYmax = entrySummary->GetNbinsY();
    for( int binX(0); binX < binXmax; ++binX ){
      for( int binY(0); binY < binYmax; ++binY ){
        entrySummary->SetBinContent( binX+1, binY+1, 0 );
      }
    }
    std::multimap<double, BinContent> belowThreshold;
    for (auto it = bins.cbegin(); it != bins.cend() ; ++it )
    {
      entrySummary->SetBinContent( it->second.binX+1, it->second.binY+1 , it->second.totalWeight );
      if( it->second.totalWeight  < threshold ) belowThreshold.insert( {-it->second.totalWeight,it->second} );
    }



    for( int mM(0); mM < mNbins; ++mM ){  
      // For each az1 az2 pair
      std::cout << "********mM " << mM << std::endl;
      for( int binXa(0); binXa < axNbins; ++binXa ){
        for( int binYa(0); binYa < ayNbins; ++binYa ){
          if(entrySummary->GetBinContent( mM*axNbins + binXa, mM*ayNbins + binYa)>=threshold)
          {
            std::cout <<  "\x1B[34;42m" << Form(" %15.4f", entrySummary->GetBinContent( mM*axNbins + binXa, mM*ayNbins + binYa)) <<  "\x1B[0m" ;
          } else {
            std::cout << Form(" %15.4f", entrySummary->GetBinContent( mM*axNbins + binXa, mM*ayNbins + binYa))  ;          
          }
         }
        std::cout << std::endl;
      }
    }


/*    auto itA = belowThreshold.begin();
    while( itA != belowThreshold.end() ){
      //auto addToBin = bins.find( {itA->second.binX,itA->second.binY} );
      //while
      std::cout << itA->first << std::endl;
      ++itA;
    }
    */

  }


  // Final bin clean up
  for (auto it = bins.cbegin(); it != bins.cend() ; )
  {
    if ( it->second.totalWeight < threshold )
    {
      it = bins.erase(it);
    } else {
      ++it;
    }
  } 
}


bool BinMergeHelper::mergeBinToIndex( BinIter theBin, int binX, int binY, bool deleteBin  ){
  auto addToBin = bins.find( {binX,binY} );
  if( addToBin == bins.end() ){
    if(verbose) std::cout << " Not Found! " << std::endl;  
    if(deleteBin) bins.erase(theBin);
    return false;
  }else{
    if(verbose) std::cout << " Found! " << std::endl;  
    addToBin->second.addBinContent(theBin->second);
    bins.erase(theBin);
    return true;
  }  
}

void BinMergeHelper::mergeBin(BinIter theBin, double threshold){
  int binXa = theBin->second.binX % axNbins;
  int binYa = theBin->second.binY % ayNbins;
  int binXm = theBin->second.binX / axNbins;
  int binYm = theBin->second.binY / ayNbins;
  if(verbose) 
    std::cout << Form("MB*** (%3i, %3i) (%3i, %3i) (%3i, %3i) %9.0f  ",  
                      theBin->second.binX, theBin->second.binY, 
                      binXa, binYa, binXm, binYm, theBin->second.totalWeight) << std::endl;

  if(  theBin->second.totalWeight < threshold  ){
    if( binXm < binYm  ){
      if(verbose) std::cout << " |*** Checking X+1:  ";  
      mergeBinToIndex( theBin, theBin->second.binX + axNbins, theBin->second.binY );
    } else if ( binXm > binYm ) {
      if(verbose) std::cout << " |*** Checking Y+1:  ";  
      mergeBinToIndex( theBin, theBin->second.binX, theBin->second.binY + ayNbins);
    } else {
      if( binXm < mNbins-1 ){
        if(verbose) std::cout << " |*** Checking X+1, Y+1:  ";  
        if( mergeBinToIndex( theBin, theBin->second.binX + axNbins, theBin->second.binY + ayNbins) )
          return;
      } 
      if( binXm > 0 ){
        if(verbose) std::cout << " |*** Checking X-1, Y-1:  ";  
        mergeBinToIndex( theBin, theBin->second.binX - axNbins, theBin->second.binY - ayNbins ); 
      }
    }
  }
}   

