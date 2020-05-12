#ifndef HELPER_CLASSES_H
#define HELPER_CLASSES_H

#include "TH3.h"
#include "TH2.h"
#include "TH1.h"
#include <vector>
#include <map>
#include <set>

#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooCategory.h"



struct BinContent{
  BinContent();

  BinContent( int x, int y, double w);

  BinContent( const BinContent& bin );

  // Add a bin to list of bins
  void addBinContent( BinContent& binContent);

  // Comparator for ordering the bins
  bool operator<(BinContent other) const;

  // Dump bin content
  void printBin() const;

  // Dump bins that have been merged in this bin
  void printBinsToMerge() const;

  // Get 1D projection from a 3D hist of all bins associated to this bin
  TH1* getProjection(const TH3* hist3D);

  // Get weighted average of all bins associated to this bin
  double getValue(const TH2* hist2D);

  // Bin x value (offset -1 wrt to root)
  int binX;
  // Bin y value (offset -1 wrt to root)
  int binY;
  // Weight of current bin
  double weight; 
  // Weight of all bins associated to this bin
  double totalWeight; 
  // All of the bins associated 
  std::set<BinContent> binsToMerge;
};

struct BinMergeHelper{

  //Constructor
  BinMergeHelper(){};
  BinMergeHelper( int ax, int ay, int m);

  int axNbins;
  int ayNbins;
  int mNbins;
  bool verbose = true;
  std::map< std::vector<int> , BinContent> bins;
  TH2* entrySummary;

  typedef std::map< std::vector<int> , BinContent>::iterator BinIter;


  // Merge bins in  a histogram according to what is 
  // choice between direct addition or weighted average
  void mergeHistBins( TH2* hist, bool add = true );

  //Print all active bins and their content;
  void printBins() const;

  // Merge all bins. If bin is still below threshold after all possible merging routes delete it. 
  void mergeBins(double threshold);

private:
  // Merge theBin with a bin at index binX binY and optionally delete if it fails to merge
  bool mergeBinToIndex( BinIter theBin, int binX, int binY, bool deleteBin = false  );


  // Merge theBin if is below threshold 
  void mergeBin(BinIter theBin, double threshold);
};

struct ABins{
  ABins( double Ar_t, double Az_t):
  Ar( Ar_t ),
  Az( Az_t )
  {};

  double Ar;
  double Az;

};


struct VariableAxis{
  VariableAxis( double xMin_t, double xMax_t, int nBins_t):
  xMin( xMin_t),
  xMax( xMax_t),
  nBins( nBins_t)
  {};
  
  double xMin;
  double xMax;
  int nBins;

  double range() const{
    return xMax-xMin;
  }
  
  double binCenter(int bin){
    return xMin + ( double(bin) + 0.5 )*this->range() / double( nBins );
  }
};

struct ScaleFitData{
  ScaleFitData(  TH3* data, TH3* mc, TH2* As1, TH2* As2, TH2* Az1, TH2* Az2, int mAx, int mAy, int mM, bool Z=true):
   data3D(data),
   mc3D(mc),
   As1_2D(As1),
   As2_2D(As2),
   Az1_2D(Az1),
   Az2_2D(Az2),
   nAx(mAx),
   nAy(mAy),
   nM(mM),
   isZ(Z),
   entryThreshold(100000),
   considerThreshold(30000)
  {};
  TH3* data3D;
  TH3* mc3D;
  TH2* As1_2D;
  TH2* As2_2D;
  TH2* Az1_2D;
  TH2* Az2_2D;
  int nAx; 
  int nAy; 
  int nM;
  bool isZ;
  double entryThreshold;
  double considerThreshold;
  
};

struct SimFitConfig{
  SimFitConfig():
   categories("AllCategories","AllCategories"),
   fitBackground(false),
   fitData(false),
   doPrint(false),
   doFloatNP(false){};

  SimFitConfig(  bool fitBkg, bool fitD, bool doP, bool adoFloatNP):
   categories("AllCategories","AllCategories"),
   fitBackground(fitBkg),
   fitData(fitD),
   doPrint(doP),
   doFloatNP(adoFloatNP)
  {};



  RooCategory categories;//("AllCategories","AllCategories");
  std::map< std::string, RooDataHist *>   dhistMap;
  std::map< std::string, RooAbsPdf *>     pdfMap;
  bool fitBackground;
  bool fitData; 
  bool doPrint;
  bool doFloatNP;
};

struct ScaleFit{

  ScaleFit(  RooWorkspace* w, TString name, ScaleFitData& data);
  void constructModel( RooRealVar* x,RooArgList& list, SimFitConfig& sfc);

 private:
  TString        m_name;
  RooWorkspace*  m_w;
  ScaleFitData   m_sfd;
  BinMergeHelper m_bhm;
};


#endif
