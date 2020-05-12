#ifndef ROOT_Hfitter_RooMyPoly
#define ROOT_Hfitter_RooMyPoly

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"
#include "RooAbsReal.h"
#include "RooArgSet.h"

#include <vector>

class RooRealVar;


class RooMyPoly : public RooAbsPdf {
public:
   RooMyPoly() ;
   RooMyPoly(const char* name, const char* title, RooAbsReal& x) ;
   RooMyPoly(const char *name, const char *title,
       RooAbsReal& _x, const RooArgList& _coefList, Int_t lowestOrder=1) ;
 
   RooMyPoly(const RooMyPoly& other, const char* name = 0);
   virtual TObject* clone(const char* newname) const { return new RooMyPoly(*this, newname); }
   virtual ~RooMyPoly() ;
 

   Double_t evaluate() const;

   Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
   Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;
 
 protected:
 
   RooRealProxy _x;
   RooListProxy _coefList ;
   Int_t _lowestOrder ;
 
   mutable std::vector<Double_t> _wksp; //! do not persist
 
 

private:

  ClassDef(RooMyPoly, 1)
};


#endif


