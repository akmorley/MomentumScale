
#include "RooMyPoly.h"
#include "RooAbsReal.h"
#include "RooArgList.h"
#include "RooMsgService.h"

#include "TError.h"

#include <cmath>
#include <cassert>
#include <vector>
using namespace std;

ClassImp(RooMyPoly);
 
 ////////////////////////////////////////////////////////////////////////////////
 /// coverity[UNINIT_CTOR]
 
 RooMyPoly::RooMyPoly()
 {
 }
 
 ////////////////////////////////////////////////////////////////////////////////
 /// Create a polynomial in the variable `x`.
 /// \param[in] name Name of the PDF
 /// \param[in] title Title for plotting the PDF
 /// \param[in] x The variable of the polynomial
 /// \param[in] coefList The coefficients \f$ a_i \f$
 /// \param[in] lowestOrder [optional] Truncate the sum such that it skips the lower orders:
 /// \f[
 ///     1. + \sum_{i=0}^{\mathrm{coefList.size()}} a_{i} * x^{(i + \mathrm{lowestOrder})}
 /// \f]
 ///
 /// This means that
 /// \code{.cpp}
 /// RooMyPoly pol("pol", "pol", x, RooArgList(a, b), lowestOrder = 2)
 /// \endcode
 /// computes
 /// \f[
 ///   \mathrm{pol}(x) = 1 * x^0 + (0 * x^{\ldots}) + a * x^2 + b * x^3.
 /// \f]
 
 
 RooMyPoly::RooMyPoly(const char* name, const char* title,
               RooAbsReal& x, const RooArgList& coefList, Int_t lowestOrder) :
   RooAbsPdf(name, title),
   _x("x", "Dependent", this, x),
   _coefList("coefList","List of coefficients",this),
   _lowestOrder(lowestOrder)
 {
   // Check lowest order
   if (_lowestOrder<0) {
     coutE(InputArguments) << "RooMyPoly::ctor(" << GetName()
            << ") WARNING: lowestOrder must be >=0, setting value to 0" << endl ;
     _lowestOrder=0 ;
   }
 
   RooFIter coefIter = coefList.fwdIterator() ;
   RooAbsArg* coef ;
   while((coef = (RooAbsArg*)coefIter.next())) {
     if (!dynamic_cast<RooAbsReal*>(coef)) {
       coutE(InputArguments) << "RooMyPoly::ctor(" << GetName() << ") ERROR: coefficient " << coef->GetName()
              << " is not of type RooAbsReal" << endl ;
       R__ASSERT(0) ;
     }
     _coefList.add(*coef) ;
   }
 }
 
 ////////////////////////////////////////////////////////////////////////////////
 
 RooMyPoly::RooMyPoly(const char* name, const char* title,
                            RooAbsReal& x) :
   RooAbsPdf(name, title),
   _x("x", "Dependent", this, x),
   _coefList("coefList","List of coefficients",this),
   _lowestOrder(1)
 { }
 
 ////////////////////////////////////////////////////////////////////////////////
 /// Copy constructor
 
 RooMyPoly::RooMyPoly(const RooMyPoly& other, const char* name) :
   RooAbsPdf(other, name),
   _x("x", this, other._x),
   _coefList("coefList",this,other._coefList),
   _lowestOrder(other._lowestOrder)
 { }
 
 ////////////////////////////////////////////////////////////////////////////////
 /// Destructor
 
 RooMyPoly::~RooMyPoly()
 { }
 
 ////////////////////////////////////////////////////////////////////////////////
 
 Double_t RooMyPoly::evaluate() const
 {
   // Calculate and return value of polynomial
 
   const unsigned sz = _coefList.getSize();
   const int lowestOrder = _lowestOrder;
   if (!sz) return lowestOrder ? 1. : 0.;
   _wksp.clear();
   _wksp.reserve(sz);
   {
     const RooArgSet* nset = _coefList.nset();
     RooFIter it = _coefList.fwdIterator();
     RooAbsReal* c;
     while ((c = (RooAbsReal*) it.next())) _wksp.push_back(c->getVal(nset));
   }
   const Double_t x = _x;
   Double_t retVal = _wksp[sz - 1];
   for (unsigned i = sz - 1; i--; ) retVal = _wksp[i] + x * retVal;
   retVal =  retVal * std::pow(x, lowestOrder) + (lowestOrder ? 1.0 : 0.0);
   if(retVal<0) retVal = 0;
   return retVal;
 }
 
 
 
 ////////////////////////////////////////////////////////////////////////////////
 /// Advertise to RooFit that this function can be analytically integrated.
 Int_t RooMyPoly::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
 {
   if (matchArgs(allVars, analVars, _x)) return 1;
   return 0;
 }
 
 ////////////////////////////////////////////////////////////////////////////////
 /// Do the analytical integral according to the code that was returned by getAnalyticalIntegral().
 Double_t RooMyPoly::analyticalIntegral(Int_t code, const char* rangeName) const
 {
   R__ASSERT(code==1) ;
 
   const Double_t xmin = _x.min(rangeName), xmax = _x.max(rangeName);
   const int lowestOrder = _lowestOrder;
   const unsigned sz = _coefList.getSize();
   if (!sz) return xmax - xmin;
   _wksp.clear();
   _wksp.reserve(sz);
   {
     const RooArgSet* nset = _coefList.nset();
     RooFIter it = _coefList.fwdIterator();
     unsigned i = 1 + lowestOrder;
     RooAbsReal* c;
     while ((c = (RooAbsReal*) it.next())) {
       _wksp.push_back(c->getVal(nset) / Double_t(i));
       ++i;
     }
   }
   Double_t min = _wksp[sz - 1], max = _wksp[sz - 1];
   for (unsigned i = sz - 1; i--; )
     min = _wksp[i] + xmin * min, max = _wksp[i] + xmax * max;
   return max * std::pow(xmax, 1 + lowestOrder) - min * std::pow(xmin, 1 + lowestOrder) +
       (lowestOrder ? (xmax - xmin) : 0.);
 }