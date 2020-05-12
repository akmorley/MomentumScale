#!/usr/bin/env python
# coding: utf-8

# In[1]:


import ROOT
import math
from ROOT import TLorentzVector
from BiasP import *


# In[2]:



class VariableAxis:
    '''Class for axis variables'''
    def __init__(self, xMin, xMax, nBins):
        self.xMin = xMin
        self.xMax = xMax
        self.nBins = nBins
    
    def getBin( self, x ):
        value = (x-self.xMin)/self.binSize()
        abin = math.ceil( value )
        if abin < 0 :
          abin = 0
        elif  abin > self.nBins :
          abin = self.nBins + 1
        return abin 
        
    def binSize(self):
        return self.range() / ( self.nBins )

    def range(self):
        return self.xMax - self.xMin

class MultiDimensionObject:
    def __init__(self):
        self.elements = {}
    
    def addValue(self, key, value):
        self.elements[key] = value
    
    def readValue(self, key):
        try:
          value = self.elements[key]
        except KeyError:
          print( "Failed to get: ", key )
          # could also be 0.0 if using floats...
          value = None
        #print( key, "-->",  value.GetName() )
        return value
    def printObjects(self):
        for key,val in self.elements.items():
          print ( key, "=>", val.GetName() )


 
def calcBin( axis, valA ):
    scaledAF = (valA - axis.xMin)/axis.range()*axis.nBins
    scaledAI = math.floor(scaledAF)
    if scaledAI < 0: 
        scaledAI = 0
    elif scaledAI >= axis.nBins:
        scaledAI = axis.nBins-1
    return scaledAI

def binXtoArAz( axisR, axisZ, binX ):
    binR = binX // axisZ.nBins
    binZ = binX % axisZ.nBins
    Ar = axisR.xMin+(float(binR) + 0.5) * axisR.range() / float(axisR.nBins) 
    Mz = axisZ.xMin+(float(binZ) + 0.5) * axisZ.range() / float(axisZ.nBins)
    Az = 0.5*Mz-Ar
    return (Ar,Az)

def binYtoArAz( axisR, axisZ, binX, Ar1 ):
    binR = binX // axisZ.nBins
    binZ = binX % axisZ.nBins
    dAr = axisR.xMin+(float(binR) + 0.5) * axisR.range() / float(axisR.nBins) 
    Mz = axisZ.xMin+(float(binZ) + 0.5) * axisZ.range() / float(axisZ.nBins)
    Ar = dAr-Ar1
    Az = 0.5*Mz-Ar
    return (Ar,Az)

def main():
 #MyTreeMC=ROOT.TChain("ID_InDetTrackParticle")
#MyTreeMC.Add("/eos/user/m/martis/data/ZmumuNtuples/Jpsi_MCAll.root")
#MyTreeMC.Add("/eos/user/m/martis/data/ZmumuNtuples/Run2Data17_JPsi.root")
#MyTreeMC.Add("/eos/user/m/martis/data/ZmumuNtuples/Run2Data17_JPsi.root")

#MyTreeMC.Add("/eos/user/m/martis/data/ZmumuNtuples/Zmumu_MCAll_n308.root")
#MyTreeMC.Add("/eos/user/m/martis/data/ZmumuNtuples/Run2Data17_Zmumu_ftag.root")
#MyTreeMC.Add("/eos/user/m/martis/data/ZmumuNtuples/Run2Data18_Zmumu_ftag.root")

#MyTreeMC=ROOT.TChain("Refit1_SiAndTRT")

#MyTreeMC.Add("/eos/user/n/navarrjo/public/qt/weakmodes/athena/InDetAlignExample/InDetPerformanceMonitoring/WeakModes/v3_samples/merged/zmumu_nominal_josep_v3.root")
#
#MyTreeMC.Add("/eos/user/n/navarrjo/public/qt/weakmodes/athena/InDetAlignExample/InDetPerformanceMonitoring/WeakModes/v3_samples/merged/zmumu_11_n400_josep_v3.root")


    MyTreeMC=ROOT.TChain("DecaysToMuons")
    MyTreeMC.Add("RefitTrack_FullTrack_003.root")




# In[3]:


    MyTreeMC.Show(0)


# In[ ]:

   
    
    #massBinning = [600,6,18]
    massBinning = [648,3600,14400]
    pt1_min  = 30
    pt2_min  = 10

    outfile  =  ROOT.TFile("JPsiMass.root","Recreate")
    h_Mass2EtaPhi = ROOT.TH3D("Mass2EtaPhi","#eta#phi m_{#mu#mu}^{2} [GeV^{2}]",50,-2.50,2.50,50,-ROOT.TMath.Pi(),ROOT.TMath.Pi(),massBinning[0],massBinning[1],massBinning[2])
    h_Mass2EtaEta = ROOT.TH3D("Mass2EtaEta","#eta_{1}#eta_{2} m_{#mu#mu}^{2} [GeV^{2}]",50,-2.50,2.50,50,-2.50,2.50,massBinning[0],massBinning[1],massBinning[2])
    
    
    h_Az   = ROOT.TH1D("Az"," Az",1200,-5.1,5.1)
    h_Ar   = ROOT.TH1D("Ar"," Ar",1200,-5.1,5.1)
    h_SAz   = ROOT.TH1D("SAz"," Az",1200,-5.1,5.1)
    h_SAr   = ROOT.TH1D("SAr"," Ar",1200,-5.1,5.1)
    h_DAz   = ROOT.TH1D("DAz"," Az",1200,-5.1,5.1)
    h_DAr   = ROOT.TH1D("DAr"," Ar",1200,-5.1,5.1)
    h_SA   = ROOT.TH1D("SA"," Ar",1200,-5.1,5.1)
    
    h_Az_T   = ROOT.TH1D("Az_T"," Az",120,-5.1,5.1)
    h_Ar_T   = ROOT.TH1D("Ar_T"," Ar",120,-5.1,5.1)
    h_dAz_T   = ROOT.TH1D("dAz_T"," dAz",1000,-0.1,0.1)
    h_dAr_T   = ROOT.TH1D("dAr_T"," dAr",1000,-0.1,0.1)
    
    h_ArAr   = ROOT.TH2D("ArAr"," Ar Mr",120,-5.1,5.1,120,-5.1,5.1)
    h_AzAz   = ROOT.TH2D("AzAz"," Az Az",120,-5.1,5.1,120,-5.1,5.1)
    h_ArAz   = ROOT.TH2D("ArAz"," Ar Az",120,-5.1,5.1,120,-5.1,5.1)
    h_ArzArz   = ROOT.TH2D("ArzArz"," Arz Arz",120,0.4,0.6,120,0.4,0.6)    
    h_AzT   = ROOT.TH1D("AzT"," Az",120,-5.1,5.1)
    h_ArT   = ROOT.TH1D("ArT"," Ar",120,-5.1,5.1)
   

    h_Mass2T = ROOT.TH1D("Mass2T","m_{#mu#mu}^{2} [GeV^{2}]",massBinning[0],massBinning[1],massBinning[2])
    h_Mass2nT = ROOT.TH1D("Mass2nT","m_{#mu#mu}^{2} [GeV^{2}]",massBinning[0],massBinning[1],massBinning[2])
    
    i=0
    maxEvents=2e5
    

    axisS = VariableAxis( 0.491, 0.501, 1)
    axisDz = VariableAxis( -8.1, 8.1, 25)
    axisDDz = VariableAxis( -1, 1, 200)
    
    
    #axisA  = VariableAxis( -5, 5, aNbins)
    #axisMA = VariableAxis( 0.95125, 1.00125, mNbins)
    #axisDA = VariableAxis( 0., 1.0, mNbins)
    
    
    etaAxis = VariableAxis( -2.5, 2.5, 5)
    phiAxis = VariableAxis( -ROOT.TMath.Pi(), ROOT.TMath.Pi(), 5)
    
    
    histBinsX = axisDz.nBins * axisS.nBins
    histBinsY = axisDDz.nBins * axisS.nBins
    
    h_Mass2AA = ROOT.TH3D("Mass2AA","#eta_{1}#eta_{2} m_{#mu#mu}^{2} [GeV^{2}]",histBinsX,0,histBinsX,histBinsY,0,histBinsY,massBinning[0],massBinning[1],massBinning[2])
    h_Mass2AAUp1 = ROOT.TH3D("Mass2AAUp1","#eta_{1}#eta_{2} m_{#mu#mu}^{2} [GeV^{2}]",histBinsX,0,histBinsX,histBinsY,0,histBinsY,massBinning[0],massBinning[1],massBinning[2])
    h_Mass2AADn1 = ROOT.TH3D("Mass2AADn1","#eta_{1}#eta_{2} m_{#mu#mu}^{2} [GeV^{2}]",histBinsX,0,histBinsX,histBinsY,0,histBinsY,massBinning[0],massBinning[1],massBinning[2])
    h_Mass2AAUp2 = ROOT.TH3D("Mass2AAUp2","#eta_{1}#eta_{2} m_{#mu#mu}^{2} [GeV^{2}]",histBinsX,0,histBinsX,histBinsY,0,histBinsY,massBinning[0],massBinning[1],massBinning[2])
    h_Mass2AADn2 = ROOT.TH3D("Mass2AADn2","#eta_{1}#eta_{2} m_{#mu#mu}^{2} [GeV^{2}]",histBinsX,0,histBinsX,histBinsY,0,histBinsY,massBinning[0],massBinning[1],massBinning[2])
    
    #Histograms for testing
    h_Mass2AA_P = ROOT.TH3D("Mass2AA_P","#eta_{1}#eta_{2} m_{#mu#mu}^{2} [GeV^{2}]",histBinsX,0,histBinsX,histBinsY,0,histBinsY,massBinning[0],massBinning[1],massBinning[2])
    h_Mass2AA_Pt = ROOT.TH3D("Mass2AA_Pt","#eta_{1}#eta_{2} m_{#mu#mu}^{2} [GeV^{2}]",histBinsX,0,histBinsX,histBinsY,0,histBinsY,massBinning[0],massBinning[1],massBinning[2])
    h_Mass2AA_PtP = ROOT.TH3D("Mass2AA_PtP","#eta_{1}#eta_{2} m_{#mu#mu}^{2} [GeV^{2}]",histBinsX,0,histBinsX,histBinsY,0,histBinsY,massBinning[0],massBinning[1],massBinning[2])
    
    # Profiles that hold the constants
    h_Ar1 = ROOT.TProfile2D("Ar1","#eta_{1}#eta_{2} Ar1",histBinsX,0,histBinsX,histBinsY,0,histBinsY,-20,20)
    h_Ar2 = ROOT.TProfile2D("Ar2","#eta_{1}#eta_{2} Ar2",histBinsX,0,histBinsX,histBinsY,0,histBinsY,-20,20)
    h_Az1 = ROOT.TProfile2D("Az1","#eta_{1}#eta_{2} Az1",histBinsX,0,histBinsX,histBinsY,0,histBinsY,-20,20)
    h_Az2 = ROOT.TProfile2D("Az2","#eta_{1}#eta_{2} Az2",histBinsX,0,histBinsX,histBinsY,0,histBinsY,-20,20) 
    h_As1 = ROOT.TProfile2D("As1","#eta_{1}#eta_{2} As1",histBinsX,0,histBinsX,histBinsY,0,histBinsY,-20,20)
    h_As2 = ROOT.TProfile2D("As2","#eta_{1}#eta_{2} As2",histBinsX,0,histBinsX,histBinsY,0,histBinsY,-20,20)
    h_Adz2 = ROOT.TProfile2D("Adz2","#eta_{1}#eta_{2} Adz2",histBinsX,0,histBinsX,histBinsY,0,histBinsY,-20,20)
    
    h_Ar1_T = ROOT.TProfile2D("Ar1_T","#eta_{1}#eta_{2} Ar1",histBinsX,0,histBinsX,histBinsY,0,histBinsY,-20,20)
    h_Ar2_T = ROOT.TProfile2D("Ar2_T","#eta_{1}#eta_{2} Ar2",histBinsX,0,histBinsX,histBinsY,0,histBinsY,-20,20)
    h_Az1_T = ROOT.TProfile2D("Az1_T","#eta_{1}#eta_{2} Az1",histBinsX,0,histBinsX,histBinsY,0,histBinsY,-20,20)
    h_Az2_T = ROOT.TProfile2D("Az2_T","#eta_{1}#eta_{2} Az2",histBinsX,0,histBinsX,histBinsY,0,histBinsY,-20,20) 
    
   
    histHolderMass = MultiDimensionObject()
    histHolderAs1  = MultiDimensionObject()
    histHolderAs2  = MultiDimensionObject()
    histHolderAz1 = MultiDimensionObject()
    histHolderAz2 = MultiDimensionObject()
    histHolderAdz2 = MultiDimensionObject()

    for eta1 in range( 1, etaAxis.nBins+1 ) :
      for eta2 in range( eta1, etaAxis.nBins+1 ) :
        for phi1 in range(1, phiAxis.nBins+1 ) :
          for phi2 in range( phi1, phiAxis.nBins+1 ) :
            nameSub = str( eta1 ) + "_" + str( eta2 ) + "_" +  str( phi1 ) + "_" + str( phi2 )  
            histHolderMass.addValue( (eta1,eta2,phi1,phi2), ROOT.TH3D("Mass2AA_"+nameSub,"A_{1}A_{2} m_{#mu#mu}^{2} [GeV^{2}]",histBinsX,0,histBinsX,histBinsY,0,histBinsY,massBinning[0],massBinning[1],massBinning[2]) )
            histHolderAs1.addValue( (eta1,eta2,phi1,phi2), ROOT.TProfile2D("As1_"+nameSub,"#eta_{1}#eta_{2} As1",histBinsX,0,histBinsX,histBinsY,0,histBinsY,-20,20) )
            histHolderAs2.addValue( (eta1,eta2,phi1,phi2), ROOT.TProfile2D("As2_"+nameSub,"#eta_{1}#eta_{2} As2",histBinsX,0,histBinsX,histBinsY,0,histBinsY,-20,20) )
            histHolderAz1.addValue( (eta1,eta2,phi1,phi2), ROOT.TProfile2D("Az1_"+nameSub,"#eta_{1}#eta_{2} Az1",histBinsX,0,histBinsX,histBinsY,0,histBinsY,-20,20) )
            histHolderAz2.addValue( (eta1,eta2,phi1,phi2), ROOT.TProfile2D("Az2_"+nameSub,"#eta_{1}#eta_{2} Az2",histBinsX,0,histBinsX,histBinsY,0,histBinsY,-20,20) )
            histHolderAdz2.addValue( (eta1,eta2,phi1,phi2), ROOT.TProfile2D("Adz2_"+nameSub,"#eta_{1}#eta_{2} Adz2",histBinsX,0,histBinsX,histBinsY,0,histBinsY,-20,20) )
            #print( eta1, eta2, phi1, phi2 )

    
    l1 = TLorentzVector()
    l2 = TLorentzVector()
    
    
    biasPT = ROOT.TVector3(1.005,1.005,1)
    biasP = ROOT.TVector3(1.01,1.01,1.01)
    
    
    
    for event in MyTreeMC:
        if(i%50000==0):
            print( "Events processed : ", i )
        if(i>maxEvents):
            break
        i+=1


        l1.SetXYZM(event.Negative_T_Px*1e-3,event.Negative_T_Py*1e-3,event.Negative_T_Pz*1e-3,105.6583*1e-3)
        l2.SetXYZM(event.Positive_T_Px*1e-3,event.Positive_T_Py*1e-3,event.Positive_T_Pz*1e-3,105.6583*1e-3)
        Ar1_T,Az1_T,Mrz1_T = calculateAS(l1,l2)
        Ar2_T,Az2_T,Mrz2_T = calculateAS(l2,l1)
        isTruth = True
        if( event.Negative_T_Px == -999 or event.Positive_T_Px==-999):
          isTruth = False
          
 
        l1.SetXYZM(event.Negative_Px*1e-3,event.Negative_Py*1e-3,event.Negative_Pz*1e-3,105.6583*1e-3)
        l2.SetXYZM(event.Positive_Px*1e-3,event.Positive_Py*1e-3,event.Positive_Pz*1e-3,105.6583*1e-3)
        if l1.Pt() < pt1_min and l2.Pt() < pt1_min:
          continue
        if l1.Pt() < pt2_min or l2.Pt() < pt2_min:
          continue


        Ar1,Az1,Mrz1 = calculateAS(l1,l2)
        Ar2,Az2,Mrz2 = calculateAS(l2,l1)
        
        #print( Mrz1*0.5-Ar1-Az1  )
        h_SAr.Fill( Ar1 + Ar2 )
        h_SAz.Fill( Az1 + Az2 ) 
        h_DAr.Fill( Ar1 - Ar2 )
        h_DAz.Fill( Az1 - Az2 ) 
        h_SA.Fill( Az1 + Az2 +  Ar1 + Ar2 ) 
        
        h_Ar.Fill( Ar1 )
        h_Az.Fill( Az1 ) 
        h_Ar.Fill( Ar2 )
        h_Az.Fill( Az2 )
        if isTruth:        
          h_dAr_T.Fill( Ar1 - Ar1_T )
          h_dAz_T.Fill( Az1 - Az1_T ) 
          h_dAr_T.Fill( Ar2 - Ar2_T )
          h_dAz_T.Fill( Az2 - Az2_T ) 
 
          h_Ar_T.Fill( Ar1_T )
          h_Az_T.Fill( Az1_T ) 
          h_Ar_T.Fill( Ar2_T )
          h_Az_T.Fill( Az2_T )
      
        h_ArzArz.Fill( Ar1 + Az1, Ar2 + Az2)
        h_ArAr.Fill( Ar1, Ar1+Ar2 )
        h_AzAz.Fill( Az1, Az2)
        h_ArAz.Fill( Ar1, Az1)
        h_ArAz.Fill( Ar2, Az2)
    
        #binX = calcBin( axisA, Ar1 )*axisMA.nBins + calcBin( axisMA, Mrz1 )
        #binY = calcBin( axisDA, Ar1+Ar2 )*axisMA.nBins  + calcBin( axisMA, Mrz2 ) 
        binX = calcBin( axisS, Ar1+Az1 )*axisDz.nBins + calcBin( axisDz, Az1-Ar1 )
        binY = calcBin( axisS, Ar2+Az2 )*axisDDz.nBins + calcBin( axisDDz, Az2-Ar2+(Az1-Ar1) )
        #print( Ar1,Az1, binX, binXtoArAz(axisA,axisMA,binX) )
        
        Jpsi = l1 + l2
        JpsiM2 = Jpsi.M2()
        h_Mass2EtaPhi.Fill(l1.Eta(),l1.Phi(),JpsiM2)
        if l1.Eta() > l2.Eta():
            h_Mass2EtaEta.Fill( l2.Eta(), l1.Eta(), JpsiM2 )
        else :
            h_Mass2EtaEta.Fill( l1.Eta(), l2.Eta(), JpsiM2 )
    
        if isTruth:
          h_Mass2T.Fill( JpsiM2 )
        else:
          h_Mass2nT.Fill( JpsiM2 )



        h_Mass2AA.Fill( binX, binY, JpsiM2 )
        h_Mass2AAUp1.Fill( binX, binY, JpsiM2*1.004 )
        h_Mass2AAUp2.Fill( binX, binY, JpsiM2*1.040 )
        h_Mass2AADn1.Fill( binX, binY, JpsiM2*0.996 )
        h_Mass2AADn2.Fill( binX, binY, JpsiM2*0.960 )
        h_Ar1.Fill( binX,binY, Ar1 )
        h_Ar2.Fill( binX,binY, Ar2 )
        h_Az1.Fill( binX,binY, Az1 )
        h_Az2.Fill( binX,binY, Az2 )
        h_As1.Fill( binX,binY, Ar1+Az1 )
        h_As2.Fill( binX,binY, Ar2+Az2 )
        h_Adz2.Fill( binX,binY, Az2+Az1 )
        h_Ar1_T.Fill( binX,binY, Ar1_T )
        h_Ar2_T.Fill( binX,binY, Ar2_T )
        h_Az1_T.Fill( binX,binY, Az1_T )
        h_Az2_T.Fill( binX,binY, Az2_T )
    
        biasMomentum( l1,  biasP) 
        biasMomentum( l2,  biasP) 
        h_Mass2AA_P.Fill( binX, binY, (l1 + l2).M2() )
        
        l1.SetXYZM(event.Negative_Px*1e-3,event.Negative_Py*1e-3,event.Negative_Pz*1e-3,105.6583*1e-3)
        l2.SetXYZM(event.Positive_Px*1e-3,event.Positive_Py*1e-3,event.Positive_Pz*1e-3,105.6583*1e-3)
     
        biasMomentum( l1,  biasPT) 
        biasMomentum( l2,  biasPT) 
        h_Mass2AA_Pt.Fill( binX, binY, (l1 + l2).M2() )
    
     
        biasMomentum( l1,  biasP) 
        biasMomentum( l2,  biasP) 
        h_Mass2AA_PtP.Fill( binX, binY, (l1 + l2).M2() )
      

        binEta1 = etaAxis.getBin(l1.Eta())  
        binEta2 = etaAxis.getBin(l2.Eta())
        binPhi1 = phiAxis.getBin(l1.Phi())
        binPhi2 = phiAxis.getBin(l2.Phi())
        if( binEta1 < 1 or binEta2 < 1 or binEta1 > etaAxis.nBins or binEta2 > etaAxis.nBins  ):
          continue
        if( binPhi1 < 1 or binPhi2 < 1 or binPhi1 > phiAxis.nBins or binPhi2 > phiAxis.nBins  ):
          continue
        if( binEta1 < binEta2 ):
          if( binPhi1 < binPhi2 ):
            abin = (binEta1, binEta2, binPhi1, binPhi2)
          else:
            abin = (binEta1, binEta2, binPhi2, binPhi1)
        else:
          if( binPhi1 < binPhi2 ):
            abin = (binEta2, binEta1, binPhi1, binPhi2)
          else:
            abin = (binEta2, binEta1, binPhi2, binPhi1)
        #print( abin )  
        histHolderMass.readValue( abin ).Fill( binX, binY, JpsiM2 )
        histHolderAs1.readValue(  abin ).Fill( binX, binY, Ar1+Az1 )
        histHolderAs2.readValue(  abin ).Fill( binX, binY, Ar2+Az2 )
        histHolderAz1.readValue(  abin ).Fill( binX, binY, Az1 )
        histHolderAz2.readValue(  abin ).Fill( binX, binY, Az2 )
        histHolderAdz2.readValue( abin ).Fill( binX, binY, Az2 + Az1 )
    
    outfile.Write()


# In[ ]:

if __name__ == "__main__":
    main()
