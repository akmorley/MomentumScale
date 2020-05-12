from ROOT import TLorentzVector

def biasMomentum(  tlvec, scale ):
  vec3 = tlvec.Vect()
  vec3[0] *= scale(0)
  vec3[1] *= scale(1)
  vec3[2] *= scale(2)
  tlvec.SetXYZM( vec3.x(), vec3.y(), vec3.z(), tlvec.M() )


def calculateASM( pos, neg, M2 ):
  Ep = pos.Energy()
  En = neg.Energy()
  M2p = pos.M2()
  M2n = neg.M2()
  px_p = pos.Vect().x()
  px_n = neg.Vect().x()
  py_p = pos.Vect().y()
  py_n = neg.Vect().y()
  pz_p = pos.Vect().z()
  pz_n = neg.Vect().z()
  
  eRatio = En/Ep
  AT = ( eRatio * ( px_p*px_p + py_p*py_p ) - px_p*px_n - py_p*py_n ) / M2;
  AZ = ( eRatio * ( pz_p*pz_p ) - pz_p*pz_n ) / M2;
  ATPAZ = (M2-M2p*(1+2*eRatio)-M2n)/M2 
  return (AT, AZ, ATPAZ)



def calculateAM_PT( pos, neg, M2 ):
  Ep = pos.Energy()
  En = neg.Energy()
  px_p = pos.Vect().x()
  px_n = neg.Vect().x()
  py_p = pos.Vect().y()
  py_n = neg.Vect().y()

  A = ( En/Ep * ( px_p*px_p + py_p*py_p ) - px_p*px_n - py_p*py_n ) / M2;
  return A


def calculateAM_PZ(pos, neg, M2):
  Ep = pos.Energy()
  En = neg.Energy()
 
  pz_p = pos.Vect().z()
  pz_n = neg.Vect().z()
 
  A = ( En/Ep * ( pz_p*pz_p ) - pz_p*pz_n ) / M2
  return A

def calculateAS( pos, neg ):
  M2 = (pos + neg).M2()
  return calculateASM( pos, neg, M2) 

def calculateA_PT( pos, neg ):
  M2 = (pos + neg).M2()
  return calculateAM_PT( pos, neg, M2) 


def calculateA_PZ(pos, neg):
  M2 = (pos + neg).M2()
  return calculateAM_PZ( pos, neg, M2)




