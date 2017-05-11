class EoS ;

class TransportCoeff
{
  double etaS_min, etaS_slope_QGP, etaS_slope_HRG, zetaS, taupi, tauPi ;
  bool isvis;
 EoS *eos ;
 public:
 TransportCoeff(double _etaS, double _etaS_slope_QGP, double _etaS_slope_HRG, double _zetaS, EoS *_eos) ;
 ~TransportCoeff() {} ;
 void getEta(double e, double T, double &_etaS, double &_zetaS) ;
 void getTau(double T, double &_taupi, double &_tauPi) ;
 inline bool isViscous() { return isvis;}
} ;
