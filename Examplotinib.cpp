$PROB
  
# Model: `Examplotinib`
  
  - Author: "Oliver"
  - Date: "02.11.2019"
  - Source: "My Mind"
  
$SET delta=0.25
  
$PARAM @annotated

SEX : 0 : 0 = male, 1 = female
WT : 75 : Body surface area (kg)
CLCR : 115 : creatinine clearance (mL/min)
    
$THETA @annotated
0.22   : typical KA (1/h)
0.097  : typical KE (1/h)
35     : typical V (L)
0.8    : typical fraction of female clearance
0.65   : typical fraction of female central volume
0.75   : exponent for wt on volume
450    : volume of peripheral compartment (L)
5     : intercompartmental clearance (L/hr)

  
$MAIN


double KA_ind = THETA1*exp(ETA_KA);
double KE_ind = THETA2*(CLCR/115)*pow(THETA4, SEX) *exp(ETA_KE);

double V_ind = THETA3*pow(WT/75,THETA6)*pow(THETA5, SEX)*exp(ETA_V);

double CL_ind = KE_ind * V_ind;  

double VP_ind = THETA7 * exp(ETA_V2);
double Q_ind = THETA8;

F_GUT = 1.0;
ALAG_GUT = 0.75;
  
  
$OMEGA @block @annotated
ETA_V  : 0.15             : ETA on V
ETA_KE : 0.15 0.36        : ETA on ke
ETA_KA : 0.00 0.00 0.75   : ETA on ka
ETA_V2 : 0.00 0.00 0.00 0.32 : ETA on VP
    
$SIGMA @annotated
PROP : 0.06   : proportional error
ADD  : 0.005   : additive error
    
    
$INIT
GUT = 0, CENT = 0, PER = 0
  
$ODE
dxdt_GUT  = - KA_ind * GUT;
dxdt_CENT = + KA_ind * GUT - KE_ind*CENT-(Q_ind/V_ind)*CENT+(Q_ind/VP_ind)*PER;
dxdt_PER = +(Q_ind/V_ind)*CENT-(Q_ind/VP_ind)*PER;


$TABLE
capture IPRED = CENT/V_ind;

double DV = IPRED * (1+PROP) + ADD;
  while (DV < 0) {
    simeps();
    DV = IPRED*(1+PROP) + ADD;
  }
  
$CAPTURE @annotated
DV    : Plasma concentration (mg/L)
KA_ind : individual ka (1/hr)
CL_ind: individual Clearance (L/hr)
V_ind : individual central volume (L)
KE_ind: individual ke (1/hr)
VP_ind : individual V_per (L)
Q_ind : individual Q (L/hr)
WT : individual WT (kg)
CLCR: individual CLCR (mL/min)
SEX : indivdual sex