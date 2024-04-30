$PARAM
TVCL=17.9, TVV1=22.9, TVQ=56.8, TVV2=34.3, R=0.764, TF=0.655, LBW=64.3, FM=41.5, TVBMAX=247, KD=65.3

$CMT
CENT P AUCPL AUCISF

$MAIN
double CL=TVCL*exp(ETA(1));
double V1=TVV1*pow(LBW/64.3, 1)*exp(ETA(2));
double Q=TVQ*pow(LBW/64.3, 0.75)*exp(ETA(3));
double V2=TVV2*((LBW/64.3)*R+(FM/41.5)*(1-R))*exp(ETA(4));

double BMAX=TVBMAX*exp(ETA(5));

double k10 = CL/V1;
double k12 = Q/V1;
double k21 = Q/V2;

$OMEGA 0 0 0 0 0

$SIGMA @labels PROP ADD
0 0

$ODE
dxdt_CENT= -k12*CENT -k10*CENT +k21*P;
dxdt_P= -k21*P +k12*CENT;

dxdt_AUCPL=CENT/V1;
dxdt_AUCISF=(P/V2)*TF;

$TABLE
double CU=CENT/V1;
double CT=CU+(CU*BMAX/(CU+KD));
double C2=P/V2;
double CISF=(P/V2)*TF;

double AUCPLO=AUCPL;
double AUCISFO=AUCISF;
double PI=AUCISF/AUCPL;

$CAPTURE
CU CT CISF C2 AUCPLO AUCISFO PI LBW FM 
