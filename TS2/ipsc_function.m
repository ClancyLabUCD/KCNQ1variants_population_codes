function [dY, dati] = ipsc_function(time, Y, model_parameter_inputs)

%global time_AP V_AP run_time

%% State variable definitions:

% 1: Vm (millivolt)
% Ionic Flux: -------------------------------------------------------------
% 2: Ca_SR (millimolar)
% 3: Cai (millimolar)
% Current Gating (dimensionless):------------------------------------------
% 4: g     (inactivation in i_rel)
% 5: d     (activation in i_CaL)
% 6: f1    (inactivation in i_CaL)
% 7: f2    (not used - in i_CaL)
% 8: fCa   (calcium-dependent inactivation in i_CaL)
% 9: Xr1   (activation in i_Kr)
% 10: Xr2  (inactivation in i_Kr
% 11: Xs   (activation in i_Ks)
% 12: h    (inactivation in i_Na)
% 13: j    (slow inactivation in i_Na)
% 14: m    (activation in i_Na)
% 15: Xf   (inactivation in i_f)
% 16: s    (inactivation in i_to)
% 17: r    (activation in i_to)
% 18: dCaT (activation in i_CaT)
% 19: fCaT (inactivation in i_CaT)
% Ionic Flux Cont.:--------------------------------------------------------
% 20: Nai (millimolar)
% 21: Ki (millimolar)
% 22: Ca_ligand (millimolar)

% -------------------------------------------------------------------------------
%% Parameter inputs:

%Current parameter values:
sensitivity_conductance=model_parameter_inputs(1:16);
x_K1=model_parameter_inputs(17:22);
x_KR=model_parameter_inputs(23:33);
x_IKS= model_parameter_inputs(34:39);
xTO= model_parameter_inputs(40:50);
x_cal=model_parameter_inputs(51:61);
x_cat=model_parameter_inputs(62:72);
x_NA=model_parameter_inputs(73:86);
x_F=model_parameter_inputs(87:92);

% -------------------------------------------------------------------------------
%%  Constants for flag protocols:
cyclelength=1000; %1000ms = 1hz beating
R_clamp = 0.02; %clamps only


% for stim:
stim_flag=0;
i_stim_Amplitude = -5;   % pA/pF (in stim_mode)
i_stim_End = 10e10;   % milisecond (in stim_mode)
i_stim_PulseDuration = 3;   % milisecond (in stim_mode)
i_stim_Start = 0;   % milisecond (in stim_mode)


%--------------------------------------------------------------------------
%% Conductance Factors:
Cm=60; %pF
scalegpk=0;
scale_Ical_Fca_Cadep=1.2;
scale_nonCa_CaL=1;
scaleKup=0.702;

x_scale_conductance=[ 1      %1: g_k1
    1      %2: g_kr
    1      %3: g_Ks
    1      %4: g_to
    1      %5: g_CaL
    1      %6: g_CaT
    1      %7: g_Na
    1      %8: g_f
    1.1    %9: g_naca
    0.26   %10: i_up
    1      %11: i_rel
    0.02   %12: i_leak
    1.818  %13: i_nak
    1.5    %14: i_bna
    0.62   %15: i_bca
    10.5];  %16: i_pca


size_input=size(sensitivity_conductance);
if  size_input(2)==1
    x_scale_conductance=x_scale_conductance.*sensitivity_conductance;
elseif size_input(1)==1
    x_scale_conductance=x_scale_conductance.*sensitivity_conductance';
end


dY=zeros(size(Y));


% -------------------------------------------------------------------------------
%% Cell geometry
%V_total data from Hwang et al., V_c and V_SR  proportionally scaled from Ten Tusscher 2004 values

V_tot=3960; %um^3 from hwang et al.
Vc_tenT=16404; VSR_tenT=1094; V_tot_tenT=Vc_tenT+VSR_tenT;
Vc=V_tot*(Vc_tenT/V_tot_tenT); %=3712.4 um^3 (93.7% total volume)
V_SR=V_tot*(VSR_tenT/V_tot_tenT);%=247.6 um^3 (6.3% total volume)

% -------------------------------------------------------------------------------
%% Constants
T = 310.0;   % kelvin (in model_parameters)'

R = 8.314472;   % joule_per_mole_kelvin (in model_parameters)
F = 96.4853415;   % coulomb_per_mmole (in model_parameters)

Ko = 5.4;   % millimolar (in model_parameters)
Cao = 1.8;   % millimolar (in model_parameters
Nao = 140.0;   % millimolar (in model_parameters)

% -------------------------------------------------------------------------------
%% Reversal Potentials:

E_Ca = 0.5*R*T/F*log(Cao/Y(3)); % millivolt
E_Na = R*T/F*log(Nao/Y(20));  % millivolt
E_K = R*T/F*log(Ko/Y(21));  % millivolt

% -------------------------------------------------------------------------------
%% Inward Rectifier K+ current (Ik1):

%define parameters from x_K1
xK11=x_K1(2); xK12=x_K1(3); xK13=x_K1(4);
xK14=x_K1(5); xK15=x_K1(6);

alpha_xK1=xK11.*exp((Y(1)+xK13)./xK12);
beta_xK1=exp((Y(1)+xK15)./xK14);
XK1_inf=alpha_xK1./(alpha_xK1+beta_xK1);

%Current:
g_K1=x_K1(1);
i_K1 = g_K1*XK1_inf*(Y(1)-E_K)*sqrt(Ko/5.4);

%-------------------------------------------------------------------------------
%% Rapid Delayed Rectifier Current (Ikr):

%define parameters from x_KR
Xr1_1=x_KR(2); Xr1_2=x_KR(3); Xr1_5=x_KR(4); Xr1_6=x_KR(5);
Xr2_1=x_KR(6); Xr2_2=x_KR(7); Xr2_5=x_KR(8); Xr2_6=x_KR(9);

%parameter-dependent values:
Xr1_3=Xr1_5*Xr1_1; Xr2_3=Xr2_5*Xr2_1;
Xr1_4=1/((1/Xr1_2)+(1/Xr1_6)); Xr2_4=1/((1/Xr2_2)+(1/Xr2_6));

% 9: Xr1 (dimensionless) (activation in i_Kr_Xr1)
alpha_Xr1=Xr1_1.*exp((Y(1))./Xr1_2);
beta_Xr1=Xr1_3.*exp((Y(1))./Xr1_4);
Xr1_inf=alpha_Xr1./(alpha_Xr1+ beta_Xr1);
tau_Xr1=((1./(alpha_Xr1+ beta_Xr1))+x_KR(10));

dY(9) = (Xr1_inf-Y(9))./tau_Xr1;

% 10: Xr2 (dimensionless) (inactivation in i_Kr_Xr2)
alpha_Xr2=Xr2_1.*exp((Y(1))./Xr2_2);
beta_Xr2=Xr2_3.*exp((Y(1))./Xr2_4);
Xr2_inf=alpha_Xr2./(alpha_Xr2+beta_Xr2);
tau_Xr2=((1./(alpha_Xr2+beta_Xr2))+x_KR(11));
dY(10) = (Xr2_inf-Y(10))./tau_Xr2;




%Current:
g_Kr = x_KR(1); % nS_per_pF (in i_Kr)
i_Kr = g_Kr*(Y(1)-E_K)*Y(9)*Y(10)*sqrt(Ko/5.4);

%----------------------------------------------------------------------------
%% Slow delayed rectifier current (IKs):

%define parameters from x_IKS:
ks1=x_IKS(2); ks2=x_IKS(3); ks5=x_IKS(4); ks6=x_IKS(5);
tauks_const=x_IKS(6);

%parameter-dependent values:
ks3=ks5*ks1; ks4=1/((1/ks2)+(1/ks6));

% 11: Xs (dimensionless) (activation in i_Ks)
alpha_Xs=ks1.*exp((Y(1))./ks2);
beta_Xs=ks3.*exp((Y(1))./ks4);
Xs_inf=alpha_Xs./(alpha_Xs+beta_Xs);
tau_Xs=(1./(alpha_Xs+beta_Xs))+ tauks_const;
dY(11) = (Xs_inf-Y(11))./tau_Xs;

E_Ks=E_K; %original baseline
%PKNa=0.01833;
%E_Ks=(R*T/F)*log((Ko+PKNa*Nao)/(Y(21)+PKNa*Y(20))); %new IKs

%Current:
g_Ks = x_IKS(1);  % nS_per_pF (in i_Ks)
i_Ks = g_Ks*(Y(1)-E_Ks)*(Y(11).^2);

%-------------------------------------------------------------------------------
%% Transient outward current (Ito):

%define parameters from xTO
r1=xTO(2); r2=xTO(3); r5=xTO(4); r6=xTO(5);
s1=xTO(6); s2=xTO(7); s5=xTO(8); s6=xTO(9);
tau_r_const=xTO(10);
tau_s_const=xTO(11);

%parameter-dependent values:
r3=r5*r1; r4=1/((1/r2)+(1/r6));
s3=s5*s1; s4=1/((1/s2)+(1/s6));

% 16: s (dimensionless) (inactivation in i_to)
% q gate in Paci
alpha_s=s1.*exp((Y(1))./s2);
beta_s=s3.*exp((Y(1))./s4);
s_inf=alpha_s./(alpha_s+beta_s);
tau_s=((1./(alpha_s+beta_s))+tau_s_const);
dY(16) = (s_inf-Y(16))./tau_s;


% 17: r (dimensionless) (activation in i_to)
alpha_r=r1.*exp((Y(1))./r2);
beta_r=r3.*exp((Y(1))./r4);
r_inf=alpha_r./(alpha_r+ beta_r);
tau_r=(1./(alpha_r+ beta_r))+tau_r_const;
dY(17) = (r_inf-Y(17))./tau_r;

%Current:
g_to = xTO(1);% nS_per_pF (in i_to)
i_to = g_to*(Y(1)-E_K)*Y(16)*Y(17);

%-------------------------------------------------------------------------------
%% L-type Ca2+ current (ICaL):

%define parameters from x_cal
d1=x_cal(2);d2=x_cal(3);d5=x_cal(4);d6=x_cal(5);
f1=x_cal(6);f2=x_cal(7);f5=x_cal(8); f6=x_cal(9);
taud_const=x_cal(10);
tauf_const=x_cal(11);

%parameter-dependent values:
d3=d5*d1;d4=1/((1/d2)+(1/d6));
f3=f5*f1;f4=1/((1/f2)+(1/f6));

% 5: d (dimensionless) (activation in i_CaL)
alpha_d=d1.*exp(((Y(1)))./d2);
beta_d=d3.*exp(((Y(1)))./d4);
d_inf=alpha_d./(alpha_d+ beta_d);
tau_d=((1./(alpha_d+ beta_d))+taud_const);
dY(5)= (d_inf-Y(5))/tau_d;

% 6: f1 (dimensionless) (inactivation  i_CaL)
alpha_f=f1.*exp(((Y(1)))./f2);
beta_f=f3.*exp(((Y(1)))./f4);
f_inf=alpha_f./(alpha_f+beta_f);
tau_f=((1./(alpha_f+beta_f)) + tauf_const);
dY(6)= (f_inf-Y(6))/tau_f;

% 7: EMPTY (not used - in i_CaL)
dY(7)=0;

% 8: fCa (dimensionless) (calcium-dependent inactivation in i_CaL)
% from Ten tusscher 2004
alpha_fCa = 1.0/(1.0+((scale_Ical_Fca_Cadep.*Y(3))/.000325)^8.0);
beta_fCa = 0.1/(1.0+exp((scale_Ical_Fca_Cadep.*Y(3)-.0005)/0.0001));
gamma_fCa = .2/(1.0+exp((scale_Ical_Fca_Cadep.*Y(3)-0.00075)/0.0008));

fCa_inf =  ((alpha_fCa+beta_fCa+gamma_fCa+.23)/(1.46));
tau_fCa=2; %ms
if fCa_inf>Y(8) && Y(1)>-60
    k_fca=0;
else
    k_fca=1;
end
dY(8) = k_fca.*(fCa_inf-Y(8))/tau_fCa;


%Current:
p_CaL =  x_cal(1); % nS_per_pF (in i_CaL)
p_CaL_shannonCa=5.4e-4;
p_CaL_shannonNa=1.5e-8;
p_CaL_shannonK=2.7e-7;
p_CaL_shannonTot=p_CaL_shannonCa + p_CaL_shannonNa + p_CaL_shannonK;
p_CaL_shannonCap=p_CaL_shannonCa/p_CaL_shannonTot;
p_CaL_shannonNap=scale_nonCa_CaL*p_CaL_shannonNa/p_CaL_shannonTot;
p_CaL_shannonKp=scale_nonCa_CaL*p_CaL_shannonK/p_CaL_shannonTot;


p_CaL_Ca=p_CaL_shannonCap*p_CaL;
p_CaL_Na=p_CaL_shannonNap*p_CaL;
p_CaL_K=p_CaL_shannonKp*p_CaL;

ibarca= p_CaL_Ca*4.0*Y(1)*F^2.0/(R*T)* (.341*Y(3)*exp(2.0*Y(1)*F/(R*T))-0.341*Cao)/(exp(2.0*Y(1)*F/(R*T))-1.0);
i_CaL_Ca =  ibarca *Y(5)*Y(6)*Y(8);

ibarna= p_CaL_Na*Y(1)*F^2.0/(R*T)* (.75*Y(20)*exp(Y(1)*F/(R*T))-0.75*Nao)/(exp(Y(1)*F/(R*T))-1.0);
i_CaL_Na=  ibarna *Y(5)*Y(6)*Y(8);

ibark= p_CaL_K*Y(1)*F^2.0/(R*T)* (.75*Y(21)*exp(Y(1)*F/(R*T))-0.75*Ko)/(exp(Y(1)*F/(R*T))-1.0);
i_CaL_K = ibark *Y(5)*Y(6)*Y(8);


i_CaL=i_CaL_Ca+i_CaL_Na+i_CaL_K;
%-------------------------------------------------------------------------------
%% T-type Calcium Current (ICaT):

%SAN T-TYPE CA2+ model (Demir et al., Maltsev-Lakatta ), G_CaT determined by fit to Kurokawa IV:
dcat_inf= 1./(1+exp(-((Y(1)) +26.3)./6));
tau_dcat=1./(1.068*exp(((Y(1))+26.3)./30)+ 1.068*exp(-((Y(1))+26.3)./30));
dY(18) = (dcat_inf-Y(18))/tau_dcat;

fcat_inf= 1./(1+exp(((Y(1)) +61.7)./5.6));
tau_fcat=1./(.0153*exp(-((Y(1))+61.7)./83.3)+ 0.015*exp(((Y(1))+61.7)./15.38));
dY(19) = (fcat_inf-Y(19))/tau_fcat;

g_CaT=x_cat(1); % nS_per_pF (in i_CaT)
i_CaT= g_CaT*(Y(1)-E_Ca)*Y(18)*Y(19);

% -------------------------------------------------------------------------------
%% Sodium Current (INa):

%define parameters from x_Na
m1=x_NA(2); m2=x_NA(3); m5=x_NA(4); m6= x_NA(5);
h1=x_NA(6); h2=x_NA(7); h5=x_NA(8); h6=x_NA(9);
j1=x_NA(10); j2=x_NA(11);
tau_m_const=x_NA(12);
tau_h_const=x_NA(13);
tau_j_const=x_NA(14);

%parameter-dependent values:
m3=m5*m1; m4=1/((1/m2)+(1/m6));
h3=h5*h1; h4=1/((1/h2)+(1/h6));
j5=h5; j6=h6;
j3=j5*j1; j4=1/((1/j2)+(1/j6));

% % 12: h (dimensionless) (inactivation in i_Na)
alpha_h=h1.*exp((Y(1))./h2);
beta_h=h3.*exp((Y(1))./h4);
h_inf=(alpha_h./(alpha_h+beta_h));
tau_h=((1./(alpha_h+beta_h))+tau_h_const);
dY(12)=(h_inf-Y(12))./tau_h;

% % 13: j (dimensionless) (slow inactivation in i_Na)
alpha_j=j1.*exp((Y(1))./j2);
beta_j=j3.*exp((Y(1))./j4);
j_inf=(alpha_j./(alpha_j+beta_j));
tau_j=((1./(alpha_j+beta_j))+tau_j_const);
dY(13)=(j_inf-Y(13))./tau_j;

% % 14: m (dimensionless) (activation in i_Na)
alpha_m=m1.*exp((Y(1))./m2);
beta_m=m3.*exp((Y(1))./m4);
m_inf=alpha_m./(alpha_m+beta_m);
tau_m=((1./(alpha_m+beta_m))+tau_m_const);
dY(14) =(m_inf-Y(14))./tau_m;

%Current:
g_Na=x_NA(1); % nS_per_pF (in i_Na)
i_Na = g_Na*Y(14)^3.0*Y(12)*(Y(1)-E_Na)*Y(13);


%-------------------------------------------------------------------------------%-------------------------------------------------------------------------------
%% Funny/HCN current (If):

%define parameters from x_F
xF1=x_F(2); xF2=x_F(3); xF5=x_F(4); xF6=x_F(5);
xF_const=x_F(6);

%parameter-dependent values:
xF3=xF5*xF1;  xF4=1/((1/xF2)+(1/xF6));

% 15: Xf (dimensionless) (inactivation in i_f)
alpha_Xf=xF1.*exp((Y(1))./xF2);
beta_Xf=xF3.*exp((Y(1))./xF4);
Xf_inf=alpha_Xf./(alpha_Xf+beta_Xf);
tau_Xf=((1./(alpha_Xf+beta_Xf))+xF_const);
dY(15) = (Xf_inf-Y(15))./tau_Xf;

%Current:
g_f=x_F(1); % nS_per_pF (in i_f)
NatoK_ratio=.491; %Verkerk et al. 2013
Na_frac=NatoK_ratio./(NatoK_ratio+1);
i_fNa=Na_frac*g_f*Y(15)*(Y(1)-E_Na);
i_fK=(1-Na_frac)*g_f*Y(15)*(Y(1)-E_K);
i_f=i_fNa+i_fK;


%-------------------------------------------------------------------------------
%% Na+/Ca2+ Exchanger current (INaCa):
% Ten Tusscher formulation

KmCa = 1.38;   % Cai half-saturation constant millimolar (in i_NaCa)
KmNai = 87.5;   % Nai half-saturation constnat millimolar (in i_NaCa)
Ksat = 0.1;   % saturation factor dimensionless (in i_NaCa)

gamma = 0.35*2;   % voltage dependence parameter dimensionless (in i_NaCa)
alpha = 2.5*1.1;    %factor to enhance outward nature of inaca dimensionless (in i_NaCa)

    kNaCa = 1000*1.1;  % maximal inaca pA_per_pF (in i_NaCa)
i_NaCa = kNaCa*((exp(gamma*Y(1)*F/(R*T))*(Y(20)^3.0)*Cao)-(exp((gamma-1.0)*Y(1)*F/(R*T))*(Nao^3.0)*Y(3)*alpha))/(((KmNai^3.0)+(Nao^3.0))*(KmCa+Cao)*(1.0+Ksat*exp((gamma-1.0)*Y(1)*F/(R*T))));


%-------------------------------------------------------------------------------
%% Na+/K+ pump current (INaK):
% Ten Tusscher formulation

Km_K = 1.0;   % Ko half-saturation constant millimolar (in i_NaK)
Km_Na = 40.0;   %  Nai half-saturation constant millimolar (in i_NaK)
PNaK = 1.362*1.818;   % maxiaml nak pA_per_pF (in i_NaK)
i_NaK = PNaK*((Ko*Y(20))/((Ko+Km_K)*(Y(20)+Km_Na)*(1.0+0.1245*exp(-0.1*Y(1)*F/(R*T))+0.0353*exp(-Y(1)*F/(R*T)))));

%-------------------------------------------------------------------------------
%% SR Uptake/SERCA (i_up):
% Ten Tusscher formulation

Kup = 0.00025*scaleKup;   % millimolar (in calcium_dynamics)

    VmaxUp =0.000425;   % millimolar_per_milisecond (in calcium_dynamics)


i_up =0.26*VmaxUp/(1.0+Kup^2.0/Y(3)^2.0);

%-------------------------------------------------------------------------------
%% SR Leak (i_leak):
% Ten Tusscher formulation

V_leak = 0.02*0.00008;   % per_millisecond (in calcium_dynamics)

i_leak = (Y(2)-Y(3))*V_leak;

%-------------------------------------------------------------------------------
%% SR Release/RYR (i_rel):


ks = 12.5*x_scale_conductance(11); % [1/ms]
koCa = 56320*11.43025;                  % [mM^-2 1/ms]   %default 10   modified 20
kiCa = 54*0.3425;                 % [1/mM/ms]
kom = 1.5*0.1429;                      % [1/ms]
kim = 0.001*0.5571;                    % [1/ms]
ec50SR = 0.45;
MaxSR = 15;
MinSR = 1;


kCaSR = MaxSR - (MaxSR-MinSR)/(1+(ec50SR/Y(2))^2.5);
koSRCa = koCa/kCaSR;
kiSRCa = kiCa*kCaSR;
RI = 1-Y(23)-Y(24)-Y(25);

dY(23) =(kim*RI-kiSRCa*Y(3)*Y(23))-(koSRCa*Y(3)^2*Y(23)-kom*Y(24));   % R
dY(24) =(koSRCa*Y(3)^2*Y(23)-kom*Y(24))-(kiSRCa*Y(3)*Y(24)-kim*Y(25));% O
dY(25) =(kiSRCa*Y(3)*Y(24)-kim*Y(25))-(kom*Y(25)-koSRCa*Y(3)^2*RI);   % I

i_rel= ks*Y(24)*(Y(2)-Y(3))*(V_SR/Vc);


%% Background Sodium (I_bNa):
% Ten Tusscher formulation

g_b_Na = .00029*1.5;   % nS_per_pF (in i_b_Na)
i_b_Na = g_b_Na*(Y(1)-E_Na);

%-------------------------------------------------------------------------------
%% Background Calcium (I_bCa):
% Ten Tusscher formulation

g_b_Ca = .000592*0.62;   % nS_per_pF (in i_b_Ca)
i_b_Ca = g_b_Ca*(Y(1)-E_Ca);
%-------------------------------------------------------------------------------
%% Calcium SL Pump (I_pCa):
% Ten Tusscher formulation

g_PCa = 0.025*10.5;   % pA_per_pF (in i_PCa)
KPCa = 0.0005;  % millimolar (in i_PCa)
i_PCa = g_PCa*Y(3)/(Y(3)+KPCa);


%-------------------------------------------------------------------------------
%% 22: Calcium Ligand Buffering
%Shannon-Bers Buffering for Fluo-3

dY(22)=0;

%-------------------------------------------------------------------------------
%% 3: Cai (millimolar)
%Sobie TT/ Paci form (rapid equilibrium approximation equations) -- not as formulated in ten Tusscher 2004 text

Buf_C =.06;%.15*.45;   % (.15mM in TT) millimolar (in calcium_dynamics)
Kbuf_C = .0006;%.001*1.6;   % (.001mM in TT) millimolar (in calcium_dynamics)

Cai_bufc = 1/(1.0+Buf_C*Kbuf_C/(Y(3)+Kbuf_C)^2.0);

dY(3) = (Cai_bufc)*(i_leak-i_up+i_rel-dY(22)-(i_CaL_Ca+i_CaT+i_b_Ca+i_PCa-2*i_NaCa)*Cm/(2.0*Vc*F));

%-------------------------------------------------------------------------------
%% 2: CaSR (millimolar)
%Sobie TT/ Paci form (rapid equilibrium approximation equations) -- not as formulated in ten Tusscher 2004 text

Buf_SR = 10.0*1.2;%*(Buf_C/.15);   % (10mM in TT) millimolar (in calcium_dynamics)
Kbuf_SR = 0.3;%*(Kbuf_C/.001); %(.3 mM in TT) millimolar (in calcium_dynamics)
Ca_SR_bufSR = 1/(1.0+Buf_SR*Kbuf_SR/(Y(2)+Kbuf_SR)^2.0);

dY(2) = Ca_SR_bufSR*Vc/V_SR*(i_up-(i_rel+i_leak));


%-------------------------------------------------------------------------------
%% 20: Nai (millimolar) (in sodium_dynamics)

dY(20) = -Cm*(i_Na+i_b_Na+i_fNa+3.0*i_NaK+3.0*i_NaCa +i_CaL_Na)/(F*Vc);
%dY(20)=0;



%-------------------------------------------------------------------------------
%% 21: Ki (millimolar) (in potatssium_dynamics)
%note: no i_ax, no i_pK
GpK= scalegpk*0.00146;
i_pK = GpK*(Y(1)-E_K)/(exp(-(Y(1)-25)/5.98) + 1) ;

dY(21) = -Cm*(i_K1+i_to+i_Kr+i_Ks+i_fK+i_pK -2.*i_NaK + i_CaL_K)/(F*Vc); %+ i_stim
%dY(21)=0;


%-------------------------------------------------------------------------------
%% 1: Vm (volt) (in Membrane)

% I_stim:
if (time >= i_stim_Start) && (time <= i_stim_End) && (mod(time, cyclelength)<i_stim_PulseDuration)
    i_stim = stim_flag*i_stim_Amplitude;%/Cm;
else
    i_stim = 0.0;
end

v_clamp=Y(1); %set i_voltageclamp to 0

i_voltageclamp=(v_clamp-Y(1))/R_clamp;


dY(1, 1) = -(i_K1+i_to+i_Kr+i_Ks+i_CaL+i_CaT+i_NaK+i_Na+i_NaCa+i_PCa+i_f+i_b_Na+i_b_Ca+i_stim- i_voltageclamp);

%% -------------------------------------------------------------------------------
dati = [i_K1, i_to, i_Kr, i_Ks, i_CaL, i_NaK, i_Na, i_NaCa, i_PCa, i_f, i_b_Na, i_b_Ca, i_rel, i_up, i_leak, i_stim, i_CaT, E_K, E_Na];


end
%===============================================================================
% End of file
%===============================================================================
