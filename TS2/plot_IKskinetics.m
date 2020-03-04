clear
close all
clc

mutations_expinputs=readtable('Vanoye_TS2.xlsx');
mutation_names=table2array(mutations_expinputs(:,1));

N_mutations= length(mutation_names);
%N_mutations=5;
muts_to_analyze=1:length(mutation_names);

Vm_plot=-100:5:100;
Vm_plot_IV=-100:10:100;

load('KCNQ1_mutation_allVariables_WT');
x_ks_fit = wt_parameters(1,34:39);
act_inf_fit=inf_gate(x_ks_fit(2:5), Vm_plot);
[Iks_IV_fit] = run_IKs_IVcurve(x_ks_fit,Vm_plot_IV);

x_ks_fit_all=zeros(length(mutation_names)+1,6);
 x_ks_fit_all(1,:)=x_ks_fit;

plot_IV=1;

figure(1),set(gcf,'color','w')
set(gca,'box','off','tickdir','out','fontsize',24, 'LineWidth', 2)
hold on
plot(Vm_plot, act_inf_fit,'k','LineWidth',3);

figure(2),set(gcf,'color','w')
set(gca,'box','off','tickdir','out','fontsize',24, 'LineWidth', 2)
hold on
plot(Vm_plot_IV, Iks_IV_fit,'k','LineWidth',3);

for mt=muts_to_analyze
    name_test=mutation_names(mt);
    
    filename=char(strcat('KCNQ1_mutation_allVariables_', name_test));
    load(filename)
   
    x_ks_fit = mut_parameters(1,34:39);
    x_ks_fit_all(mt+1,:)=x_ks_fit;
    %% Calculations for plotting
    act_inf_fit=inf_gate(x_ks_fit(2:5), Vm_plot);

    %plot steady state
    figure(1);plot(Vm_plot, act_inf_fit,'LineWidth',3)
    hold on
      
    if plot_IV==1
    % Plot IV
    [Iks_IV_fit] = run_IKs_IVcurve(x_ks_fit,Vm_plot_IV);
    figure(2);plot(Vm_plot_IV, Iks_IV_fit,'LineWidth',3);
    hold on
    end
    
end

figure(1)
legend(['WT';mutation_names(muts_to_analyze)]);
legend boxoff
xlabel('Voltage (mV)','FontSize',24); ylabel('steady state activation (X_s^2','FontSize',24);


figure(2)
legend(['WT';mutation_names(muts_to_analyze)]);
legend boxoff
 xlabel('Voltage (mV)','FontSize',24); ylabel('I_{Ks} (pA/pF)','FontSize',24);

 
 %%
 
 function [ iks_IV_plot, Time, IKs ] = run_IKs_IVcurve(x,Vm_step)%, protocol)

options = odeset('MaxStep',1,'InitialStep',2e-5);

Y=[-80; 0 ];

period=[ 1000, 1990, 1000];

Vm_mat=[-80.*ones(size(Vm_step')), Vm_step', -80.*ones(size(Vm_step')), ];

xy=length(Y); %variables in y

values=zeros(sum(period)*100, length(Vm_step')*xy);
Time=zeros(sum(period)*100, length(Vm_step'));
IV_place=zeros(size(Vm_step));
for j=1:length(Vm_step')

    Time_forVmi=0;
    Y_init=Y;
    for i=1:length(period)
        Y_init(1)=(Vm_mat(j,i))';
        
        [Time_hold, values_hold] = ode15s(@IKs_ODE,[Time_forVmi(end), Time_forVmi(end)+period(i)],Y_init, options, x);%, protocol);
        Y_init=values_hold(end,:)';
        
        if i==1
            Time_forVmi=Time_hold;
            values_forVmi=values_hold;
        else
            values_forVmi=[values_forVmi; values_hold];
            Time_forVmi=[Time_forVmi; Time_hold];
        end
        if i==2
            IV_place(j)=length(Time_forVmi);
        end
    end
    
    Time(1:length(Time_forVmi),j)= Time_forVmi;
    values(1:length(Time_forVmi),((j*xy)-(xy-1)):(j*(xy)))= values_forVmi;
    
end
Time(all(values==0,2),:)=[];
values(all(values==0,2),:)=[];
Yc=values;
t=Time;


IKs=zeros(size(t));

for j=1:size(t,2)
    for i= 1:size(Yc,1)
        [~, IKs_temp] =  IKs_ODE(t(i,j), Yc(i,((j*xy)-(xy-1)):(j*(xy))),x);%, protocol);
        IKs(i,j) = IKs_temp;
    end
end

iks_IV_plot=zeros(size(Vm_step));
for i=1:length(Vm_step)
    iks_IV_plot(i)=IKs(IV_place(i),i);
end

 end


 function [dY, i_Ks, alpha_ks, beta_ks] = IKs_ODE(time, Y, x_iks, protocol)

%-------------------------------------------------------------------------------
% State variables
%-------------------------------------------------------------------------------

% 1: Vm (volt) (in Membrane)
% 2: Xs1 (dimensionless) 


%-------------------------------------------------------------------------------
%Constants
%-------------------------------------------------------------------------------

g_Ks = x_iks(1);   % nS_per_pF (in i_Kr)

R = 8.314472;   % joule_per_mole_kelvin (in model_parameters)
T = 310.0;   % kelvin (in model_parameters)'
F = 96.4853415;   % coulomb_per_mmole (in model_parameters)

ks1=x_iks(2);
ks2=x_iks(3);
ks5=x_iks(4);
ks6=x_iks(5);
ks3=ks5*ks1;
ks4=1/((1/ks2)+(1/ks6));

alpha_ks=ks1.*exp((Y(1)+20)./ks2);
beta_ks=ks3.*exp((Y(1)+20)./ks4);


Ko=4;
Ki=110;
Nao=140;
Nai=10;


    tau_ks=(1./(alpha_ks+beta_ks))+x_iks(6);
    ks_inf=alpha_ks./(alpha_ks+beta_ks);

    PKNa=0.01833;
    E_Ks=(R*T/F)*log((Ko+PKNa*Nao)/(Ki+PKNa*Nai));

dY(2, 1) = (ks_inf-Y(2))./tau_ks;

 
i_Ks = g_Ks*(Y(1)-E_Ks)*(Y(2).^2);



dY(1, 1) = 0; %voltage clamp

 end

 
 function [ x_inf ] = inf_gate(  var, Vm )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


x_inf=1./(1+(var(3).*exp((Vm)./var(4))));
x_inf=x_inf.^2;


 end


