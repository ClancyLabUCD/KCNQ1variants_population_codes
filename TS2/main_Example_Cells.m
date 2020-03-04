clear
clc
close all

options = odeset('MaxStep',1,'InitialStep',2e-5);
%load AP_abbrev

load('KCNQ1_mutation_allVariables_WT');
mutations_expinputs=readtable('Vanoye_TS2.xlsx');
mutation_names=table2array(mutations_expinputs(:,1));


run_time = 7e3;
len_plot=floor(run_time/2);


mt= 1; %T104I
name_test_mt1=mutation_names(mt);
filename=char(strcat('KCNQ1_mutation_allVariables_', name_test_mt1));
load(filename)
all_ICs_mt1=mut_ICs;
all_parameters_mt1=mut_parameters;
all_outputs_mt1=mut_outputs;
Time_trials_mt1=mut_Time_APtrace;
Vm_trials_mt1=mut_Vm_APtrace;

flag_mt1=flag;
color_mt1=[0 0.45 0.74];

mt= 20; %i227L
name_test_mt20=mutation_names(mt);
filename=char(strcat('KCNQ1_mutation_allVariables_', name_test_mt20));
load(filename)
all_ICs_mt20=mut_ICs;
all_parameters_mt20=mut_parameters;
all_outputs_mt20=mut_outputs;
Time_trials_mt20=mut_Time_APtrace;
Vm_trials_mt20=mut_Vm_APtrace;
flag_mt20=flag;
color_mt20=[0.85 0.33 0.1];

for j=[2649 10420]    

        y0n = wt_ICs(j,:);
        model_parameter_inputs_trial = wt_parameters(j,:);
        [t_wt,y_wt] = ode15s(@ipsc_function,[0 wt_outputs(j,10)+200],y0n,options,model_parameter_inputs_trial);
        [ Ik1_wt, Ito_wt, Ikr_wt, Iks_wt, ICaL_wt, INaK_wt, INa_wt, INaCa_wt ] = Current_outputs( t_wt, y_wt, model_parameter_inputs_trial);
        dvdt=(y_wt(2:end,1)-y_wt(1:end-1,1))./(t_wt(2:end)-t_wt(1:end-1));
        [~,i]=max(dvdt); t_offset_wt=t_wt(i);
        
        y0n = all_ICs_mt1(j,:);
        model_parameter_inputs_trial = all_parameters_mt1(j,:);
        [t_mt1,y_mt1] = ode15s(@ipsc_function,[0 all_outputs_mt1(j,10)+200],y0n,options,model_parameter_inputs_trial);
        [ Ik1_mt1, Ito_mt1, Ikr_mt1, Iks_mt1, ICaL_mt1, INaK_mt1, INa_mt1, INaCa_mt1 ] = Current_outputs( t_mt1, y_mt1, model_parameter_inputs_trial);
        dvdt=(y_mt1(2:end,1)-y_mt1(1:end-1,1))./(t_mt1(2:end)-t_mt1(1:end-1));
        [~,i]=max(dvdt); t_offset_mt1=t_mt1(i);
        
        y0n = all_ICs_mt20(j,:);
        model_parameter_inputs_trial = all_parameters_mt20(j,:);
        [t_mt20,y_mt20] = ode15s(@ipsc_function,[0 all_outputs_mt20(j,10)+200],y0n,options,model_parameter_inputs_trial);
        [ Ik1_mt20, Ito_mt20, Ikr_mt20, Iks_mt20, ICaL_mt20, INaK_mt20, INa_mt20, INaCa_mt20 ] = Current_outputs( t_mt20, y_mt20, model_parameter_inputs_trial);
        dvdt=(y_mt20(2:end,1)-y_mt20(1:end-1,1))./(t_mt20(2:end)-t_mt20(1:end-1));
        [~,i]=max(dvdt); t_offset_mt20=t_mt20(i);
        
        figure,set(gcf,'color','w')
        subplot(4,1,1)
        hold on
        set(gca,'box','off','tickdir','out','fontsize',30, 'LineWidth', 2)
        plot(t_wt-t_offset_wt,y_wt(:,1), 'k', 'LineWidth',2);
        plot(t_mt1-t_offset_mt1,y_mt1(:,1), 'Color', color_mt1, 'LineWidth',2);
        plot(t_mt20-t_offset_mt20,y_mt20(:,1), 'Color', color_mt20, 'LineWidth',2);
        ylabel('Voltage (mV)')
        xlabel('Time (s)')
        %title_hold="Mutation " + name_test_mt1 + " and "+ name_test_mt20 +", cell number "+ num2str(j) ;
        title_hold="Cell #"+ num2str(j) ;
        title(title_hold)
        legend(['Control', name_test_mt1, name_test_mt20])
        legend boxoff
        
        subplot(4,1,2)
        hold on
        set(gca,'box','off','tickdir','out','fontsize',30, 'LineWidth', 2)
        plot(t_wt-t_offset_wt,Iks_wt, 'k', 'LineWidth',2);
        plot(t_mt1-t_offset_mt1,Iks_mt1, 'Color', color_mt1, 'LineWidth',2);
        plot(t_mt20-t_offset_mt20,Iks_mt20, 'Color', color_mt20, 'LineWidth',2);
        ylabel('I_{Ks} (pA/pF)')
        xlabel('Time (s)')
        
        subplot(4,1,3)
        hold on
        set(gca,'box','off','tickdir','out','fontsize',30, 'LineWidth', 2)
        plot(t_wt-t_offset_wt,Ikr_wt, 'k', 'LineWidth',2);
        plot(t_mt1-t_offset_mt1,Ikr_mt1, 'Color', color_mt1, 'LineWidth',2);
        plot(t_mt20-t_offset_mt20,Ikr_mt20, 'Color', color_mt20, 'LineWidth',2);
        ylabel('I_{Kr} (pA/pF)')
        xlabel('Time (s)')
        hold off
        
        subplot(4,1,4)
        hold on
        set(gca,'box','off','tickdir','out','fontsize',30, 'LineWidth', 2, 'ylim', [-1 0.1])
        plot(t_wt-t_offset_wt,INa_wt, 'k', 'LineWidth',2);
        plot(t_mt1-t_offset_mt1,INa_mt1, 'Color', color_mt1, 'LineWidth',2);
        plot(t_mt20-t_offset_mt20,INa_mt20, 'Color', color_mt20, 'LineWidth',2);
        ylabel('I_{CaL} (pA/pF)')
        xlabel('Time (s)')
        hold off

end






