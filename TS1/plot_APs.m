
clear all
close all
clc

%% Output selection
load KCNQ1_mutation_allVariables_WT
[N_trials, ~]=size(wt_ICs);

cells_to_plot=round(N_trials./6);

mutations_expinputs=readtable('Vanoye_TS1.xlsx');
mutation_names=table2array(mutations_expinputs(:,1));

select_muts=[1 5 6];

[~, N_outputs] = size(wt_outputs);

for mt=select_muts
    name_test=mutation_names(mt);
    filename=char(strcat('KCNQ1_mutation_allVariables_', name_test));
    load(filename)
    
     figure,set(gcf,'color','w')
    hold on
    
    for j= 1:cells_to_plot
        v_mut=mut_Vm_APtrace(:,j);
        if sum(mut_outputs(j,:))~=0  && (max(v_mut)-min(v_mut))>70  && mut_outputs(j,3)< mut_outputs(j,8)+320 %&& ab_repol_noise(j)==0   
        
        t_mut=mut_Time_APtrace(:,j);
        v_mut=mut_Vm_APtrace(:,j);
        v_mut=v_mut(t_mut~=0); t_mut=t_mut(t_mut~=0);
        v_mut=v_mut(1:end-10); t_mut=t_mut(1:end-10);
        
        t_wt=wt_Time_APtrace(:,j);
        v_wt=wt_Vm_APtrace(:,j);
        v_wt=v_wt(t_wt~=0);t_wt=t_wt(t_wt~=0);
         v_wt=v_wt(1:end-10); t_wt=t_wt(1:end-10);
        
        if mt==select_muts(1)
            color_select=[0.64 0.08 0.18];
        elseif mt==select_muts(2)
            color_select=[0.85 0.33 0.1];
        else
            color_select=[0.49 0.18 0.56];
        end
        
        
        plot(t_mut,v_mut, 'Color', color_select, 'LineWidth',1.5);
        plot(t_wt,v_wt, 'Color', [0.5 0.5 0.5], 'LineWidth',1.5);
        set(gca,'box','off','tickdir','out','fontsize',30, 'LineWidth', 2, 'xlim', [-100 2000], 'ylim', [-85 50])
        
        end
    end
    
    title(name_test)
    ylabel('Voltage (mV)')
    xlabel('Time (s)')
    %hold off
    name_test
  
end

