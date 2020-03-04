clear
clc
close all

options = odeset('MaxStep',1,'InitialStep',2e-5);


load('KCNQ1_mutation_allVariables_WT');

mutations_etripinputs=readtable('Vanoye_TS2.xlsx');
mutation_names=table2array(mutations_etripinputs(:,1));

%threshold percentage:
percen_cut=0.04;


%%
pick_color=0;
col=hsv(5);
coeffs=zeros(length(mutation_names),2);
percent_gcutoff=zeros(size(mutation_names));
percent_exclude=zeros(size(mutation_names));

muts_to_analyze=1:length(mutation_names);
for mt= muts_to_analyze

 pick_color=pick_color+1; 
    
name_test=mutation_names(mt);
filename=char(strcat('KCNQ1_mutation_allVariables_', name_test));
load(filename)


tri= (mut_triangulation-wt_triangulation)./wt_triangulation ;
z=isnan(mut_triangulation)+wt_ab_repol';
tri=tri(z==0);
b2b=(mut_b2b-wt_b2b)./wt_b2b ;
apd=((mut_outputs(:,3)-wt_outputs(:,3))./wt_outputs(:,3)) ;

b2b=b2b(z==0);
apd=apd(z==0);

apd_thresh=apd(tri>=percen_cut);
b2b_thresh=b2b(tri>=percen_cut);
tri_thresh=tri(tri>=percen_cut);

tri_thresh=tri_thresh(apd_thresh>=percen_cut);
b2b_thresh=b2b_thresh(apd_thresh>=percen_cut);
apd_thresh=apd_thresh(apd_thresh>=percen_cut);

tri_thresh=tri_thresh(b2b_thresh>=percen_cut);
apd_thresh=apd_thresh(b2b_thresh>=percen_cut);
b2b_thresh=b2b_thresh(b2b_thresh>=percen_cut);
if isempty(tri)==0
percent_gcutoff(mt)=100*length(tri_thresh)/length(tri);
percent_exclude(mt)=100*(1-length(z(z==0))/length(wt_ab_repol(wt_ab_repol==0)));

else
   percent_gcutoff(mt)=0;
   percent_exclude(mt)=0;
end
  

    if mt==1 %0Gks
        figure,set(gcf,'color','w')
        hold on
        set(gca,'box','off','tickdir','out','fontsize',28, 'LineWidth', 2)%, 'trilim', [-50 300], 'ylim', [-50 150])
        plot(tri.*100,b2b.*100,'.', 'Color', [0.5 0.5 0.5])
        plot(tri_thresh.*100,b2b_thresh.*100,'r.')
        plot([-20 60],[percen_cut percen_cut].*100,'--k')
        plot([percen_cut percen_cut].*100,[-100 200],'--k')
        ylabel('B2B')
        xlabel('triang')
        title(name_test)
        
        figure,set(gcf,'color','w')
        hold on
        set(gca,'box','off','tickdir','out','fontsize',30, 'LineWidth', 2)%, 'trilim', [-50 300], 'ylim', [-50 150])
        plot(tri.*100,apd.*100,'.', 'Color', [0.5 0.5 0.5])
        plot(tri_thresh.*100,apd_thresh.*100,'r.')
        plot([-20 60],[percen_cut percen_cut].*100,'--k')
        plot([percen_cut percen_cut].*100,[-20 60],'--k')
        ylabel('APD')
        xlabel('triang')
        title(name_test)
    end
    
    if mt==8 %F2791
        
        figure,set(gcf,'color','w')
        hold on
        set(gca,'box','off','tickdir','out','fontsize',28, 'LineWidth', 2)%, 'trilim', [-50 300], 'ylim', [-50 150])
        plot(tri.*100,b2b.*100,'.', 'Color', [0.5 0.5 0.5])
        plot(tri_thresh.*100,b2b_thresh.*100,'.','Color', [0.93 0.63 0.13])
        plot([-20 60],[percen_cut percen_cut].*100,'--k')
        plot([percen_cut percen_cut].*100,[-100 200],'--k')
        ylabel('B2B')
        xlabel('triang')
        title(name_test)
        
        figure,set(gcf,'color','w')
        hold on
        set(gca,'box','off','tickdir','out','fontsize',30, 'LineWidth', 2)%, 'trilim', [-50 300], 'ylim', [-50 150])
        plot(tri.*100,apd.*100,'.', 'Color', [0.5 0.5 0.5])
        plot(tri_thresh.*100,apd_thresh.*100,'.','Color', [0.93 0.63 0.13])
        plot([-20 60],[percen_cut percen_cut].*100,'--k')
        plot([percen_cut percen_cut].*100,[-20 60],'--k')
        ylabel('APD')
        xlabel('triang')
        title(name_test)
    end
    
    if mt==6 %V207M
        figure,set(gcf,'color','w')
        hold on
        set(gca,'box','off','tickdir','out','fontsize',28, 'LineWidth', 2)%, 'trilim', [-50 300], 'ylim', [-50 150])
        plot(tri.*100,b2b.*100,'.', 'Color', [0.5 0.5 0.5])
        plot(tri_thresh.*100,b2b_thresh.*100,'g.')
        plot([-20 60],[percen_cut percen_cut].*100,'--k')
        plot([percen_cut percen_cut].*100,[-100 200],'--k')
        ylabel('B2B')
        xlabel('triang')
        title(name_test)
        
        figure,set(gcf,'color','w')
        hold on
        set(gca,'box','off','tickdir','out','fontsize',30, 'LineWidth', 2)%, 'trilim', [-50 300], 'ylim', [-50 150])
        plot(tri.*100,apd.*100,'.', 'Color', [0.5 0.5 0.5])
        plot(tri_thresh.*100,apd_thresh.*100,'g.')
        plot([-20 60],[percen_cut percen_cut].*100,'--k')
        plot([percen_cut percen_cut].*100,[-20 60],'--k')
        ylabel('APD')
        xlabel('triang')
        title(name_test)
    end
    

end

percent_gcutoff=percent_gcutoff(muts_to_analyze);
percent_exclude=percent_exclude(muts_to_analyze);

figure,set(gcf,'color','w')
hold on
        set(gca,'box','off','tickdir','out','fontsize',28, 'LineWidth', 2)
plot(ones(1,length(percent_gcutoff)), percent_gcutoff,'k.','Markersize', 20)
text(ones(1,length(percent_gcutoff)),percent_gcutoff,mutation_names(muts_to_analyze),'VerticalAlignment','bottom','HorizontalAlignment','right')
ylabel('% population meeting APD90, b2b, and Tri criteria')
title(percen_cut)

%%
figure,set(gcf,'color','w')
hold on
        set(gca,'box','off','tickdir','out','fontsize',28, 'LineWidth', 2)
plot(ones(1,length(percent_gcutoff)), percent_gcutoff,'k.','Markersize', 20)
plot([0 2], [4 4],'--k','LineWidth', 2)
plot([0 2], [1 1],'--k','LineWidth', 2)
%text(ones(1,length(percent_gcutoff)),percent_gcutoff,mutation_names(muts_to_analyze),'VerticalAlignment','bottom','HorizontalAlignment','right')
ylabel('% population meeting APD90, b2b, and Tri criteria')
title(percen_cut)

%% percent meeting threshold criteria
percent_gcutoff

%% percent excluded
percent_exclude

mutation_names(muts_to_analyze)

