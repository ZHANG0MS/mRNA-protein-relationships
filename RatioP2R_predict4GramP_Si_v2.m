clear
%% condition 1: use ratio of all negative
%% condition 2: use ratio of all bacteria 
str_ori = {'Si', 'Sg', 'Sa', 'Pg', 'Aa Y4', 'VT1169BY', 'VT1169TY', 'Aa 624', 'Pa'};
gramN = [9, 8,5, 6:7,4,2,3,1];
mat_mdl_ori = [   
    0.5126    0.9437    0.6050    0.9329    0.7721    0.9083    0.9537    0.7736    0.7016
    8.2005    2.4294    6.6596    2.9145    4.3209    2.7567    2.4671    4.6491    5.0745];
mat_mdl = mat_mdl_ori(:,gramN);
Str_gramN = {'Pa','Aa 624','Aa Y4','VT1169BY','Vt1169TY','Pg','Sg','Sa','Si'};


id_species=8;
% Gram+
if id_species<9
load('Species=6_table_Num.mat');
tag1 = table_allNum.PA_tag;
PTR_one = mean(table_allNum.mat_ortholog_PTR,2);
load('Species=8_table_Num.mat');
tag2 = table_allNum.PA_tag;
end
% archaea
if id_species>8
load('Species=8_table_Num.mat');
tag1 = table_allNum.PA_tag;
PTR_one = mean(table_allNum.mat_ortholog_PTR,2);
load('Species=9_table_Num.mat');
tag2 = table_allNum.PA_tag;
end
%% plot for which bacteria species 
% indicate i

RNA_Si = table_allNum.mat_ortholog_RNA(:,id_species);
Protein_Si = table_allNum.mat_ortholog_protein(:,id_species);

[~,ia, ib] = intersect(tag1, tag2);
PTR_Si = PTR_one(ia);
XX_Si = RNA_Si(ib);
YY_Si = Protein_Si(ib);

YY_c = XX_Si + PTR_Si;
XX_c = (XX_Si + PTR_Si-mat_mdl(2,id_species))./mat_mdl(1,id_species);%
[array_after,~]=outputcoeff_mat2(YY_Si, YY_c, 1);
[array_before,~]=outputcoeff_mat2(XX_Si,YY_Si, 1);
[array_PTR,~]=outputcoeff_mat2(YY_Si-XX_Si,PTR_Si,1);
sprintf(['The correlation for archaea is ',num2str(array_before),'before and ', num2str(array_after),'after correction'])

% figure(1),subplot(1,2,1),hold on, title("Protein"), scatter(YY_Si,YY_c,10,'r');axis square
% subplot(1,2,2), title("PTR"),scatter(YY_Si-XX_Si,PTR_Si,10,'b');axis square
% hold on
for idx = 1:2
    addpath('D:/Programs/2020Matlab/PlotPub-master/lib')
    plt=Plot([100 100],[100 100]);
    plt.BoxDim = [2, 2];
    hold on, axis square , title(Str_gramN(id_species))
    if idx==1, scatter(XX_Si, YY_Si,20,'k');ylim([4 20]);
        xlabel('log_2(mRNA)'),ylabel('log_2(observed protein)');
        legend({['N=',num2str(length(YY_c))],['\rho_{before}=',num2str(round(100*array_before)/100)]} ,'Location','southeast','Box','off') ;
    end
    if idx==2, scatter(YY_c,YY_Si,20,'r');ylim([4 20]);xlim([4 20]);
        ylabel('log_2(observed protein)'),xlabel('log_2(predicted protein)');
        legend({['N=',num2str(length(YY_c))],['\rho_{after}=',num2str(round(100*array_after)/100)]} ,'Location','southeast','Box','off') ;
    end
    saveas(gcf,['Scatter_UseBacteria_PredictFor= ',Str_gramN{id_species},'PLOT',num2str(idx),'.pdf'])  
end
%%

%% function
function [array_coeff,array_pval]=outputcoeff_mat2(protein, rna, type)
% protein=esemblepro';
% rna=esemblerna';
% type=0;
%
sz= size(protein);
array_coeff=zeros(1, sz(2));
array_pval=zeros(1, sz(2));
    for i=1:sz(2)
        temp_pro=protein(:,i);
        temp_rna=rna(:,i);
        if type==1
            [array_coeff(i),array_pval(i)]=corr( temp_pro, temp_rna,'Type','Spearman');
        end
        if type==0
            [array_coeff(i),array_pval(i)]=corr( temp_pro, temp_rna,'Type','Pearson');
        end
    end
end