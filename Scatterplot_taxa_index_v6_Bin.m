
clearvars -except RNA_totalcount Protein_totalcount lim_max lim_min temp_Protein_totalcount R_square mat_b1 mat_b0

densityplot=0;
essentialplot=1;    
varianceplot=0; 
barplot=0;
residualplot=0; % set densityplot=1
shapeplot=0;
gibbsplot=0;
tigrfamplot=0;

% % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % data load  % % % % % % % % %  
% % % % % % % % % % % % % % % % % % % % % % % % % 
taxa = 3
% must bigger than zero
if taxa ==1, str_taxa = 'Si';end
if taxa ==2, str_taxa = 'Sg';end
if taxa ==3, str_taxa = 'LAC';end
if taxa ==4, str_taxa = 'Pg';end
if taxa ==5, str_taxa = 'Y4';end
if taxa ==6, str_taxa = 'BY';end
if taxa ==7, str_taxa = 'TY';end
if taxa ==8, str_taxa = 'Aa624';end
if taxa ==9, str_taxa = 'PA14';end
if taxa ==10, str_taxa = 'PA14hr25';end % check protein ID is correct or NOT
if taxa<5
[table_RNA_Protein, table_RNA]=load_RNA_protein_mapped(str_taxa);
whos
end
if taxa>4
[table_RNA_Protein, table_RNA]=load_RNA_protein_mapped_AA(str_taxa);
whos
end

% % % % % % deal with RNA hotspot % % % % % % % % % 
id_hotspot = 2;
if id_hotspot==1 % % remove genes with 0.1 M reads genes , detected 3 out of 4 replicates (we don't want consistency)
    id_hotspot = table_RNA_Protein.RNA>1E5;
    id_hotspot_NOTdetected = sum(id_hotspot,2)<3;
end
if id_hotspot==2 % % remove genes with average more than 0.1 M reads genes
    id_hotspot = mean(table_RNA_Protein.RNA,2)>1E5;
    id_hotspot_NOTdetected = id_hotspot<1;
end
table_RNA_Protein_ORI = table_RNA_Protein;
table_RNA_Protein=table_RNA_Protein(id_hotspot_NOTdetected,:);
disp(['There are ', num2str(sum(id_hotspot_NOTdetected==0)),' GENE with RNA > 0.1M detected in 3 out of 4 replicates'])

% % with RNA and protein mapped
RNA = table_RNA_Protein.RNA;
RNA_all = table_RNA.RNA; 
Protein = table_RNA_Protein.Protein;
locustag= table_RNA_Protein.locustag;

% % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % data normalization  % % % % % % % % %  
% % % % % % % % % % % % % % % % % % % % % % % % % 
%% scale across replicates all counts to same depth, RNA ONLY
% scale RNA to 1 million;
RNA_size=size(RNA);
Protein_size=size(Protein);
RNA_proportional=zeros(RNA_size);
if exist('RNA_totalcount','var')==0,RNA_totalcount = zeros(10,4);end
if exist('Protein_totalcount','var')==0,Protein_totalcount = zeros(10,4);end
% RNA_totalcount = zeros(10,4);
% Protein_totalcount = zeros(10,4);
RNA_totalcount(taxa,1:RNA_size(2))=sum(RNA_all)';% normalized based on the all proteins
Protein_totalcount(taxa,1:Protein_size(2))=sum(Protein)';
for i=1:RNA_size(2)
    RNA_proportional(:,i)=log2(1+RNA(:,i)*1e6/RNA_totalcount(taxa,i));
end
%% Scale Protein to median of the total amount = 10M
Protein_proportional=zeros(Protein_size);
for j=1:Protein_size(2)
    Protein_proportional(:,j)=log2( 1+Protein(:,j)*1e7/Protein_totalcount(taxa,1) );% based on the first replicates
    temp_Protein(:,j) = Protein(:,j)*1e7/Protein_totalcount(taxa,1);
end 
Protein_totalcount(taxa,1)
%% total RNA after normalization (test)
if exist('temp_Protein_totalcount','var')==0,temp_Protein_totalcount = zeros(10,4);end
temp_Protein_totalcount(taxa,1:Protein_size(2)) = sum(temp_Protein)';

%% get average for each growth rate
% for mRNA
RNA_mean_logcounts=mean(RNA_proportional,2);

% for protein
Protein_mean_logcounts=mean(Protein_proportional,2);

%% OVERALL CORRELATION
% using Spearman's rank correlation or Pearson's correlation
% INPUT: x, y, id_correlation
% OUTPUT: correlation coefficent
% 0=Pearson, 1=Rank, Spearman
[array_coeff1(taxa),array_pval1(taxa)]=outputcoeff_mat2(Protein_mean_logcounts, RNA_mean_logcounts, 1)

%% Overall Regression
md1 = fitlm(RNA_mean_logcounts, Protein_mean_logcounts);
mat_b1(taxa)= md1.Coefficients.Estimate(2);
mat_b0(taxa)= md1.Coefficients.Estimate(1);
p_valF(taxa) = coefTest(md1);
R_square(taxa) = md1.Rsquared.Ordinary;
md_overall = md1; 



% % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % scatter & density scatter % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % 
%% OVERALL REGRESSION
%% OVERALL PICUTRE
%% scatter plot
%% overall correlation of mRNA and protein (included all growth rates)
%% density plot=overall correlation of mRNA and protein (included all growth rates)
blocksz=20;plotid=1;
if densityplot==1
    close all
    tempRR = RNA_mean_logcounts;
    tempPP = Protein_mean_logcounts;
    %% plot for the density scatter plot
    % % % % % % density scatter % % % % % % 
    [values_all, centers_all] = hist3([tempRR,tempPP],[blocksz blocksz]);
    % % % % % colorbar setup
    maxN = max(max(values_all));
    minN = 7;% threshold 1
    cdata = [...
        minN   0 0 0  0  255  255   255
           0 255 255 255 maxN   255   0   0];
    if plotid==1
        addpath('D:\Programs\2020Matlab\Plot_cptcmap');
        dlmwrite('D:\Programs\2020Matlab\Plot_cptcmap\mycmap.cpt', cdata, ' ');
        axes;
        imagesc(centers_all{:}, values_all.');
        set(gca,'YDir','normal')
        cptcmap('D:\Programs\2020Matlab\Plot_cptcmap\mycmap', 'mapping', 'direct'); 
        c = colorbar;c.LineWidth = 1; c.TickLength = 0.01; c.TickDirection='out';
    %         title('all growth rates');
        set(gca,'FontSize',20),axis square,
        ax=gca;ax.TickDir = 'out';ax.TickLength = [0.02 0.035];%ax.XTick=xinter(1:labelstep(1):end);ax.YTick=Time0(1:labelstep(2):end)+4;
        ax.LineWidth = 1;
        xlabel('log_2mRNA','FontSize',25),ylabel('log_2Protein','FontSize',25);
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.9]);
        % set(gcf,'PaperOrientation','portrait')
        print(['density_scatterplot', str_taxa,'.pdf'],'-dpdf','-r600');
    end
    
    % % % % % % scatter plot % % % % % %  
    addpath('D:/Programs/2020Matlab/PlotPub-master/lib')
    % start to plot
    plt=Plot([-2 -1],[-1 -1]);
    % size of plots
    plt.BoxDim = [2.5, 2.5]; % for scatter plot, [width, height] in inches
    % plt.XLim = [0 Inf]; % [min, max]
    % plt.YLim = [0 Inf]; %
    plt.Colors = {                 % three colors for three data set
    [0,   0,    0]     % data set 2
    [0,      0,       0]        % data set 1
    [0,      0,       1]        % data set 3
    [0,      1,       0] 
    };
    % plot with GRAY or BLACK
    %hold on, scatter(tempRR,tempPP,'k','MarkerEdgeAlpha',.1);
    hold on, scatter(tempRR,tempPP,20,'k');

    %ax=gca;ax.TickDir = 'out'
    xlabel('log_2(mRNA)'),ylabel('log_2(Protein)') 
    % title(['rpm=',num2str(2^(i-1))]),pbaspect([1 1 1])
    legend({['\rho=',num2str(round(100*array_coeff1(taxa))/100)],[ 'N=',num2str(length(locustag))]} ,'Location','southeast','Box','off') ;

    lim_max(taxa,1)=max(centers_all{1,1})+( max(centers_all{1,1})-min(centers_all{1,1}))/(blocksz-1);
    lim_max(taxa,2)=max(centers_all{1,2})+( max(centers_all{1,2})-min(centers_all{1,2}))/(blocksz-1);
    
    lim_min(taxa,1)=min(centers_all{1,1})-( max(centers_all{1,1})-min(centers_all{1,1}))/(blocksz-1);
    lim_min(taxa,2)=min(centers_all{1,2})-( max(centers_all{1,2})-min(centers_all{1,2}))/(blocksz-1);
    % lower for min and higher for max
    xy_max = ceil(lim_max(taxa,:));
    xy_min = floor(abs(lim_min(taxa,:)));
%     lim_max = ceil(lim_max);
    xlim([xy_min(1) xy_max(1)]);
    ylim([xy_min(2) xy_max(2)]);
    plt.YTick = [0:5:20]; 
    plt.XTick = [0:5:20]; 
    title(str_taxa);
    saveas(gcf,['scatterplot_', str_taxa,'.pdf'])
end

lim_min =[
     1     4
     0     0
     0     2
     0     1
     0     1
     2     2
     2     3
     0     1
     0     2
     ];
lim_max =[
    14    18
    15    18
    15    18
    14    18
    14    20
    16    21
    16    21
    16    20
    13    18
    16    18
    ];
    
% % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % essential gene analysis % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % %
if essentialplot==1


    % essential gene by GEPTOP method
    currentFolder = pwd;
    load([pwd,'\datafiles06\output_EssentialGene_ortholog_',str_taxa,'.mat']); 

    [~,ia,~] = intersect(locustag, locustagE);%S3=434; overlapped=404
    whos ia
    disp(['There are ', num2str(length(locustagE)),'essential gene present, with', num2str(length(ia)),'detected with both RNA and protein'])
    % calculate the correlation in certain ways
    yy=Protein_mean_logcounts(ia);
    xx=RNA_mean_logcounts(ia);
    locustagR_essential=locustag(ia);
    ED1 = [0.3 0.20];
    % get nonessential part
    ssRR=size(Protein_mean_logcounts); 
    ia_not = setdiff(1:ssRR(1),ia);
    yn=Protein_mean_logcounts(ia_not);
    xn=RNA_mean_logcounts(ia_not);
    ssxx=size(xx);
    % get the xlabel ready
    xinter=0:1:20;
    xinter_idx=[lim_max(taxa,1)+1  lim_max(taxa,2)+1];
    xinter_dt = [1,1];
    % plot for the distribution of essential VS all genes
        xxyy=[RNA_mean_logcounts, Protein_mean_logcounts];
        xeye=[xx, yy];
        xnyn=[xn, yn];
        ddplot=1;% if tt=1, loop into the distribution 
        while (ddplot>0)&&(ddplot<3)
        dd=1;distribution=['pdf';'cdf'];
        % start to plot things
        addpath('D:/Programs/2020Matlab/PlotPub-master/lib')
        plt=Plot([-2 -1],[-1 -1]);
        % size of plots
        plt.BoxDim = [2, 1.5]; % for scatter plot, [width, height] in inches
        % plt.XLim = [0 Inf]; % [min, max]
        % plt.YLim = [0 Inf]; %
        plt.Colors = {                 % three colors for three data set
        [0,   0,    0]     % data set 2
        [0,      0,       0]        % data set 1
        [0,      0,       1]        % data set 3
        [0,      1,       0] 
        };
        hold on, 
        % stairs plot
        histogram(xeye(:,ddplot), 0:xinter_dt(ddplot):xinter_idx(ddplot), 'LineWidth' ,2.5,'FaceAlpha',0,'EdgeColor',[0.8,   0,   0],'DisplayStyle','stairs',"Normalization",'probability'); %16=round(max(max(RNA_mean_logcounts)))
        histogram(xnyn(:,ddplot), 0:xinter_dt(ddplot):xinter_idx(ddplot), 'LineWidth' ,2.5,'FaceAlpha',0,'EdgeColor',[0,   0,   0.6],'DisplayStyle','stairs',"Normalization",'probability');
        plot([-1 25], [0 0],'k', 'LineWidth', 2.5)
        hold on, xlim([lim_min(taxa,ddplot) xinter_idx(ddplot)])
        plt.YTick = [0:0.05:ED1(ddplot)]; 
        ax=gca;ax.TickDir = 'out';box off
        %legend({"  ","Essential","Non-essential"},'Box','off','Location','northwest' )
        ylabel('frequency');
        if ddplot==1,  xlabel('log_2(mRNA)');ylim([0 ED1(ddplot)]);end
        if ddplot==2,  xlabel('log_2(Protein)');ylim([0 ED1(ddplot)]);end
%         a = get(gca,'XTickLabel');
%         set(gca,'XTickLabel',a,'fontsize',15)
        saveas(gcf,['Essential_',num2str(ddplot),'type_histogram', str_taxa,'.pdf'])
        ddplot=ddplot+1;% tt=2 for protein
        end
       %% calculate the correlations
        if ssxx(1)>0
            [array_coeffW,array_pvalW]=outputcoeff_mat2([xx;xn], [yy; yn],0);
            [array_coeffNE,array_pvalNE]=outputcoeff_mat2(yn, xn, 1);
            [array_coeffE,array_pvalE]=outputcoeff_mat2(yy, xx, 1); % spearman correlations=1
            mat_coeffE(1,:)=array_coeffE;
            mat_pvalE(1,:)=array_pvalE;
            mat_coeffE(2,:)=array_coeffNE;
            mat_pvalE(2,:)=array_pvalNE;
        end  
        
        if ssxx(1)>0
            %% Overall Regression
            ee = 1; 
            md1 = fitlm(xx, yy);
            mat_b1(ee)= md1.Coefficients.Estimate(2);
            mat_b0(ee)= md1.Coefficients.Estimate(1);
            p_valF(ee) = coefTest(md1);
            R_square(ee) = md1.Rsquared.Ordinary;
            mdE = md1;
            ee = 2; 
            md1 = fitlm(xn, yn);
            mat_b1(ee)= md1.Coefficients.Estimate(2);
            mat_b0(ee)= md1.Coefficients.Estimate(1);
            p_valF(ee) = coefTest(md1);
            R_square(ee) = md1.Rsquared.Ordinary;
            mdNE = md1;
        end
        
       %% some statistics
        aa=[median(xnyn), mean(xnyn), var(xnyn)];    bb=[median(xeye), mean(xeye),var(xeye)];
        % mRNA
        [p(1),h(1)] = ranksum(xeye(:,1),xxyy(:,1),'Tail','right'); % right: x>y 
        % Protein
        [p(2),h(2)] = ranksum(xeye(:,2),xxyy(:,2),'Tail','right'); % right: x>y
               
       %% plot for scatter plot 
        % start to plot
        plt=Plot([-2 -1],[-1 -1],[-3 -1],[-3 -1]);
        % size of plots
        sf = 1.3; 
        %plt.BoxDim = [2.5*sf, 2.5*sf]; % for scatter plot, [width, height] in inches
        plt.BoxDim = [2, 1.5]; % for scatter plot, [width, height] in inches
        plt.Colors = {                 % three colors for three data set
        [0.8,   0,    0]     % data set 2
        [0,      0,       0.6]        % data set 1
        [0,      0,       1]        % data set 3
        [0,      1,       0] 
        };
    sf=0.5;
        %hold on, scatter(xxyy(:,1),xxyy(:,2),'k');
        hold on, scatter(xnyn(:,1),xnyn(:,2),20*sf,[0,   0,   0.6]);
        scatter(xeye(:,1),xeye(:,2),20*sf,[0.8,0,0],'filled');
        %ax=gca;ax.TickDir = 'out'
        xlabel('log_2(mRNA)'),ylabel('log_2(Protein)') 
        % title(['rpm=',num2str(2^(i-1))]),pbaspect([1 1 1])
        legend({['\rho_E=',num2str(round(100*mat_coeffE(1))/100)], ['\rho_{NE}=',num2str(round(100*mat_coeffE(2))/100)], ['N=',num2str( length(ia))]},'Location','southeast','Box','off') ;
%         xlim([0 lim_max(taxa,1)])
%         ylim([0 lim_max(taxa,2)])
         xlim([0 Inf]); ylim([0 Inf])
        saveas(gcf,['essential_scatterplot', str_taxa,'.pdf'])
end


% % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % variance plot analysis % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % %
if varianceplot==1
%     RNA_var=log2(var(RNA_proportional,0,2));
%     Protein_var=log2(var(Protein_proportional,0,2));
    boxplot=1;
    % input of essential gene list
    currentFolder = pwd;
    load([pwd,'\datafiles06\output_EssentialGene_ortholog_',str_taxa,'.mat']); 
    [~,ia,~] = intersect(locustag, locustagE);%S3=434; overlapped=404
    whos ia
    % calculate the correlation in certain ways
    yy=Protein_proportional(ia,:);
    xx=RNA_proportional(ia,:);
    locustag_essential=locustag(ia);
    % get nonessential part
    ssRR=size(Protein_mean_logcounts); 
    ia_not = setdiff(1:ssRR(1),ia);
    yn=Protein_proportional(ia_not,:);
    xn=RNA_proportional(ia_not,:);
    ssxx=size(xx);
    % calculate the essential gene for differe growth rates
    xx_var=std(xx,0,2);
    yy_var=std(yy,0,2);
    xn_var=std(xn,0,2);
    yn_var=std(yn,0,2);
    % get some statistics
    all4_mean = [mean(xx_var), mean(xn_var),mean(yy_var), mean(yn_var)]; 
    all4_sem = [std(xx_var), std(xn_var),std(yy_var), std(yn_var)]./[length(xx_var), length(xn_var), length(yy_var), length(yn_var)]; 
    all4_var = [std(xx_var), std(xn_var),std(yy_var), std(yn_var)];
    %% some statistics
    addpath('D:\Programs\2020Matlab\Statistical_test')
    [px, hx] = ranksum(xn_var,xx_var,'Tail','right'); % right: x>y
    [py, hy] = ranksum(yn_var,yy_var,'Tail','right'); % right: x>y
    
    if taxa == 9 
        yinter=0:0.05:0.75;
        xinter=0:0.1:1.0;
    end
    if taxa == 3
        yinter=0:0.1:1.5;
        xinter=0:0.1:1;
    end
%% move the outlier to the biggest bins
        % essential
        id_big = yy_var>yinter(end-1);
        yy_var(id_big) = (yinter(end) + yinter(end-1))/2;
        % Non-essential
        id_big = yn_var>yinter(end-1);
        yn_var(id_big) = (yinter(end) + yinter(end-1))/2;
% move the outlier to the biggest bins
        % essential
        id_big = xx_var>xinter(end-1);
        xx_var(id_big) = (xinter(end) + xinter(end-1))/2;
        % Non-essential
        id_big = xn_var>xinter(end-1);
        xn_var(id_big) = (xinter(end) + xinter(end-1))/2;
        
 %% plot for the histogram  
    for ddplot=2:-1:1
    % stairs plot
    if ddplot==2
        % start to plot things
        distribution=['pdf';'cdf'];dd=1;
        addpath('D:/Programs/2020Matlab/PlotPub-master/lib')
        plt=Plot([-2 -1],[-1 -1]);
        % size of plots
        plt.BoxDim = [2, 1.5]; % for scatter plot, [width, height] in inches
        plt.Colors = {                 % three colors for three data set
        [0,   0,    0]     % data set 2
        [0,      0,       0]        % data set 1
        [0,      0,       1]        % data set 3
        [0,      1,       0] 
        };
         
        hold on, 
        histogram( yn_var,yinter, 'LineWidth' ,2.5,'FaceAlpha',0,'EdgeColor',[0,   0,   0.6],'DisplayStyle','stairs',"Normalization",'probability');
        histogram( yy_var,yinter, 'LineWidth' ,2.5,'FaceAlpha',0,'EdgeColor',[0.8,   0,   0],'DisplayStyle','stairs',"Normalization",'probability'); %16=round(max(max(RNA_mean_logcounts)))
%         histogram(log(yy_var),yinter, 'LineWidth' ,2.5,'FaceAlpha',0,'EdgeColor',[0.8,   0,   0],'DisplayStyle','stairs',"Normalization",distribution(dd,:)); %16=round(max(max(RNA_mean_logcounts)))
%         histogram(log(yn_var),yinter, 'LineWidth' ,2.5,'FaceAlpha',0,'EdgeColor',[0,   0,   0.6],"Normalization",distribution(dd,:),'DisplayStyle','stairs');
       plot([yinter(1), yinter(end)], [0 0],'k', 'LineWidth', 2)%title(str_taxa)
       hold on, xlabel('Standard deviation of log_2(Protein)');
       plt.XTick = [0:0.2:yinter(end)]; 
        ax=gca;ax.TickDir = 'out';
        box off,ylim([0 0.51])
        plt.YTick = [0:0.1:1.0];
        %legend({"  ","Non-essential", "Essential"},'Box','off','Location','northeast' ), 
        ylabel('frequency');
    end
    if ddplot==1
            % start to plot things
        distribution=['pdf';'cdf'];dd=1;
        plt=Plot([-2 -1],[-1 -1]);
        % size of plots
        plt.BoxDim = [2, 1.5]; % for scatter plot, [width, height] in inches
        plt.Colors = {                 % three colors for three data set
        [0,   0,    0]     % data set 2
        [0,      0,       0]        % data set 1
        [0,      0,       1]        % data set 3
        [0,      1,       0] 
        };
        hold on, 
        
        histogram(abs(xn_var),xinter, 'LineWidth' ,2.5,'FaceAlpha',0,'EdgeColor',[0,   0,   0.6],'DisplayStyle','stairs',"Normalization",'probability');       
        histogram(abs(xx_var),xinter, 'LineWidth' ,2.5,'FaceAlpha',0,'EdgeColor',[0.8,   0,   0],'DisplayStyle','stairs',"Normalization",'probability'); 

        plot([xinter(1), xinter(end)], [0 0],'k', 'LineWidth', 2)
        hold on, xlabel('Standard deviation of log_2(mRNA)'); 
         plt.XTick = [0:0.2:1.6]; 
        ax=gca;ax.TickDir = 'out';box off;
        plt.YTick = [0:0.1:0.9]; 
        ylim([0 0.305]);
        %legend({"  ","Non-essential","Essential"},'Box','off','Location','northeast' ), 
        ylabel('frequency');
    end
      saveas(gcf,[num2str(ddplot),'_std_4replicates', str_taxa,'.pdf'])
    end
        
    %% boxplot for the variance
    if boxplot==1
        Xbox=[xx_var;xn_var];
        Ybox=[yy_var;yn_var];
        g1 = repmat({'Essential'},size(xx_var));g2 = repmat({'Non-essential'},size(xn_var));
        Gbox=[g1;g2];
        clear boxplot
        figure,subplot(1,2,1),boxplot(Xbox,Gbox);
        subplot(1,2,2),boxplot(Ybox,Gbox);
        ax = gca;ax.YAxis.Scale ="log";
    end
end


% % % % % % % % % % % % % % % % % % % %  
% % % % % % FUNCTION  % % % % % % % % %  
% % % % % % % % % % % % % % % % % % % % 
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
