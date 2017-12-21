%============================================
%
% Plot results of Decaying Dimerization
%
%============================================
% Read in file with histogram from SSM and visualise it



% function plotDimerizationALL

% 0) load data
bSave=1;
myBase    = '../../Simulations/Dimerization/NonStiff/';


for eps = [0.05, 0.03, 0.01]
    
    inputSSA  = ['SSA/DimerizationSSA_histogram.txt'];
    dataSSA   = load([myBase,inputSSA]);
    
    inputR  = ['RLeapingJana/eps_',num2str(eps),'/Dimerization-RLeapingJana_histogram.txt'];
    dataR   = load([myBase,inputR]);
    
    inputT  = ['TauLeaping/eps_',num2str(eps),'/Dimerization-TauLeaping_histogram.txt'];
    dataT = load([myBase,inputT]);
    
    inputS  = ['SLeaping/eps_',num2str(eps),'/Dimerization-SLeaping_histogram.txt'];
    dataS = load([myBase,inputS]);
    
    % select species to be plotted
    dataSSA = dataSSA(:,2);
    dataR   = dataR(:,2);
    dataT   = dataT(:,2);
    dataS   = dataS(:,2);
    
    % get optimal number of bins (data must be a row vector)
    minNbins = 1;   % range of
    maxNbins = 100;
    optNbins = optBINS(dataSSA',minNbins,maxNbins)
    
    % histogram set up
    bNomralisation  = 1;
    bBarPlot        = 1;      % plot histogram as bar plot, else as markers
    marker          = 'o';
    
    % plotting setups
    cmap    = get(0, 'defaultaxescolororder');
    fs      = 20;    % font size
    lw      = 3;     % line width
    ms      = 13;    % marker size
    ga      = 0.4;   % grid opacity (alpha)
    opacity = 0.3;   % bar opacity
    
    figure
    hold on
    
    plotHistogramSSM(dataSSA,optNbins, 1 , cmap(5,:), opacity, lw-1);
    plotHistogramSSM(dataR,optNbins, 0, cmap(1,:),  1,         lw, 'o', ms);
    plotHistogramSSM(dataT,optNbins, 0, cmap(4,:),  1,         lw, 'x', ms);
    plotHistogramSSM(dataS,optNbins, 0, cmap(2,:),  1,         lw, '*', ms);
    legend('SSA','R','T','S')
    title(['eps=',num2str(eps)])
    set(gca,'FontSize',fs);
    set(gca,'GridLineStyle', '--','LineWidth',lw);
    set(gca,'Box','on');
    set(gca,'LineWidth',lw);
    % set(gca,'YGrid','on');
    % set(gca,'GridAlpha',ga);
    
    set(gcf,'pos',[300 300 700 400])
    set(gcf,'papersize',[10,10]);
    
    if(bSave)
        print(gcf,[myBase,'Histograms_eps_',num2str(eps),'.pdf'],'-dpdf');
    end;
    
end;


