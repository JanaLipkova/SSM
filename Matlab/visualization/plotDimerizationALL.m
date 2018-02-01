%============================================
%
% Plot results of Decaying Dimerization
%
%============================================
% Read in file with histogram from SSM and visualise it



% function plotDimerizationALL

addpath('../lib/tightfig')
addpath('../lib/')


set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

% 0) load data
bSave=1;
myBase    = '../../../Simulations/Dimerization/NonStiff/';
im=1;

% figure
for eps = [0.05, 0.03, 0.01]
    
    inputSSA  = ['SSA/Dimerization_SSA_histogram.txt'];
    dataSSA   = load([myBase,inputSSA]);
    
    inputR  = ['RLeapingJana/eps_',num2str(eps),'/Dimerization_RLeapingJana_histogram.txt'];
    dataR   = load([myBase,inputR]);
    
%     inputT  = ['SLeaping/eps_',num2str(eps),'/Dimerization_SLeaping_histogram1.txt'];
    inputT  = ['AdaptiveTau/eps_',num2str(eps),'/Dimerization_AdaptiveTau_histogram.txt'];
    dataT = load([myBase,inputT]);
    
    inputS  = ['AdaptiveS/eps_',num2str(eps),'/Dimerization_AdaptiveS_histogram1.txt'];
    %      inputS  = ['SLeaping/eps_',num2str(eps),'/Dimerization_SLeaping_histogram.txt'];
    dataS = load([myBase,inputS]);
    
    % select species to be plotted
    id=2;
    dataSSA = dataSSA(:,id);
    dataR   = dataR(:,id);
    dataT   = dataT(:,id);
    dataS   = dataS(:,id);
    
    % get optimal number of bins (data must be a row vector)
    minNbins = 1;   % range of
    maxNbins = 100;
    optNbins = optBINS(dataSSA',minNbins,maxNbins)
    optNbins=20;
    
    % histogram set up
    bNomralisation  = 1;
    bBarPlot        = 1;      % plot histogram as bar plot, else as markers
    marker          = 'o';
    
    % plotting setups
    cmap    = get(0, 'defaultaxescolororder');
    fs      = 25;    % font size
    lw      = 3;     % line width
    ms      = 13;    % marker size
    ga      = 0.4;   % grid opacity (alpha)
    opacity = 0.5;   % bar opacity
    
    fig=figure; clf
    % subplot(2,2,im)
    hold on
    plotHistogramSSM(dataSSA,optNbins, 1 , cmap(6,:), opacity, lw-1);
    plotHistogramSSM(dataT,optNbins, 0, [1,0,0],  1,         lw, 'x', ms);
    plotHistogramSSM(dataR,optNbins, 0, cmap(3,:),  1,         lw, 'o', ms);
    plotHistogramSSM(dataS,optNbins, 0, [0,0,0],  1,         lw, 's ', ms);
    legend('SSA','$\tau$-leap','r-leap','s-leap')
    title(['$\varepsilon$=',num2str(eps)])
    set(gca,'FontSize',fs);
    set(gca,'GridLineStyle', '--','LineWidth',lw);
    set(gca,'Box','on');
    set(gca,'LineWidth',lw);
    % set(gca,'YGrid','on');
    % set(gca,'GridAlpha',ga);
    xlabel('$S_2(10)$')
    im=im+1;
    
%     set(gcf,'papersize',[7.3,6.3]);
    set(gcf,'pos',[300 300 600 450])

    axis tight
    if(bSave)
%         print(gcf,[myBase,'non-stiff-dim-histogram-conv-',num2str(eps),'.pdf'],'-dpdf');
%       saveas(gcf,[myBase,'non-stiff-dim-histogram-conv-',num2str(eps),'.eps'],'epsc');
      print(gcf,[myBase,'non-stiff-dim-histogram-conv1-',num2str(eps),'.eps'],'-depsc');

    end;
    %       saveas(gcf,'time_lac_dim', 'epsc' );


    
end;



