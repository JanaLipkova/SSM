clear; clc


eps = '0.05';

load('SSA/hist.mat');
f_ssa=frq;
e_ssa=edges;

load([ 'AdaptiveS/eps_' eps '_hist.mat']);
f_sleap=frq;
e_sleap=edges;

load([ 'AdaptiveTau/eps_' eps '_hist.mat']);
f_tleap=frq;
e_tleap=edges;

load([ 'RLeapingJana/eps_' eps '_hist.mat']);
f_rleap=frq;
e_rleap=edges;

Sp = 1;

N = size(e_ssa,1);

%% ========================================================================


 col = get(groot,'DefaultAxesColorOrder');



for i=2:N-1
  
    clf
    
    centers = (e_ssa{i,Sp}(1:end-1) + e_ssa{i,Sp}(2:end))/2;
    bh=bar(centers,f_ssa{i,Sp});
    bh.BarWidth=1;
    bh.FaceColor = col(1,:);
    
    hold on
    
    centers = (e_tleap{i,Sp}(1:end-1) + e_tleap{i,Sp}(2:end))/2;
    bh=bar(centers,f_tleap{i,Sp});
    bh.BarWidth=1;
    bh.FaceColor = col(2,:);
    bh.FaceAlpha = 0.9;
    
    
    centers = (e_rleap{i,Sp}(1:end-1) + e_rleap{i,Sp}(2:end))/2;
    bh=bar(centers,f_rleap{i,Sp});
    bh.BarWidth=1;
    bh.FaceColor = col(3,:);
    bh.FaceAlpha = 0.9;
    
    
    centers = (e_sleap{i,Sp}(1:end-1) + e_sleap{i,Sp}(2:end))/2;
    bh=bar(centers,f_sleap{i,Sp});
    bh.BarWidth=1;
    bh.FaceColor = col(4,:);
    bh.FaceAlpha = 0.9;
    
    legend('SSA', 'tau-leap', 'r-leap', 's-leap')
    
    
    
    pause
    
    
    
    
    
    
end
   




