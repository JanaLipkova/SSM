clear; clc


eps = {'0.01', '0.03', '0.05'};

load('SSA/hist.mat');
f_ssa=frq;
e_ssa=edges;

N = size(e_ssa,1);
M = size(e_ssa,2);
%% ========================================================================

SpList = [6 13 14];

col = get(groot,'DefaultAxesColorOrder');

err_t = zeros(1,length(eps));
err_r = zeros(1,length(eps));
err_s = zeros(1,length(eps));
epsn  = zeros(1,length(eps));

for k = 1:length(eps)

    epsn(k) = str2double(eps{k}); 
    
    load([ 'AdaptiveS/eps_' eps{k} '_hist.mat']);
%     load([ 'SLeaping/eps_' eps{k} '_hist.mat']);
    f_sleap=frq;
    e_sleap=edges;

    load([ 'AdaptiveTau/eps_' eps{k} '_hist.mat']);
%     load([ 'TauLeaping/eps_' eps{k} '_hist.mat']);
    f_tleap=frq;
    e_tleap=edges;

    load([ 'RLeapingJana/eps_' eps{k} '_hist.mat']);
    f_rleap=frq;
    e_rleap=edges;
    
    
    for i = 2:N-1
        for Sp = SpList
            h = e_ssa{i,Sp}(2) - e_ssa{i,Sp}(1);
            err_t(k) = err_t(k) + h*sum( abs(f_tleap{i,Sp}-f_ssa{i,Sp}) );
            err_r(k) = err_r(k) + h*sum( abs(f_rleap{i,Sp}-f_ssa{i,Sp}) );
            err_s(k) = err_s(k) + h*sum( abs(f_sleap{i,Sp}-f_ssa{i,Sp}) );
        end
    end
    
    
    
end
M = length(SpList);
err_t = err_t/(N-2)/M;
err_r = err_r/(N-2)/M;
err_s = err_s/(N-2)/M;   

%%
fig=figure(1); clf

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex'); 

loglog(epsn,err_t,'x-'); 

hold on; grid on; axis tight

loglog(epsn,err_r,'o-')
loglog(epsn,err_s,'s-');



set(findall(fig,'-property','FontSize'),'FontSize',16)


ax = gca;
set(ax.Children, 'MarkerSize', 13)
set(ax.Children, 'LineWidth', 2)
ax.XLabel.String = '$\varepsilon$';
ax.YLabel.String = 'histogram distance';

lh = legend('$\;$ $\tau$-leap','$\;$ r-leap','$\;$ s-leap');
lh.FontSize = 18;
lh.Location='best';

