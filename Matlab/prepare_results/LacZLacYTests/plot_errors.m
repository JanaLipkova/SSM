clear; clc


eps = {'0.01', '0.03', '0.05'};

load('SSA/mean.mat');
f_ssa=mn;

N = size(mn,1);
M = size(mn,2);

SpList = [6 13 14];
% SpList = 6;

%% ========================================================================


col = get(groot,'DefaultAxesColorOrder');

err_t = zeros(1,length(eps));
err_r = zeros(1,length(eps));
err_s = zeros(1,length(eps));
eps_n = zeros(1,length(eps));

 for k = 1:length(eps)

    eps_n(k) = str2double(eps{k});
    
    load([ 'SLeaping/eps_' eps{k} '_mean.mat']);
    f_sleap=mn;
  
    load([ 'AdaptiveTau/eps_' eps{k} '_mean.mat']);
    f_tleap=mn;
  
    load([ 'RLeapingJana/eps_' eps{k} '_mean.mat']);
    f_rleap=mn;
    
    
    for i = 2:N-1
        for Sp = SpList
            err_t(k) = err_t(k) + abs(f_tleap{i,Sp}-f_ssa{i,Sp});
            err_r(k) = err_r(k) + abs(f_rleap{i,Sp}-f_ssa{i,Sp});
            err_s(k) = err_s(k) + abs(f_sleap{i,Sp}-f_ssa{i,Sp});
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

p=loglog(eps_n,err_t,'x-');hold on
p.Color=[1 0 0];

grid on; axis tight

p=loglog(eps_n,err_r,'o-');
p.Color = [ 0.929, 0.694, 0.1250 ];

p=loglog(eps_n,err_s,'s-');
p.Color=[0 0 0];



set(findall(fig,'-property','FontSize'),'FontSize',20)


ax = gca;
set(ax.Children, 'MarkerSize', 13)
set(ax.Children, 'LineWidth', 2)
ax.XLabel.String = '$\varepsilon$';
ax.YLabel.String = 'histogram distance';

% lh = legend('$\;$ $\tau$-leap','$\;$ r-leap','$\;$ s-leap');
% lh = legend('$\;$ $\tau$-leap','$\;$ adaptive $\tau$-leap','$\;$ r-leap','$\;$ s-leap','$\;$ adaptive s-leap');
lh = legend('$\;$ $\tau$-leap','$\;$ r-leap','$\;$ s-leap');
lh.FontSize = 20;
lh.Location='best';


%%


saveas(gcf,'error_mean', 'epsc' );







