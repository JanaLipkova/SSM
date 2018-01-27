clear; clc


eps = {'0.01', '0.03', '0.05'};

load('SSA/hist.mat');
f_ssa=frq;
e_ssa=edges;

N = size(e_ssa,1);
M = size(e_ssa,2);
%% ========================================================================


 col = get(groot,'DefaultAxesColorOrder');

err_t = zeros(1,length(eps));
err_r = zeros(1,length(eps));
err_s = zeros(1,length(eps));
eps_n = zeros(1,length(eps));
for k = 1:length(eps)

    eps_n(k) = str2double(eps{k});
    
    load([ 'AdaptiveS/eps_' eps{k} '_hist.mat']);
    f_sleap=frq;
    e_sleap=edges;

    load([ 'AdaptiveTau/eps_' eps{k} '_hist.mat']);
    f_tleap=frq;
    e_tleap=edges;

    load([ 'RLeapingJana/eps_' eps{k} '_hist.mat']);
    f_rleap=frq;
    e_rleap=edges;
    
    
    for i = 2:N-1
        for Sp = 1:M
            h = e_ssa{i,Sp}(2) - e_ssa{i,Sp}(1);
            err_t(k) = err_t(k) + h*sum( abs(f_tleap{i,Sp}-f_ssa{i,Sp}) );
            err_r(k) = err_r(k) + h*sum( abs(f_rleap{i,Sp}-f_ssa{i,Sp}) );
            err_s(k) = err_s(k) + h*sum( abs(f_sleap{i,Sp}-f_ssa{i,Sp}) );
        end
    end
    
    
    
end

err_t = err_t/(N-2)/M;
err_r = err_r/(N-2)/M;
err_s = err_s/(N-2)/M;   

%%
fig=figure; clf

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex'); 

p=loglog(eps_n,err_t,'x-'); 
p.Color=[1 0 0];

hold on; grid on; axis tight; box on;

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

set(gca,'GridLineStyle', '--','LineWidth',1.5);
set(gca,'GridAlpha',0.4);

lh = legend('$\;$ $\tau$-leap','$\;$ r-leap', '$\;$ s-leap');
lh.FontSize = 20;
% lh.Location='best';
lh.Location='northwest';

%%

saveas(gcf,'error_nonstiff_dim', 'epsc' );

