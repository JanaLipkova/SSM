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
err_t_a = zeros(1,length(eps));
err_r = zeros(1,length(eps));
err_s = zeros(1,length(eps));
err_s_a = zeros(1,length(eps));
eps_n = zeros(1,length(eps));

 for k = 1:length(eps)

    eps_n(k) = str2double(eps{k});
    
    load([ 'SLeaping/eps_' eps{k} '_hist.mat']);
    f_sleap=frq;
    e_sleap=edges;

    load([ 'AdaptiveS/eps_' eps{k} '_hist.mat']);
    f_sleap_a=frq;
    e_sleap_a=edges;
    
    load([ 'TauLeaping/eps_' eps{k} '_hist.mat']);
    f_tleap=frq;
    e_tleap=edges;

    load([ 'AdaptiveTau/eps_' eps{k} '_hist.mat']);
    f_tleap_a = frq;
    e_tleap_a = edges;
    
    load([ 'RLeapingJana/eps_' eps{k} '_hist.mat']);
    f_rleap=frq;
    e_rleap=edges;
    
    
%     for i = 2:N-1
    for i = N-1
        for Sp = 1:M
            h = e_ssa{i,Sp}(2) - e_ssa{i,Sp}(1);
            err_t(k) = err_t(k) + h*sum( abs(f_tleap{i,Sp}-f_ssa{i,Sp}) );
            err_t_a(k) = err_t_a(k) + h*sum( abs(f_tleap_a{i,Sp}-f_ssa{i,Sp}) );
            err_r(k) = err_r(k) + h*sum( abs(f_rleap{i,Sp}-f_ssa{i,Sp}) );
            err_s(k) = err_s(k) + h*sum( abs(f_sleap{i,Sp}-f_ssa{i,Sp}) );
            err_s_a(k) = err_s_a(k) + h*sum( abs(f_sleap_a{i,Sp}-f_ssa{i,Sp}) );
        end
    end
    
    
    
end

err_t = err_t/(N-2)/M;
err_t_a = err_t_a/(N-2)/M;
err_r = err_r/(N-2)/M;
err_s = err_s/(N-2)/M;
err_s_a = err_s_a/(N-2)/M;

%%
fig=figure(1); clf

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex'); 

p=loglog(eps_n,err_t,'x-');hold on
p.Color=[1 0 0];
p=loglog(eps_n,err_t_a,'x--');
p.Color=[1 0 0];

grid on; axis tight

loglog(eps_n,err_r,'o-')

p=loglog(eps_n,err_s,'s-');
p.Color=[0 0 0];
p=loglog(eps_n,err_s_a,'s--');
p.Color=[0 0 0];



set(findall(fig,'-property','FontSize'),'FontSize',20)


ax = gca;
set(ax.Children, 'MarkerSize', 13)
set(ax.Children, 'LineWidth', 2)
ax.XLabel.String = '$\varepsilon$';
ax.YLabel.String = 'histogram distance';

% lh = legend('$\;$ $\tau$-leap','$\;$ r-leap','$\;$ s-leap');
lh = legend('$\;$ $\tau$-leap','$\;$ adaptive $\tau$-leap','$\;$ r-leap','$\;$ s-leap','$\;$ adaptive s-leap');
% lh = legend('$\;$ adaptive $\tau$-leap','$\;$ r-leap','$\;$ adaptive s-leap');
lh.FontSize = 20;
lh.Location='best';


%%


saveas(gcf,'error_stiff_dim', 'epsc' );







