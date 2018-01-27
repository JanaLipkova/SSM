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

%% load timings

s=load('Timings/Dimerization-NonStiff-SSA_times.txt');

% s1=load('Timings/Dimerization-NonStiff-TauLeaping_times.txt');
% s1(:,2)=s1(:,2)/s;

s1a=load('Timings/Dimerization-NonStiff-AdaptiveTau_times.txt');
s1a(:,2)=s1a(:,2)/s;


s2=load('Timings/Dimerization-NonStiff-RLeapingJana_times.txt');
s2(:,2)=s2(:,2)/s;

s3=load('Timings/Dimerization-NonStiff-SLeaping_times.txt');
s3(:,2)=s3(:,2)/s;

% s3a=load('Timings/Dimerization-NonStiff-AdaptiveS_times.txt');
% s3a(:,2)=s3a(:,2)/s;


%%
fig=figure; clf

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex'); 

p=loglog(err_t, 1./s1a(end:-1:1,2),'x-'); 
p.Color=[1 0 0];

hold on; grid on; axis tight; box on;

p=loglog(err_r,1./s2(end:-1:1,2),'o-');
p.Color = [ 0.929, 0.694, 0.1250 ];

p=loglog(err_s,1./s3(end:-1:1,2),'s-');
p.Color=[0 0 0];

set(findall(fig,'-property','FontSize'),'FontSize',20)

axis([0.06,0.6,10,1000])

ax = gca;
set(ax.Children, 'MarkerSize', 13)
set(ax.Children, 'LineWidth', 2)
ax.YLabel.String = 'speed up';
ax.XLabel.String = 'histogram distance';

set(gca,'GridLineStyle', '--','LineWidth',1.5);
set(gca,'GridAlpha',0.5);

lh = legend('$\;$ $\tau$-leap','$\;$ r-leap', '$\;$ s-leap');
lh.FontSize = 20;
% lh.Location='best';
lh.Location='northwest';

set(gca,'GridLineStyle', '--','LineWidth',1.8);
set(gca,'GridAlpha',0.4);

%%

saveas(gcf,'error_vs_time_nonstiff_dim', 'epsc' );

