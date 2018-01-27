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
err_s_b = zeros(1,length(eps));
err_s_c = zeros(1,length(eps));
eps_n = zeros(1,length(eps));

 for k = 1:length(eps)

    eps_n(k) = str2double(eps{k});
    
    load([ 'SLeaping/eps_' eps{k} '_hist.mat']);
    f_sleap=frq;
    e_sleap=edges;

    load([ 'AdaptiveS/eps_' eps{k} '_hist.mat']);
    f_sleap_a=frq;
    e_sleap_a=edges;
    
    load([ 'AdaptiveS_a0_x_star/eps_' eps{k} '_hist.mat']);
    f_sleap_b=frq;
    e_sleap_b=edges;
    
    load([ 'AdaptiveS_a0_x_t/eps_' eps{k} '_hist.mat']);
    f_sleap_c=frq;
    e_sleap_c=edges;
    
    load([ 'TauLeaping/eps_' eps{k} '_hist.mat']);
    f_tleap=frq;
    e_tleap=edges;

    load([ 'AdaptiveTau/eps_' eps{k} '_hist.mat']);
    f_tleap_a = frq;
    e_tleap_a = edges;
    
    load([ 'RLeapingJana/eps_' eps{k} '_hist.mat']);
    f_rleap=frq;
    e_rleap=edges;
    
    
    for i = 2:N-1
%     for i = N-1
        for Sp = 1:M
            h = e_ssa{i,Sp}(2) - e_ssa{i,Sp}(1);
            err_t(k) = err_t(k) + h*sum( abs(f_tleap{i,Sp}-f_ssa{i,Sp}) );
            err_t_a(k) = err_t_a(k) + h*sum( abs(f_tleap_a{i,Sp}-f_ssa{i,Sp}) );
            err_r(k) = err_r(k) + h*sum( abs(f_rleap{i,Sp}-f_ssa{i,Sp}) );
            err_s(k) = err_s(k) + h*sum( abs(f_sleap{i,Sp}-f_ssa{i,Sp}) );
            err_s_a(k) = err_s_a(k) + h*sum( abs(f_sleap_a{i,Sp}-f_ssa{i,Sp}) );
            err_s_b(k) = err_s_b(k) + h*sum( abs(f_sleap_b{i,Sp}-f_ssa{i,Sp}) );
            err_s_c(k) = err_s_c(k) + h*sum( abs(f_sleap_c{i,Sp}-f_ssa{i,Sp}) );
        end
    end
    
    
    
end

err_t = err_t/(N-2)/M;
err_t_a = err_t_a/(N-2)/M;
err_r = err_r/(N-2)/M;
err_s = err_s/(N-2)/M;
err_s_a = err_s_a/(N-2)/M;
err_s_b = err_s_b/(N-2)/M;
err_s_c = err_s_c/(N-2)/M;


%% data for timing
s=load('Timings/Dimerization-Stiff-SSA_times.txt');

s1=load('Timings/Dimerization-Stiff-TauLeaping_times.txt');
s1(:,2)=s1(:,2)/s;

s1a=load('Timings/Dimerization-Stiff-AdaptiveTau_times.txt');
s1a(:,2)=s1a(:,2)/s;


s2=load('Timings/Dimerization-Stiff-RLeapingJana_times.txt');
s2(:,2)=s2(:,2)/s;

s3=load('Timings/Dimerization-Stiff-SLeaping_times.txt');
s3(:,2)=s3(:,2)/s;

s3a=load('Timings/Dimerization-Stiff-AdaptiveS_times.txt');
s3a(:,2)=s3a(:,2)/s;


%% plotting
fig=figure; clf

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex'); 

subplot(2,1,1)
p=loglog(err_t, 1./s1(end:-1:1,2),'x-'); 
p.Color=[1 0 0];

hold on; grid on; axis tight; box on;

% p=loglog(err_t_a, 1./s1a(end:-1:1,2),'x--'); 
% p.Color=[1 0 0];

p=loglog(err_r,1./s2(end:-1:1,2),'o-');
p.Color = [ 0.929, 0.694, 0.1250 ];

p=loglog(err_s,1./s3(end:-1:1,2),'s-');
p.Color=[0 0 0];

% p=loglog(err_s_a,1./s3a(end:-1:1,2),'s--');
% p.Color=[0 0 0];

set(findall(fig,'-property','FontSize'),'FontSize',20)

% axis([0.06,0.6,10,1000])

ax = gca;
set(ax.Children, 'MarkerSize', 13)
set(ax.Children, 'LineWidth', 2)
ax.YLabel.String = 'speed up';
ax.XLabel.String = 'histogram distance';

set(gca,'GridLineStyle', '--','LineWidth',1.5);
set(gca,'GridAlpha',0.5);

lh = legend('$\;$ $\tau$-leap','$\;$ r-leap','$\;$ s-leap');
% lh = legend('$\;$ $\tau$-leap','$\;$ adaptive $\tau$-leap','$\;$ r-leap','$\;$ s-leap','$\;$ adaptive s-leap');
lh.FontSize = 20;
% lh.Location='best';
lh.Location='southeast';

set(gca,'GridLineStyle', '--','LineWidth',1.8);
set(gca,'GridAlpha',0.4);



subplot(2,1,2)

p=loglog(err_t_a, 1./s1a(end:-1:1,2),'x--'); 
p.Color=[1 0 0];
hold on
p=loglog(err_s_a,1./s3a(end:-1:1,2),'s--');
p.Color=[0 0 0];

set(findall(fig,'-property','FontSize'),'FontSize',20)

% axis([0.06,0.6,10,1000])

ax = gca;
set(ax.Children, 'MarkerSize', 13)
set(ax.Children, 'LineWidth', 2)
ax.YLabel.String = 'speed up';
ax.XLabel.String = 'histogram distance';

set(gca,'GridLineStyle', '--','LineWidth',1.5);
set(gca,'GridAlpha',0.5);

lh = legend('$\;$ adaptive $\tau$-leap','$\;$ adaptive s-leap');
lh.FontSize = 20;
% lh.Location='best';
lh.Location='southeast';

set(gca,'GridLineStyle', '--','LineWidth',1.8);
set(gca,'GridAlpha',0.4);
box on; grid on;
%%

saveas(gcf,'error_vs_time_stiff_dim', 'epsc' );
