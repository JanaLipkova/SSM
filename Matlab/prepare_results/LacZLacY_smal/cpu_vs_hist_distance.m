% clear; clc

eps = { '0.5',  '0.1', '0.05', '0.03', '0.01'};

load('SSA/hist.mat');
f_ssa=frq;
e_ssa=edges;

N = size(e_ssa,1);
M = size(e_ssa,2);

SpList = [6 13 14];

%% ========================================================================


col = get(groot,'DefaultAxesColorOrder');

err_t = zeros(1,length(eps));
err_t_a = zeros(1,length(eps));
err_r = zeros(1,length(eps));
err_r_c = zeros(1,length(eps));
err_s_v3 = zeros(1,length(eps));
err_s_v4 = zeros(1,length(eps));
err_s_v5 = zeros(1,length(eps));
err_s_v6 = zeros(1,length(eps));
err_s_c  = zeros(1,length(eps));
eps_n = zeros(1,length(eps));

 for k = 1:length(eps)

    eps_n(k) = str2double(eps{k});
    
    
    load([ 'TauLeap/eps_' eps{k} '_hist.mat']);
    f_tleap=frq;
    e_tleap=edges;
    
    load([ 'AdaptiveTau/eps_' eps{k} '_hist.mat']);
    f_tleap_a = frq;
    e_tleap_a = edges;
    
    
    
    load([ 'RLeaping/eps_' eps{k} '_hist.mat']);
    f_rleap=frq;
    e_rleap=edges;
    
    load([ 'RLeaping_c/eps_' eps{k} '_hist.mat']);
    f_rleap_c=frq;
    e_rleap_c=edges;
    
    
    
    load([ 'SLeaping_v3/eps_' eps{k} '_hist.mat']);
    f_sleap_v3=frq;
    e_sleap_v3=edges;
    
    
    load([ 'SLeaping_v4/eps_' eps{k} '_hist.mat']);
    f_sleap_v4=frq;
    e_sleap_v4=edges;

    load([ 'SLeaping_v5/eps_' eps{k} '_hist.mat']);
    f_sleap_v5=frq;
    e_sleap_v5=edges;
    
    load([ 'SLeaping_v6/eps_' eps{k} '_hist.mat']);
    f_sleap_v6=frq;
    e_sleap_v6=edges;
    
    load([ 'SLeaping_c/eps_' eps{k} '_hist.mat']);
    f_sleap_c=frq;
    e_sleap_c=edges;
    
    
    list = 1:N;
    n = length(list);
    m = length(SpList);
    
    for i = list
        for Sp = SpList
            h = e_ssa{i,Sp}(2) - e_ssa{i,Sp}(1);
            err_t(k)     = err_t(k)      +  h*sum( abs( f_tleap{i,Sp}    - f_ssa{i,Sp}) );
            err_t_a(k)   = err_t_a(k)    +  h*sum( abs( f_tleap_a{i,Sp}  - f_ssa{i,Sp}) );
            err_r(k)     = err_r(k)      +  h*sum( abs( f_rleap{i,Sp}    - f_ssa{i,Sp}) );
            err_r_c(k)   = err_r_c(k)    +  h*sum( abs( f_rleap_c{i,Sp}  - f_ssa{i,Sp}) );
            err_s_v3(k)  = err_s_v3(k)   +  h*sum( abs( f_sleap_v3{i,Sp} - f_ssa{i,Sp}) );
            err_s_v4(k)  = err_s_v4(k)   +  h*sum( abs( f_sleap_v4{i,Sp} - f_ssa{i,Sp}) );
            err_s_v5(k)  = err_s_v5(k)   +  h*sum( abs( f_sleap_v5{i,Sp} - f_ssa{i,Sp}) );
            err_s_v6(k)  = err_s_v6(k)   +  h*sum( abs( f_sleap_v6{i,Sp} - f_ssa{i,Sp}) );
            err_s_c(k)   = err_s_c(k)    +  h*sum( abs( f_sleap_c{i,Sp}  - f_ssa{i,Sp}) );
        end
    end
    
    
    
 end

M = length(SpList);
err_t    = err_t/n/m;
err_t_a  = err_t_a/n/m;
err_r    = err_r/n/m;
err_r_c  = err_r_c/n/m;
err_s_v3 = err_s_v3/n/m;
err_s_v4 = err_s_v4/n/m;
err_s_v5 = err_s_v5/n/m;
err_s_v6 = err_s_v6/n/m;
err_s_c  = err_s_c/n/m;

%% timing
ssa=load('Timings/LacZLacY-small-T2100-SSA_times.txt');
s = ssa(1,2);
LB=3;
UB=6;

s1=load('Timings/LacZLacY-small-T2100-TauLeaping_times.txt');
s1 = s1(LB:UB,:);
s1(:,2)=s1(:,2)/s;

s1a=load('Timings/LacZLacY-big-T2100-TauLeaping-control-noSSA.txt'); % adaptive SSA
s1a = s1a(LB:UB,:);
s1a(:,2)=s1a(:,2)/s;

s2=load('Timings/LacZLacY-small-T2100-RLeaping_times.txt');
s2 = s2(LB:UB,:);
s2(:,2)=s2(:,2)/s;

s3=load('Timings/LacZLacY-small-T2100-SLeaping_v3_times.txt');
s3 = s3(LB:UB,:);
s3(:,2)=s3(:,2)/s;

s4=load('Timings/LacZLacY-big-T2100-RLeaping-control-p-10M.txt');
s4 = s4(LB:UB,:);
s4(:,2)=s4(:,2)/s;

s5=load('Timings/LacZLacY-big-T2100-SLeaping_v6-control-p-10M.txt');
s5 = s5(LB:UB,:);
s5(:,2)=s5(:,2)/s;


%%

fig=figure; clf

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex'); 
hold on; grid on; axis tight; box on;

LB=2;
UB=5;

% p=loglog(err_t(LB:UB), 1./s1(:,2),'x-'); 
% p.Color=[1 0 0];

% p=loglog(err_t_a(LB:UB), 1./s1a(end:-1:1,2),'x--'); 
% p.Color=[1 0 0];
% 
% 
% p=loglog(err_r(LB:UB),1./s2(:,2),'o-');
% p.Color = [ 0.929, 0.694, 0.1250 ];

p=loglog(err_r_c(LB:UB),1./s4(:,2),'o--');
p.Color = [ 0.929, 0.694, 0.1250 ];

% p=loglog(err_s_v3(LB:UB),1./s3(end:-1:1,2),'s-');
% p.Color=[0 0 0];
% 
% p=loglog(err_s_c(LB:UB),1./s5(end:-1:1,2),'s--');
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

lh = legend('$\;$ r-leap $\theta=0.1$');
% lh = legend('$\;$ $\tau$-leap no control','$\;$ r-leap','$\;$ s-leap');
lh.FontSize = 20;
% lh.Location='best';
lh.Location='northwest';

set(gca,'GridLineStyle', '--','LineWidth',1.8);
set(gca,'GridAlpha',0.4);

%%

saveas(gcf,'error_vs_time_small_lac', 'epsc' );

