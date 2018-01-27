
ssa=load('Timings/LacZLacY-big-T100-SSA_times.txt');
s = ssa(1,2);
LB=4;
UB=6;

% s1=load('Timings/LacZLacY-big-T100-TauLeaping_times.txt');
% s1 = s1(LB:UB,:);
% s1(:,2)=s1(:,2)/s;

s1a=load('Timings/LacZLacY-big-T100-AdaptiveTau_times.txt');
s1a = s1a(LB:UB,:);
s1a(:,2)=s1a(:,2)/s;

s2=load('Timings/LacZLacY-big-T100-RLeaping_times.txt');
s2 = s2(LB:UB,:);
s2(:,2)=s2(:,2)/s;

s3=load('Timings/LacZLacY-big-T100-SLeaping_v3_times.txt');
s3 = s3(LB:UB,:);
s3(:,2)=s3(:,2)/s;

% s4=load('Timings/LacZLacY-big-T100-SLeaping_v4_times.txt');
% s4 = s4(LB:UB,:);
% s4(:,2)=s4(:,2)/s;

s5=load('Timings/LacZLacY-big-T100-SLeaping_v5_times.txt');
s5 = s5(LB:UB,:);
s5(:,2)=s5(:,2)/s;

% s6=load('Timings/LacZLacY-big-T100-SLeaping_v6_times.txt');
% s6 = s6(LB:UB,:);
% s6(:,2)=s6(:,2)/s;


fig = figure; clf

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex'); 


% p=semilogy(s1(:,1),1./s1(:,2),'x-'); hold on
% p.Color=[1 0 0];
p=loglog(s1a(:,1),1./s1a(:,2),'x--'); hold on
p.Color=[1 0 0];

p=loglog(s2(:,1),1./s2(:,2),'o-'); 
p.Color = [ 0.929, 0.694, 0.1250 ];

p=semilogy(s3(:,1),1./s3(:,2),'s-'); 
p.Color=[0 0 0];


% p=semilogy(s4(:,1),1./s4(:,2),'s-'); 
% p.Color=[0 0 0];

% p=semilogy(s5(:,1),1./s5(:,2),'p-'); 
% p.Color=[0 0 0];

% p=semilogy(s6(:,1),1./s6(:,2),'s--'); 
% p.Color=[0 0 0];
axis([0.01, 0.05, 1,14])

grid on; axis tight

set(findall(gcf,'-property','FontSize'),'FontSize',20)

ax = gca;
set(ax.Children, 'MarkerSize', 13)
set(ax.Children, 'LineWidth', 2)
ax.XLabel.String = '$\varepsilon$';
ax.YLabel.String = 'speed-up';

% lh = legend('$\;$ $\tau$-leap','$\;$ adaptive $\tau$-leap','$\;$ r-leap','$\;$ s-leap','$\;$adaptive s-leap');
lh = legend('$\;$ $\tau$-leap','$\;$ r-leap','$\;$ s-leap','$\;$ s-leap conservative');
lh.FontSize = 20;
lh.Location='northwest';

%%
set(gca,'GridLineStyle', '--','LineWidth',1.8);
set(gca,'GridAlpha',0.4);

saveas(gcf,'time_lac_big_dim', 'epsc' );


