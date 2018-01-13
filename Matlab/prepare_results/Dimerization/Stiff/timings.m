
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




fig = figure(1); clf

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex'); 


p=semilogy(s1(:,1),1./s1(:,2),'x-'); hold on
p.Color=[1 0 0];
p=loglog(s1a(:,1),1./s1a(:,2),'x--'); hold on
p.Color=[1 0 0];

p=loglog(s2(:,1),1./s2(:,2),'o-'); 
p.Color = [ 0.929, 0.694, 0.1250 ];

p=semilogy(s3(:,1),1./s3(:,2),'s-'); 
p.Color=[0 0 0];
p=loglog(s3a(:,1),1./s3a(:,2),'s--'); 
p.Color=[0 0 0];


grid on; axis tight

set(findall(gcf,'-property','FontSize'),'FontSize',20)

ax = gca;
set(ax.Children, 'MarkerSize', 13)
set(ax.Children, 'LineWidth', 2)
ax.XLabel.String = '$\varepsilon$';
ax.YLabel.String = 'speed-up';

lh = legend('$\;$ $\tau$-leap','$\;$ adaptive $\tau$-leap','$\;$ r-leap','$\;$ s-leap','$\;$adaptive s-leap');
% lh = legend('$\;$ adaptive $\tau$-leap','$\;$ r-leap','$\;$ adaptive s-leap');
lh.FontSize = 20;
lh.Location='best';

%%


saveas(gcf,'time_stiff_dim', 'epsc' );


