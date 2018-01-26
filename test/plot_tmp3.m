

s=load('all-times.txt');



fig=figure(1);clf

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex'); 



h=stairs(s(:,1),s(:,2:end) ); hold on;

grid on; axis tight

set(findall(fig,'-property','FontSize'),'FontSize',20)


ax = gca;
set(ax.Children, 'MarkerSize', 13)
set(ax.Children, 'LineWidth', 2)
ax.XLabel.String = '$t$';
ax.YLabel.String = 'number of molecules';



lh = legend( '$S_1$', '$S_2$', '$S_3$' );

lh.FontSize = 20;
lh.Location='best';



%%


saveas(gcf,'trajectory_bsubtilis', 'epsc' );
