
s=load('BSubtilis_SSA_Output.txt');
cmap    = get(0, 'defaultaxescolororder');



fig = figure(1); clf
hold on
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex'); 
p=plot( s(:,1), s(:,2));
p.Color=cmap(5,:);

p=plot( s(:,1), s(:,3));
p.Color=cmap(2,:);

p=plot( s(:,1), s(:,4));
p.Color=cmap(6,:);

legend('S_1','S_2','S_3')
set(gca,'LineWidth',2);
grid on; axis tight
set(findall(gcf,'-property','FontSize'),'FontSize',20)

ax = gca;
set(ax.Children, 'MarkerSize', 13)
set(ax.Children, 'LineWidth', 3)
ax.XLabel.String = '$t$';
ax.YLabel.String = 'number of molecules';
lh = legend('$\; S_1$','$\; S_2$','$\; S_3$');
lh.FontSize = 20;
lh.Location='best';
box on
%%


saveas(gcf,'BSubtilis-one-sample', 'epsc' );


