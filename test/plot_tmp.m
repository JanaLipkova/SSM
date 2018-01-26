
s=load('all-times.txt');

x=load('BSubtilis_Output.txt');

figure();clf
hold on
sp=1;

h=stairs(s(:,1),s(:,1+sp),'r'); hold on;

h.LineWidth = 2;

h=plot( x(:,1), x(:,1+sp), 'k.' );
h.MarkerSize = 10;

% compare with ssa
% s=load('all-times.txt');
% x=load('BSubtilis_Output.txt');
% 
% h=stairs(s(:,1),s(:,1+sp),'b'); hold on;
% 
% h.LineWidth = 2;
% 
% h=plot( x(:,1), x(:,1+sp), 'go' );
% h.MarkerSize = 10;