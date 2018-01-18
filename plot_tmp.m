clc;
clear

s=load('all-times.txt');

x=load('BSubtilis_Output.txt');

figure();clf
sp=2;

h=stairs(s(:,1),s(:,1+sp),'r'); hold on;

h.LineWidth = 2;

h=plot( x(:,1), x(:,1+sp), 'k.' );
h.MarkerSize = 10;