%===============
%
%  SSM vs Matlab
%
%===============

function plot_tests

close all
clear all

% NON STIFF

% load DimerizationSSA_Output.txt;
% dataS = DimerizationSSA_Output;
 N = 1;%0000;

% f=figure(2)
% hold on
% set(gca,'Fontsize',20);
% plot(dataS(:,1),dataS(:,2)*N,'k','Linewidth',1.5)
% plot(dataS(:,1),dataS(:,3)*N,'k','Linewidth',1.5)
% plot(dataS(:,1),dataS(:,4)*N,'k','Linewidth',1.5)
% legend('SSA')
% % N=100;
% load Dimerization_SLeap_testA_Output.txt;
% data1 = Dimerization_SLeap_testA_Output;
% hold on
% plot(data1(:,1),data1(:,2)*N,'m*')
% plot(data1(:,1),data1(:,3)*N,'m*')
% plot(data1(:,1),data1(:,4)*N,'m*')
% box on
% grid on


%STIFF
load DimerizationSSA_Stiff_Output.txt
dataS = DimerizationSSA_Stiff_Output;

f=figure(2)
% hold on
% set(gca,'Fontsize',20);
% plot(dataS(:,1),dataS(:,2)*N,'k','Linewidth',1.5)
% plot(dataS(:,1),dataS(:,3)*N,'k','Linewidth',1.5)
% plot(dataS(:,1),dataS(:,4)*N,'k','Linewidth',1.5)
% legend('SSA')


load DimerizationAdaptiveTAU_Output.txt;
data1 = DimerizationAdaptiveTAU_Output
hold on
plot(data1(:,1),data1(:,2)*N,'r*')
plot(data1(:,1),data1(:,3)*N,'b*')
plot(data1(:,1),data1(:,4)*N,'g*')
box on
grid on


