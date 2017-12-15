%======================
%
%  PLOT ROOT FINDER JACOBIAN TEST
%
%======================

function RootFinderPlot

close all
clear all 

% load DimerizationAdaptiveTau_STIFF_0_05_Output.txt;
% data1 = DimerizationAdaptiveTau_STIFF_0_05_Output;

load DimerizationAdaptiveTAU_Jacobian_Output.txt;
data1 = DimerizationAdaptiveTAU_Jacobian_Output;


N = 1;
figure(1)
hold on
plot(data1(:,1),data1(:,2)*N,'b','Linewidth',1.5)
plot(data1(:,1),data1(:,3)*N,'b','Linewidth',1.5)
plot(data1(:,1),data1(:,4)*N,'b','Linewidth',1.5)
legend('S_{1}','S_{2}','S_{3}')
box on
grid on

% load DimerizationASCL_stiff_Output.txt
% data1 = DimerizationASCL_stiff_Output;
% hold on
% plot(data1(:,1),data1(:,2),'g*')
% plot(data1(:,1),data1(:,3),'g*')
% plot(data1(:,1),data1(:,4),'g*')
% legend('S_{1}','S_{2}','S_{3}')
% box on
% grid on


load DimerizationAdaptiveS_stiff2a_Output.txt 
data1 = DimerizationAdaptiveS_stiff2a_Output;
N = 1;
plot(data1(:,1),data1(:,2)*N,'k*')
plot(data1(:,1),data1(:,3)*N,'k*')
plot(data1(:,1),data1(:,4)*N,'k*')
legend('S_{1}','S_{2}','S_{3}')
box on
grid on

