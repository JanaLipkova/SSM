%=============================================
%
%    TESTING PLOTTING SCRIPT
%     
%       1. load data
%       2. plot histogram ( may use optimal binning)
%       3. plot average realisations
%
%=============================================


function plotting_test

  close all;

  
% load  DimerizationAdaptiveSNonStiff_Output.txt;
% data1 =  DimerizationAdaptiveSNonStiff_Output;
% 
% R = 10000;
% 
% figure(1)
% hold on
% plot(data1(:,1),data1(:,2),'r')
% plot(data1(:,1),data1(:,3),'b')
% plot(data1(:,1),data1(:,4),'g')
% legend('S_{1}','S_{2}','S_{3}')
% box on
% grid on
% 
% 
% load  DimerizationAdaptiveTauNonStiff_Output.txt;
% data1 =  DimerizationAdaptiveTauNonStiff_Output;
% 
% % figure(2)
% hold on
% plot(data1(:,1),data1(:,2),'r*')
% plot(data1(:,1),data1(:,3),'b*')
% plot(data1(:,1),data1(:,4),'g*')
% legend('S_{1}','S_{2}','S_{3}')
% box on
% grid on
% 
% load  DimerizationARLeapNonStiif_Output.txt;
% data1 =  DimerizationARLeapNonStiif_Output;
% 
% % figure(3)
% hold on
% plot(data1(:,1),data1(:,2)*R,'rs')
% plot(data1(:,1),data1(:,3)*R,'bs')
% plot(data1(:,1),data1(:,4)*R,'gs')
% legend('S_{1}','S_{2}','S_{3}')
% box on
% grid on
% 
% load  DimerizationSLeapNonStiff_Output.txt;
% data1 =  DimerizationSLeapNonStiff_Output;
% 
% % figure(4)
% hold on
% plot(data1(:,1),data1(:,2)*R,'ro')
% plot(data1(:,1),data1(:,3)*R,'bo')
% plot(data1(:,1),data1(:,4)*R,'go')
% legend('S_{1}','S_{2}','S_{3}')
% box on
% grid on
% 
% 
% 
% load  DimerizationR_leapingJana_Output.txt;
% data1 =  DimerizationR_leapingJana_Output;
% 
% % figure(4)
% hold on
% plot(data1(:,1),data1(:,2),'k')
% plot(data1(:,1),data1(:,3),'k')
% plot(data1(:,1),data1(:,4),'k')
% legend('S_{1}','S_{2}','S_{3}')
% box on
% grid on


% HISTOGRAMS
bins=30
load DimerizationAdaptiveSNonStiff_histogram.txt;
data1 = DimerizationAdaptiveSNonStiff_histogram;

[a1,b1] = hist(data1(:,2),bins);
 
figure(10)
set(gca,'Fontsize',20);
hold on
plot(b1,a1,'b')
xlabel('X2(10)')
box on
grid on

load DimerizationAdaptiveTauNonStiff_histogram.txt;
data1 = DimerizationAdaptiveTauNonStiff_histogram;

[a1,b1] = hist(data1(:,2),bins);
 
% figure(20)
set(gca,'Fontsize',20);
hold on
plot(b1,a1,'g')
xlabel('X2(10)')
box on
grid on

load DimerizationARLeapNonStiif_histogram.txt;
data1 = DimerizationARLeapNonStiif_histogram;

[a1,b1] = hist(data1(:,2),bins);
 
% figure(20)
set(gca,'Fontsize',20);
hold on
plot(b1,a1,'r')
xlabel('X2(10)')
box on
grid on

load DimerizationSLeapNonStiff_histogram_SLeaping.txt;
data1 = DimerizationSLeapNonStiff_histogram_SLeaping;

[a1,b1] = hist(data1(:,2),bins);
 
% figure(20)
set(gca,'Fontsize',20);
hold on
plot(b1,a1,'k')
xlabel('X2(10)')
box on
grid on

