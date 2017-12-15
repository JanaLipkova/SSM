%=================
%
%   Function to plot outputs from SSM
%   Decaying Dimerization
%===================


%=================
%
%   Function to plot outputs from SSM
%   Decaying Dimerization
%===================


function plotting_stiff
close all

bins = 50;
real = 5000;
%========%========================
%    Adaptive TAU case 1 e = 0.05
%=========%=======================

load  AdaptiveTau_stiff_dimerizaition_case2_e_0_05_histogram.txt ;
dataAT5 =  AdaptiveTau_stiff_dimerizaition_case2_e_0_05_histogram;

% compute mean and variance
meanAT5 = mean(dataAT5)
varAT5 = var(dataAT5)


 %optimal binning:
%   data = dataAT5(:,2) + rand(length(dataAT5(:,2)),1);
%   bins = optBINS(data',10,100)

[AT5a,AT5b] = hist(dataAT5(:,2),bins);

% normalization
w = AT5b(2) - AT5b(1);
AT5a = AT5a./(sum(AT5a)*w);

% ---------------------

%========%========================
%    Adaptive TAU case 1 e = 0.03
%=========%=======================
 
load  AdaptiveTau_stiff_dimerizaition_case2_e_0_03_histogram.txt ;
dataAT3 =  AdaptiveTau_stiff_dimerizaition_case2_e_0_03_histogram;
 
% compute mean and variance
meanAT3 = mean(dataAT3)
varAT3 = var(dataAT3)
 
 
 %optimal binning:
%   data = dataAT3(:,2) + rand(length(dataAT3(:,2)),1);
%   bins = optBINS(data',10,100)
 
[AT3a,AT3b] = hist(dataAT3(:,2),bins);
 
% normalization
w = AT3b(2) - AT3b(1);
AT3a = AT3a./(sum(AT3a)*w);
 
% ---------------------

%========%========================
%    Adaptive TAU case 1 e = 0.01
%=========%=======================
 
load  AdaptiveTau_stiff_dimerizaition_case2_e_0_01_histogram.txt ;
dataAT1 =  AdaptiveTau_stiff_dimerizaition_case2_e_0_01_histogram;
 
% compute mean and variance
meanAT1 = mean(dataAT1)
varAT1 = var(dataAT1)
 
 
 %optimal binning:
%   data = dataAT1(:,2) + rand(length(dataAT1(:,2)),1);
%   bins = optBINS(data',10,100)
 
[AT1a,AT1b] = hist(dataAT1(:,2),bins);
 
% normalization
w = AT1b(2) - AT1b(1);
AT1a = AT1a./(sum(AT1a)*w);
 
% ---------------------

%===============================================%===============================================

%========%========================
%   R-Leap case 5 e = 0.05
%=========%=======================
 
load  R_leapingJana_stiff_case2_e0_05_histogram_RLeaping.txt ;
dataR5 =  R_leapingJana_stiff_case2_e0_05_histogram_RLeaping;
 
% compute mean and variance
meanR5 = mean(dataR5)
varR5 = var(dataR5)
 
 
 %optimal binning:
%   data = dataR5(:,2) + rand(length(dataR5(:,2)),1);
%   bins = optBINS(data',10,100)
 
[R5a,R5b] = hist(dataR5(:,2),bins);
 
% normalization
w = R5b(2) - R5b(1);
R5a = R5a./(sum(R5a)*w);

%========%========================
%   R-Leap case 1 e = 0.03
%=========%=======================
 
load  R_leapingJana_stiff_case2_e0_03_histogram_RLeaping.txt ;
dataR3 =  R_leapingJana_stiff_case2_e0_03_histogram_RLeaping;
 
% compute mean and variance
meanR3 = mean(dataR3)
varR3 = var(dataR3)
 
 
 %optimal binning:
%   data = dataR3(:,2) + rand(length(dataR3(:,2)),1);
%   bins = optBINS(data',10,100)
 
[R3a,R3b] = hist(dataR3(:,2),bins);
 
% normalization
w = R3b(2) - R3b(1);
R3a = R3a./(sum(R3a)*w);


%========%========================
%   R-Leap case 1 e = 0.01
%=========%=======================
 
load  R_leapingJana_stiff_case2_e0_01_histogram_RLeaping.txt ;
dataR1 =  R_leapingJana_stiff_case2_e0_01_histogram_RLeaping;
 
% compute mean and variance
meanR1 = mean(dataR1)
varR1 = var(dataR1)
 
 
 %optimal binning:
%   data = dataR1(:,2) + rand(length(dataR1(:,2)),1);
%   bins = optBINS(data',10,100)
 
[R1a,R1b] = hist(dataR1(:,2),bins);
 
% normalization
w = R1b(2) - R1b(1);
R1a = R1a./(sum(R1a)*w);
 
% ---------------------

%===============================================%===============================================

%---------------------

%==========%==========
% TAU LEAPING
%==========%==========

% ----- e= 0.05 --------
load  TAU_leaping_stiff_case2_e0_05_histogram_tauLeaping.txt;
dataT5 =  TAU_leaping_stiff_case2_e0_05_histogram_tauLeaping;

meanT5 = mean(dataT5)
varT5 = var(dataT5)

[T5a,T5b] = hist(dataT5(:,2),bins);

% normalization
w = T5b(2) - T5b(1);
T5a = T5a./(sum(T5a)*w);


%----------------------

% ----- e= 0.03 --------
load  TAU_leaping_stiff_case2_e0_03_histogram_tauLeaping.txt;
dataT3 =  TAU_leaping_stiff_case2_e0_03_histogram_tauLeaping;
 
meanT3 = mean(dataT3)
varT3 = var(dataT3)

[T3a,T3b] = hist(dataT3(:,2),bins);

% normalization
w = T3b(2) - T3b(1);
T3a = T3a./(sum(T3a)*w);


%----------------------


% ----- e= 0.01 --------
load  TAU_leaping_stiff_case2_e0_01_histogram_tauLeaping.txt;
dataT1 =  TAU_leaping_stiff_case2_e0_01_histogram_tauLeaping;
 
meanT1 = mean(dataT1)
varT1 = var(dataT1)

[T1a,T1b] = hist(dataT1(:,2),bins);

% normalization
w = T1b(2) - T1b(1);
T1a = T1a./(sum(T1a)*w);


%----------------------


%==================
%   PLOTTING
%==================

% ---------------------
f=figure(1)
set(gca,'Fontsize',20);
hold on
plot(T5b,T5a,'bo')
plot(R5b,R5a,'rv')
 plot(AT5b,AT5a,'gs')
legend('Tau-Leap','R-Leap','AT-Leap')
xlabel('X2(10)')
title('Dec dim \epsilon = 0.05')
box on
grid on
% saveas(f,'DecDimDistEps005.fig')
% saveas(f,'DecDimDistEps005.pdf')

% %---------------------
% 
% % ---------------------
f=figure(2)
set(gca,'Fontsize',20);
hold on
plot(T3b,T3a,'bo')
plot(R3b,R3a,'rv')
 plot(AT3b,AT3a,'gs')
legend('Tau-Leap','R-Leap','AT-Leap')
xlabel('X2(10)')
title('Dec dim \epsilon = 0.03')
box on
grid on
% saveas(f,'DecDimDistEps003.fig')
% saveas(f,'DecDimDistEps003.pdf')
% %---------------------
% 
% %---------------------
f=figure(3)
set(gca,'Fontsize',20);
hold on
plot(T1b,T1a,'bo')
plot(R1b,R1a,'rv')
 plot(AT1b,AT1a,'gs')
legend('Tau-Leap','R-Leap','AT-Leap')
xlabel('X2(0)')
title('Dec dim \epsilon = 0.01')
box on
grid on
% saveas(f,'DecDimDistEps001.fig')
% saveas(f,'DecDimDistEps001.pdf')






