clc;
clear


figure();clf

sp=6;

for i=1:10
    s=load(['lacy_lacz2_ssa' num2str(i) '.txt']);
    h=plot(s(:,1),s(:,1+sp),'r'); hold on;
    h.LineWidth = 2;
end



%%

for i=1:10
    x=load(['lacy_lacz2_rleap' num2str(i) '.txt']);
    h=plot( x(:,1), x(:,1+sp), 'k' );
    h.LineWidth = 2;
end