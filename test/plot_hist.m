s1=load('hist_ssa.txt');
s2=load('hist_r.txt');
s3=load('hist_s.txt');

sp = 13;

[f1, e1] = histcounts(s1(:,sp),   'Normalization','pdf');
[f2, e2] = histcounts(s2(:,sp),e1,'Normalization','pdf');
[f3, e3] = histcounts(s3(:,sp),e1,'Normalization','pdf');

m = e1(1:end-1)+(e1(2:end)-e1(1:end-1))/2;

plot(m,f1,'o-','LineWidth',3,'MarkerSize',10); hold on
plot(m,f2,'s-','LineWidth',3,'MarkerSize',10); hold on
plot(m,f3,'x-','LineWidth',3,'MarkerSize',10); hold on