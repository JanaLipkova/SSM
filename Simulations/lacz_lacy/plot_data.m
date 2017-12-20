clear

x1=load('out_ssa.txt');
x2=load('out_tau.txt');
x3=load('out_rleap.txt');
x4=load('out_sleap.txt');


%%

k=20;

plot(x1(1:end-1,1),x1(1:end-1,k)); hold on
plot(x2(1:end-1,1),x2(1:end-1,k))
plot(x3(1:end-1,1),x3(1:end-1,k))
plot(x4(1:end-1,1),x4(1:end-1,k))

