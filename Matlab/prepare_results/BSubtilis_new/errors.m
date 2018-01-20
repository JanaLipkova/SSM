clear
system = 'Dimerization';
% method = {'AdaptiveS','AdaptiveTau','RLeapingJana'};
method = { 'SLeaping', 'AdaptiveTau','RLeapingJana'};
% method = { 'SLeaping_SSA', 'TauLeaping_SSA'};
% method = {'AdaptiveS','SSA'};

eps = {0.01,0.03,0.05};

cfolder = pwd;


%% histogram for SSA to obtain edges

folder = 'SSA';
cd(folder)
file = 'trj.mat';
load(file);

N = size(d,2);
M = size(d,3)-1;
mn = cell(N,M);

Ns = 6000;

for i=1:N
    for j=1:M
        x = squeeze(d(1:Ns,i,j+1));
        mn{i,j} = mean(x);
    end
end
t = squeeze(d(1,:,1));

file = 'mean.mat';
save(file,'mn','t');
cd(cfolder)
%%

for k=1:length(method)    
    for l=1:length(eps)
    
        folder = method{k};
        cd(folder);
    
        insert = [ 'eps_' num2str(eps{l}) '_' ];
        
        file = [ insert 'trj' '.mat' ];
        load(file);
 
        N = size(d,2);
        M = size(d,3)-1;
        mn= cell(N,M);

        
        for i=1:N
            for j=1:M
                x = squeeze(d(1:Ns,i,j+1));
                mn{i,j} = mean(x);
            end
        end
        
        t = squeeze(d(1,:,1));
        
        file = [ insert 'mean' '.mat' ];
        save(file,'mn','t');
        
        cd(cfolder)
    end
end
