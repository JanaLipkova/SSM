clear
system = 'Dimerization';
method = {'AdaptiveS','AdaptiveTau','RLeapingJana'};
% method = {'TauLeaping'};
% method = {'AdaptiveS'};

eps = { 0.01, 0.03,0.05};

cfolder = pwd;


%% histogram for SSA to obtain edges

folder = 'SSA';
cd(folder)
file = 'trj.mat';
load(file);

N = size(d,2);
M = size(d,3)-1;
frq   = cell(N,M);
edges = cell(N,M);

Ns = 5400;

for i=1:N
    for j=1:M
        x = squeeze(d(1:Ns,i,j+1));
        [frq{i,j} edges{i,j}] = histcounts(x,'Normalization','pdf');
    end
end
t = squeeze(d(1,:,1));

file = 'hist.mat';
save(file,'frq','edges','t');
cd(cfolder)
%%
e_ssa = edges;

for k=1:length(method)    
    
    disp(method{k})
    
    for l=1:length(eps)
    
        folder = method{k};
        cd(folder);
    
        insert = [ 'eps_' num2str(eps{l}) '_' ];
        
        file = [ insert 'trj' '.mat' ];
        load(file);
 
        N = size(d,2);
        M = size(d,3)-1;
        frq   = cell(N,M);
        edges = cell(N,M);
        
        for i=1:N
            for j=1:M
                x = squeeze( d(1:Ns,i,j+1) );
                [frq{i,j} edges{i,j}] = histcounts(x, e_ssa{i,j}, 'Normalization','pdf');
            end
        end
        
        t = squeeze(d(1,:,1));
        
        file = [ insert 'hist' '.mat' ];
        save(file,'frq','edges','t');
        
        cd(cfolder)
    end
end
