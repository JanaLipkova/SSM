clear
system = 'Dimerization';
% method = {'AdaptiveS','AdaptiveTau','RLeapingJana','SLeaping','SSA'};
method = {'AdaptiveS','AdaptiveTau','RLeapingJana','SLeaping'};
% method = {'AdaptiveTau','RLeapingJana','SLeaping','SSA'};

eps = {0.05};

cfolder = pwd;

for k=1:length(method)
    
    if( ~strcmp(method{k},'SSA') )
        K = length(eps);
    else
        K = 1;
    end
    
    
    for l=1:K
    
        cd(cfolder);
        
        if( ~strcmp(method{k},'SSA') )
            insert = [ 'eps_' num2str(eps{l}) '_' ];
        else
            insert = [];
        end
            
        tmpf = [system '_' method{k} '_' insert 'Trajectories'];
        folder = fullfile( method{k}, tmpf );    
        cd(folder)
    
        N = numel(dir('*.txt'));    
        N = min(1e4,N);

        file = [system '_' method{k} '_traj_' num2str(1) '.txt'];
        x = load(file);
        M = length(x);
        
        n = size(x,1);
        m = size(x,2);
        
        d = zeros(N,n,m);
        
        for i=1:N
            file = [system '_' method{k} '_traj_' num2str(i) '.txt'];
            d(i,:,:) = load(file);
            if(mod(i,100)==0)
                disp(i)
            end
        end

        
        cd('../')
        file = [ insert 'trj' '.mat' ];
        save(file, 'd');
        
    end
end



cd(cfolder);