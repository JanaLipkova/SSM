clear
system = 'Dimerization';
% method = {'AdaptiveS','AdaptiveTau','RLeapingJana','SSA'};
method = {'AdaptiveS'};

eps = { 0.01 0.03};

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
    
        files = dir('*.txt');
        file_name = files(1).name;
        x = load(file_name);
        M = length(x);
        
        n = size(x,1);
        m = size(x,2);
        
        
        N = numel(files);    
        d = zeros(N,n,m);
        
        i=1;
        for file = files'
            
            file_name = file.name;

            d(i,:,:) = load(file_name);
            if(mod(i,100)==0)
                disp(i)
            end
            i=i+1;
        end

        
        cd('../')
        file = [ insert 'trj' '.mat' ];
        save(file, 'd');
        
    end
end



cd(cfolder);