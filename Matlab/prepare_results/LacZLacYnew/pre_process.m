clear
system = 'lacy_lacz2';
% method = {'SSA','TauLeap','AdaptiveTau','RLeaping','SLeaping_v3','SLeaping_v4','SLeaping_v5','AdaptiveS'};
% method = {'TauLeap','AdaptiveTau','RLeaping','SLeaping_v3','SLeaping_v4','SLeaping_v5','AdaptiveS'};
method = { 'SLeaping_v3_1' };


eps = {0.01 0.03 0.05};
% eps = {0.01};

cfolder = pwd;

for k=1:length(method)
    
    if( ~strcmp(method{k},'SSA') )
        K = length(eps);
    else
        K = 1;
    end
    
    
    for l=1:K
    
        cd(cfolder);
        
        cd(method{k})
        
        fprintf('%s --  %f \n',method{k},eps{l});

        
        if( ~strcmp(method{k},'SSA') )
            insert = [ 'eps_' num2str(eps{l}) '_' ];
        else
            insert = [];
        end
         
%         tmpf = [system '_' method{k} '_' insert 'Trajectories'];
        tmpf = [system '_' insert 'Trajectories'];
        
        folder =  tmpf ;    
    
        files = dir([folder '/*.txt']);
        
        file_name = files(1).name;
        x = load( [folder '/' file_name] );
        M = length(x);
        
        n = size(x,1);
        m = size(x,2);
        
        
        N = numel(files);    
        d = zeros(N,n,m);
        
        i=1;
        for file = files'
            
            file_name = file.name;

            d(i,:,:) = load( [folder '/' file_name] );
            if(mod(i,1000)==0)
                disp(i)
            end
            i=i+1;
        end

        file = [ insert 'trj' '.mat' ];
        save(file, 'd');
        
    end
end



cd(cfolder);