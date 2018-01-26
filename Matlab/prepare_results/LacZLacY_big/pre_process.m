clear
system = 'lacy_lacz2';
method = {'SSA','TauLeap','AdaptiveTau','RLeaping','SLeaping_v3','SLeaping_v4','SLeaping_v5','SLeaping_v6'};
% method = {};

% eps = { '0.05',  '0.03',  '0.01'};
eps = { '1.0', '0.5', '0.1' };

cfolder = pwd;
%%
cd('SSA');
folders = dir('*Trajectories');
data = [];

for k=1:length(folders)

    folder =  [ system '_' num2str(k) '_Trajectories' ] ;    
    
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

    
    
    file = [ 'trj_' num2str(k) '.mat' ];
    save(file, 'd'); 
    
end




%%
cd(cfolder);
for k=1:length(method)
    
    if( strcmp(method{k},'SSA') )
        continue;
    end
    
    
    K = length(eps);
    
    for l=1:K
    
        cd(cfolder);
        
        cd(method{k})
        
        fprintf('%s --  %s \n',method{k},eps{l});
       
        
        insert = [ 'eps_' eps{l} '_' ];
         
        tmpf = [system '_' method{k} '_' insert 'Trajectories'];
%         tmpf = [system '_' insert 'Trajectories'];
        
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