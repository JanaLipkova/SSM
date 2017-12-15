clear

folder = 'sbml';

% file   = 'dimerization.m';
file   = 'lacy_lacz.m';


run(fullfile(folder,file))
%%

M = length(reaction);
N = length(species);


%% make a simbio model object with the original species name
n=max( [ fix(abs(log10(abs(M))))+1  fix(abs(log10(abs(N))))+1 ] );
ID = @(i) sprintf(['%0' num2str(n) 'd'],i);

model = sbiomodel(model_name);
for i=1:M
    r_obj  = addreaction(model, reaction{i});
    kl_obj = addkineticlaw(r_obj, 'MassAction');
    set( kl_obj, 'ParameterVariablenames', ['k' ID(i)] );
    p_obj  = addparameter(kl_obj, ['k' ID(i)], rate(i));    
end

save_file_name = fullfile(folder,[model_name '.xml']);

sbmlexport(model, save_file_name);