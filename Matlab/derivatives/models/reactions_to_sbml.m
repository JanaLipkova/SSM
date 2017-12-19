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

for i=1:N
    model.species(i).InitialAmount = initial_pop(i);
    model.species(i).InitialAmountUnits = 'molecule';
end


save_file_name = fullfile(folder,[model_name '.xml']);

sbmlexport(model, save_file_name);

%% Change all IDs into corresponding names

text = fileread(save_file_name);

text = strrep(text,model.id,'model1');
text = strrep(text,model.Compartments.id,'Comp1');

for i=1:N
    id = model.Species(i).id;
    text = strrep(text,id,species{i});
end
for i=1:M
    id = model.Reactions(i).id;
    text = strrep(text,id,['R' ID(i) ]);
    
    id = model.Reactions(i).KineticLaw.Parameters.id;
    text = strrep(text,id,['k' ID(i) ]);
    
end




fid = fopen(save_file_name,'wt');
fprintf(fid,'%s',text);
fclose(fid);
