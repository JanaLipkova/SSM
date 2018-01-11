clc; clear;

sbioreset

abstkineticlawObj = sbioabstractkineticlaw('ex_mylaw1', '(k1*s)/(k2+k1+s)');

set (abstkineticlawObj, 'SpeciesVariables', {'s'});
set (abstkineticlawObj, 'ParameterVariables', {'k1', 'k2'});


sbioaddtolibrary(abstkineticlawObj);

sbiowhos -kineticlaw -userdefined
%%

model = sbiomodel('cell');
r_obj = addreaction(model, 'A + B -> B + C');
kl_obj = addkineticlaw(r_obj, 'ex_mylaw1');
set(kl_obj, 'ParameterVariablenames', {'k1','k2'} );
set(kl_obj,'SpeciesVariableNames', {'A'});

model.Reactions.ReactionRate
%%
model = sbiomodel('test');
r_obj  = addreaction(model, '2 A -> B + C' );
kl_obj = addkineticlaw(r_obj, 'MassAction');
set( kl_obj, 'ParameterVariablenames', 'k1' );
cs = getconfigset(model,'active');
cs.SolverType = 'ssa';
model.Reactions.ReactionRate

%%
model = sbiomodel('cell');
r_obj  = addreaction(model, 'A + B -> B + C' );
kl_obj = addkineticlaw(r_obj, 'Henri-Michaelis-Menten' );
set(kl_obj,'ParameterVariableNames', {'V' 'K'});
set(kl_obj,'SpeciesVariableNames', {'A'});


model.Reactions.ReactionRate




%%



