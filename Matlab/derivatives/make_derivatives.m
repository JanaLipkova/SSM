clear

%% load the model
% run dimerization.m
run lacy_lacz.m


M = length(reaction);
N = length(species);

n=max( [ fix(abs(log10(abs(M))))+1  fix(abs(log10(abs(N))))+1 ] );
ID = @(i) sprintf(['%0' num2str(n) 'd'],i);

%% make a simbio model object with the original species name
model = sbiomodel(model_name);
for i=1:M
    r_obj  = addreaction(model, reaction{i});
    kl_obj = addkineticlaw(r_obj, 'MassAction');
    set( kl_obj, 'ParameterVariablenames', ['k' ID(i)] );
    p_obj  = addparameter(kl_obj, ['k' ID(i)], rate(i));    
end
sbmlexport(model, [model_name '.xml']);



%% rename species to X01, X02, ...


for i=1:N
    model.Species(i).rename([ ' X' ID(i) ] );
end



%% evaluate the Jacobian of the propensities

% for 2017b use str2sym and remove the suppress of the warning 
warning('off', 'symbolic:sym:sym:DeprecateExpressions' );

prop = cell(M,1);
for i=1:M
    prop{i} = sym( get(model.Reactions(i), 'ReactionRate') );
end

var = cell(N,1);
for i=1:N
    var{i} = sym([ 'X' ID(i) ]);
end


Jac = cell(M,N);
for i=1:M
    for j=1:N
        Jac{i,j} = diff(prop{i},var{j});
    end
end



%% evaluate the auxialiary matrix f

S = getstoichmatrix(model);

f = cell(M,M);
syms tmp
for i=1:M
    for j=1:M
        tmp = 0;
        for k = 1:N
            tmp = tmp + Jac{i,k}*S(k,j);
        end
        f{i,j} = simplify(tmp);
    end
end





%% evaluate mu and sigma2

mu = cell(M,1);
sigma2 = cell(M,1);
syms tmp1 tmp2
for i=1:M
    tmp1 = 0;
    tmp2 = 0;
    for k=1:M
        tmp1 = tmp1 + f{i,k}*prop{k};
        tmp2 = tmp2 + f{i,k}^2*prop{k};
    end
    mu{i}     = simplify(tmp1);
    sigma2{i} = simplify(tmp2);
end




%% substitute variables in C format e.g., X1-->X[0]
mu_s = cell(M,1);
sigma2_s = cell(M,1);

for i=1:M
    mu_s{i} = char(mu{i});
    sigma2_s{i} = char(sigma2{i});
end


for i=1:M
    for j=1:N
        str_old = [ 'X' ID(j)  ];
        str_new = [ 'X[' num2str(j) ']' ];

        mu_s{i} = strrep( mu_s{i}, str_old, str_new );
        sigma2_s{i} = strrep( sigma2_s{i}, str_old, str_new );
        
    end
    
    for j=1:M
        str_old = [ 'k' ID(j) ];
        str_new = [ 'k[' num2str(j) ']' ];

        mu_s{i} = strrep( mu_s{i}, str_old, str_new );
        sigma2_s{i} = strrep( sigma2_s{i}, str_old, str_new );
    end
end


%% write function in C format


cname = [ 'mu_sigma_' model_name '.c'];
if( exist(cname,'file') )
    delete(cname);
end


fp = fopen( cname, 'w');

fprintf(fp,'//Output: \n');
fprintf(fp,'//        mu, sigma2 \n');
fprintf(fp,'//Input: \n');
fprintf(fp,'//       X: species population \n');
fprintf(fp,'//       k: reaction constant for each reaction \n');

fprintf(fp,'void mu_sigma(double *mu, double *sigma2, double *X, double *k){ \n\n');

for i=1:M
    t = ccode(sym(mu_s{i})); t(1:6)=[];
    fprintf(fp,'mu[%d] = %s \n', i-1, t );
end

fprintf(fp,'\n\n');


for i=1:M
    t = ccode(sym(sigma2_s{i})); t(1:6)=[];
    fprintf(fp,'sigma2[%d] = %s \n', i-1, t );
end



fprintf(fp,'} \n\n');

fclose(fp);














