clc; clear


% path_to_model = 'models/sbml';
% model_name = 'dimerization';

% path_to_model = '../../ReactionSystemsXML/LotkaVoltera/';
% model_name = 'LotkaVolter_R';

path_to_model = 'models/sbml/';
model_name = 'BIOMD-30.xml';




%%
model = sbmlimport( fullfile(path_to_model,model_name) );


N = length(model.Species);
M = length(model.Reactions);

%% rename species to X01, X02, ...

n=max( [ fix(abs(log10(abs(M))))+1  fix(abs(log10(abs(N))))+1 ] );
ID = @(i) sprintf(['%0' num2str(n) 'd'],i);

for i=1:N
    model.Species(i).rename([ ' X' ID(i) ] );
end



%% evaluate the Jacobian of the propensities

% for 2017b use str2sym and remove the suppress of the warning 
warning('off', 'symbolic:sym:sym:DeprecateExpressions' );

prop = cell(M,1);
for i=1:M
    s  = get(model.Reactions(i), 'ReactionRate');
    sn = strrep(s,'power',''); sn = strrep(sn,',',')^(');   % substitute power(X,2) with (X)^(2)
    prop{i} = sym( sn );
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







%% evaluate the RHS

S = getstoichmatrix(model);

t=sym('t');

F = cell(N);
for i=1:N
    
    tmp = sym(0);
    for j=1:M
        tmp = tmp + S(i,j)*prop{j};
    end
    F{i} = var{i} - t*tmp;
    
end


%% evaluate the Jacobian of RHS

J = cell(N,N);
for i=1:N
    for j=1:N
        tmp = diff(F{i},var{j});
        J{i,j} = simplify(tmp);
    end
    
end





%% substitute variables in C format e.g., X1-->X[0]
J_s = cell(N,M);

for i=1:N
    for j=1:N
        J_s{i,j} = char(J{i,j});
    end
end

for i=1:N
    for j=1:N
        
        for k=1:N
            str_old = [ 'X' ID(k)  ];
            str_new = [ 'X[' num2str(k) ']' ];

            J_s{i,j} = strrep( J_s{i,j}, str_old, str_new );
        end
        
        for k=1:M
            str_old = [ 'k' ID(k) ];
            str_new = [ 'k[' num2str(k) ']' ];

            J_s{i,j} = strrep( J_s{i,j}, str_old, str_new );
        end
    end
end


%% write function in C format


cname = [ 'computeJacobian_' model_name '.c'];
if( exist(cname,'file') )
    delete(cname);
end


fp = fopen( cname, 'w');

fprintf(fp,'/* \n');
fprintf(fp,'DESCRIPTION:  \n');
fprintf(fp,'Compute Jacobian for Decying Dimerization system in the implicit formula, \n');
fprintf(fp,'Implicit system: \n');
fprintf(fp,'X(t+tau) = X(t) + sum_j v_j * aj(X(t+tau))*tau + sum(..X(t)..) \n\n');

fprintf(fp,'Newton-Rapson is used to find roots of: \n');
fprintf(fp,'F  = X(t+tau) - sum_j v_j * aj(X(t+tau))*tau + B(X(t)) \n\n');
 
fprintf(fp,'where B is precmputed and is constant in Jacobian \n\n');
 
fprintf(fp,'Then blow is hardcoded Jacobian for F: \n');
fprintf(fp,'J = @F1/@X1  @F1/@X2  @F1/@X3  \n');
fprintf(fp,'    @F2/@X1  @F2/@X2  @F2/@X3 \n');
fprintf(fp,'    @F3/@X1  @F3/@X2  @F3/@X3 \n\n');

 fprintf(fp,'INPUT: \n');
 fprintf(fp,'Jacobian :  jacobian matrix stored as a vector\n');
 fprintf(fp,'X        :  vector systems sate X(t)\n');
 fprintf(fp,'rates    :  reaction rates\n');
 fprintf(fp,'tau      :  time step\n\n');
 
 fprintf(fp,'OUTPUT Jacobian: \n');
 fprintf(fp,'*/ \n\n\n');

 
 fprintf(fp,'#include <algorithm> \n');
 fprintf(fp,'#include <vector>    \n\n\n');
 
 fprintf(fp,'void computeJacobian_%s(vector<double>& J, vector<double> X, vector<double> k, double t ){ \n\n',model_name);
 
 % fill the vector with zeros
 fprintf(fp,'std::fill (J.begin(),J.end(),0); \n\n');
 
 

for i=1:N
    for j=1:N
        
        if ~isequal(J{i,j},sym(0))
            t = ccode(sym(J_s{i,j})); t(1:6)=[];
            fprintf(fp,'J[%d] = %s  //(%d,%d)  \n', (i-1)*N + j-1 , t,i,j);
        end
        
    end
    fprintf(fp,'\n');    
end

fprintf(fp,'\n');



fprintf(fp,'}');

fclose(fp);










