/* 
DESCRIPTION:  
Compute Jacobian for Decying Dimerization system in the implicit formula, 
Implicit system: 
X(t+tau) = X(t) + sum_j v_j * aj(X(t+tau))*tau + sum(..X(t)..) 

Newton-Rapson is used to find roots of: 
F  = X(t+tau) - sum_j v_j * aj(X(t+tau))*tau + B(X(t)) 

where B is precmputed and is constant in Jacobian 

Then blow is hardcoded Jacobian for F: 
J = @F1/@X1  @F1/@X2  @F1/@X3  
    @F2/@X1  @F2/@X2  @F2/@X3 
    @F3/@X1  @F3/@X2  @F3/@X3 

INPUT: 
Jacobian :  jacobian matrix stored as a vector
X        :  vector systems sate X(t)
rates    :  reaction rates
tau      :  time step

OUTPUT Jacobian: 
*/ 


#include <algorithm> 
#include <vector>    


void computeJacobian_dimerization(vector<double>& J, vector<double> X, vector<double> k, double t ){ 

std::fill (J.begin(),J.end(),0); 

J[0] =  t*(k[0]+X[0]*k[1]*4.0)+1.0;  //(1,1)  
J[1] =  t*k[2]*-2.0;  //(1,2)  

J[3] =  t*X[0]*k[1]*-2.0;  //(2,1)  
J[4] =  t*(k[2]+k[3])+1.0;  //(2,2)  

J[7] =  -t*k[3];  //(3,2)  
J[8] =  1.0;  //(3,3)  


}