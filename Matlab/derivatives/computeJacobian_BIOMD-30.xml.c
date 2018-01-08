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


void computeJacobian_BIOMD-30.xml(vector<double>& J, vector<double> X, vector<double> k, double t ){ 

std::fill (J.begin(),J.end(),0); 

J[0] =  t*(h_6*unnamed*X[5]+h_9*unnamed*X[5]+k1*unnamed*X[4]+k5*unnamed*X[4])+1.0;  //(1,1)  
J[4] =  t*unnamed*(k1+k5)*X[0];  //(1,5)  
J[5] =  t*unnamed*(h_6+h_9)*X[0];  //(1,6)  
J[8] =  -k_1*t*unnamed;  //(1,9)  
J[9] =  -k_5*t*unnamed;  //(1,10)  
J[16] =  -h6*t*unnamed;  //(1,17)  
J[17] =  -h9*t*unnamed;  //(1,18)  

J[19] =  t*(h7*unnamed*X[5]+h_12*unnamed*X[5]+k3*unnamed*X[4])+1.0;  //(2,2)  
J[22] =  k3*t*unnamed*X[1];  //(2,5)  
J[23] =  t*unnamed*(h7+h_12)*X[1];  //(2,6)  
J[24] =  -k_3*t*unnamed;  //(2,7)  
J[26] =  -k2*t*unnamed;  //(2,9)  
J[30] =  -h_7*t*unnamed;  //(2,13)  
J[31] =  -h12*t*unnamed;  //(2,14)  

J[38] =  t*(h4*unnamed*X[5]+h_3*unnamed*X[5]+k7*unnamed*X[4])+1.0;  //(3,3)  
J[40] =  k7*t*unnamed*X[2];  //(3,5)  
J[41] =  t*unnamed*(h4+h_3)*X[2];  //(3,6)  
J[43] =  -k_7*t*unnamed;  //(3,8)  
J[45] =  -k6*t*unnamed;  //(3,10)  
J[50] =  -h3*t*unnamed;  //(3,15)  
J[51] =  -h_4*t*unnamed;  //(3,16)  

J[57] =  t*(h1*unnamed*X[5]+h10*unnamed*X[5])+1.0;  //(4,4)  
J[59] =  t*unnamed*(h1+h10)*X[3];  //(4,6)  
J[60] =  -k4*t*unnamed;  //(4,7)  
J[61] =  -k8*t*unnamed;  //(4,8)  
J[64] =  -h_1*t*unnamed;  //(4,11)  
J[65] =  -h_10*t*unnamed;  //(4,12)  

J[72] =  t*unnamed*(k1+k5)*X[4];  //(5,1)  
J[73] =  k3*t*unnamed*X[4];  //(5,2)  
J[74] =  k7*t*unnamed*X[4];  //(5,3)  
J[76] =  t*(k1*unnamed*X[0]+k3*unnamed*X[1]+k5*unnamed*X[0]+k7*unnamed*X[2])+1.0;  //(5,5)  
J[78] =  -t*unnamed*(k4+k_3);  //(5,7)  
J[79] =  -t*unnamed*(k8+k_7);  //(5,8)  
J[80] =  -t*unnamed*(k2+k_1);  //(5,9)  
J[81] =  -t*unnamed*(k6+k_5);  //(5,10)  

J[90] =  t*unnamed*(h_6+h_9)*X[5];  //(6,1)  
J[91] =  t*unnamed*(h7+h_12)*X[5];  //(6,2)  
J[92] =  t*unnamed*(h4+h_3)*X[5];  //(6,3)  
J[93] =  t*unnamed*(h1+h10)*X[5];  //(6,4)  
J[95] =  t*(h1*unnamed*X[3]+h4*unnamed*X[2]+h7*unnamed*X[1]+h10*unnamed*X[3]+h_3*unnamed*X[2]+h_6*unnamed*X[0]+h_9*unnamed*X[0]+h_12*unnamed*X[1])+1.0;  //(6,6)  
J[100] =  -h_1*t*unnamed;  //(6,11)  
J[101] =  -h_10*t*unnamed;  //(6,12)  
J[102] =  -h_7*t*unnamed;  //(6,13)  
J[103] =  -h12*t*unnamed;  //(6,14)  
J[104] =  -h3*t*unnamed;  //(6,15)  
J[105] =  -h_4*t*unnamed;  //(6,16)  
J[106] =  -h6*t*unnamed;  //(6,17)  
J[107] =  -h9*t*unnamed;  //(6,18)  

J[109] =  -k3*t*unnamed*X[4];  //(7,2)  
J[112] =  -k3*t*unnamed*X[1];  //(7,5)  
J[114] =  t*(k4*unnamed+k_3*unnamed)+1.0;  //(7,7)  

J[128] =  -k7*t*unnamed*X[4];  //(8,3)  
J[130] =  -k7*t*unnamed*X[2];  //(8,5)  
J[133] =  t*(k8*unnamed+k_7*unnamed)+1.0;  //(8,8)  

J[144] =  -k1*t*unnamed*X[4];  //(9,1)  
J[148] =  -k1*t*unnamed*X[0];  //(9,5)  
J[152] =  t*(k2*unnamed+k_1*unnamed)+1.0;  //(9,9)  

J[162] =  -k5*t*unnamed*X[4];  //(10,1)  
J[166] =  -k5*t*unnamed*X[0];  //(10,5)  
J[171] =  t*(k6*unnamed+k_5*unnamed)+1.0;  //(10,10)  

J[183] =  -h1*t*unnamed*X[5];  //(11,4)  
J[185] =  -h1*t*unnamed*X[3];  //(11,6)  
J[190] =  t*(h2*unnamed+h_1*unnamed)+1.0;  //(11,11)  

J[201] =  -h10*t*unnamed*X[5];  //(12,4)  
J[203] =  -h10*t*unnamed*X[3];  //(12,6)  
J[209] =  t*(h11*unnamed+h_10*unnamed)+1.0;  //(12,12)  

J[217] =  -h7*t*unnamed*X[5];  //(13,2)  
J[221] =  -h7*t*unnamed*X[1];  //(13,6)  
J[228] =  t*(h8*unnamed+h_7*unnamed)+1.0;  //(13,13)  

J[235] =  -h_12*t*unnamed*X[5];  //(14,2)  
J[239] =  -h_12*t*unnamed*X[1];  //(14,6)  
J[245] =  -h11*t*unnamed;  //(14,12)  
J[247] =  h12*t*unnamed+1.0;  //(14,14)  

J[254] =  -h_3*t*unnamed*X[5];  //(15,3)  
J[257] =  -h_3*t*unnamed*X[2];  //(15,6)  
J[262] =  -h2*t*unnamed;  //(15,11)  
J[266] =  h3*t*unnamed+1.0;  //(15,15)  

J[272] =  -h4*t*unnamed*X[5];  //(16,3)  
J[275] =  -h4*t*unnamed*X[2];  //(16,6)  
J[285] =  t*(h5*unnamed+h_4*unnamed)+1.0;  //(16,16)  

J[288] =  -h_6*t*unnamed*X[5];  //(17,1)  
J[293] =  -h_6*t*unnamed*X[0];  //(17,6)  
J[303] =  -h5*t*unnamed;  //(17,16)  
J[304] =  h6*t*unnamed+1.0;  //(17,17)  

J[306] =  -h_9*t*unnamed*X[5];  //(18,1)  
J[311] =  -h_9*t*unnamed*X[0];  //(18,6)  
J[318] =  -h8*t*unnamed;  //(18,13)  
J[323] =  h9*t*unnamed+1.0;  //(18,18)  


}