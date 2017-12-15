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


void computeJacobian_egfr(vector<double>& J, vector<double> X, vector<double> k, double t ){ 

std::fill (J.begin(),J.end(),0); 

J[0] =  t*(compartment*k1+compartment*k2)+1.0;  //(1,1)  
J[1] =  -compartment*k1*t*X[25];  //(1,2)  
J[25] =  -compartment*k1*t*X[1];  //(1,26)  

J[31] =  -compartment*k2*t;  //(2,1)  
J[32] =  t*(compartment*k1+compartment*k1*X[25])+1.0;  //(2,2)  
J[56] =  compartment*k1*t*X[1];  //(2,26)  

J[62] =  -(Kcat*compartment*t*X[3])/(km+X[3]);  //(3,1)  
J[64] =  compartment*t*1.0/pow(km+X[2],2.0)*(k1*(km*km)+k1*(X[2]*X[2])+Kcat*km*X[12]+k1*km*X[2]*2.0)+1.0;  //(3,3)  
J[65] =  -Kcat*compartment*km*t*1.0/pow(km+X[3],2.0)*X[0];  //(3,4)  
J[74] =  (Kcat*compartment*t*X[2])/(km+X[2]);  //(3,13)  

J[93] =  (Kcat*compartment*t*X[3])/(km+X[3]);  //(4,1)  
J[95] =  -compartment*t*1.0/pow(km+X[2],2.0)*(k1*(km*km)+k1*(X[2]*X[2])+Kcat*km*X[12]+k1*km*X[2]*2.0);  //(4,3)  
J[96] =  1.0/pow(km+X[3],2.0)*(X[3]*X[3]+km*X[3]*2.0+km*km+Kcat*compartment*km*t*X[0]);  //(4,4)  
J[105] =  -(Kcat*compartment*t*X[2])/(km+X[2]);  //(4,13)  

J[126] =  -(Kcat*compartment*t*X[5])/(km+X[5]);  //(5,3)  
J[128] =  1.0/pow(km+X[4],2.0)*(X[4]*X[4]+km*X[4]*2.0+km*km+Kcat*compartment*km*t*X[28]);  //(5,5)  
J[129] =  -Kcat*compartment*km*t*1.0/pow(km+X[5],2.0)*X[2];  //(5,6)  
J[152] =  (Kcat*compartment*t*X[4])/(km+X[4]);  //(5,29)  

J[157] =  (Kcat*compartment*t*X[5])/(km+X[5]);  //(6,3)  
J[159] =  -Kcat*compartment*km*t*1.0/pow(km+X[4],2.0)*X[28];  //(6,5)  
J[160] =  1.0/pow(km+X[5],2.0)*(X[5]*X[5]+km*X[5]*2.0+km*km+Kcat*compartment*km*t*X[2]);  //(6,6)  
J[183] =  -(Kcat*compartment*t*X[4])/(km+X[4]);  //(6,29)  

J[190] =  -(Kcat*compartment*t*X[7])/(km+X[7]);  //(7,5)  
J[192] =  1.0/pow(km+X[6],2.0)*(X[6]*X[6]+km*X[6]*2.0+km*km+Kcat*compartment*km*t*X[16]+Kcat*compartment*km*t*X[27]);  //(7,7)  
J[193] =  -Kcat*compartment*km*t*1.0/pow(km+X[7],2.0)*X[4];  //(7,8)  
J[202] =  (Kcat*compartment*t*X[6])/(km+X[6]);  //(7,17)  
J[213] =  (Kcat*compartment*t*X[6])/(km+X[6]);  //(7,28)  

J[221] =  (Kcat*compartment*t*X[7])/(km+X[7]);  //(8,5)  
J[223] =  -Kcat*compartment*km*t*1.0/pow(km+X[6],2.0)*(X[16]+X[27]);  //(8,7)  
J[224] =  1.0/pow(km+X[7],2.0)*(X[7]*X[7]+km*X[7]*2.0+km*km+Kcat*compartment*km*t*X[4]);  //(8,8)  
J[233] =  -(Kcat*compartment*t*X[6])/(km+X[6]);  //(8,17)  
J[244] =  -(Kcat*compartment*t*X[6])/(km+X[6]);  //(8,28)  

J[254] =  -(Kcat*compartment*t*X[9])/(km+X[9]);  //(9,7)  
J[256] =  1.0/pow(km+X[8],2.0)*(X[8]*X[8]+km*X[8]*2.0+km*km+Kcat*compartment*km*t*X[26]);  //(9,9)  
J[257] =  -compartment*km*t*(Kcat*X[6]+kcat*X[23])*1.0/pow(km+X[9],2.0);  //(9,10)  
J[271] =  -(compartment*kcat*t*X[9])/(km+X[9]);  //(9,24)  
J[274] =  (Kcat*compartment*t*X[8])/(km+X[8]);  //(9,27)  

J[285] =  (Kcat*compartment*t*X[9])/(km+X[9]);  //(10,7)  
J[287] =  -Kcat*compartment*km*t*1.0/pow(km+X[8],2.0)*X[26];  //(10,9)  
J[288] =  1.0/pow(km+X[9],2.0)*(X[9]*X[9]+km*X[9]*2.0+km*km+Kcat*compartment*km*t*X[6]+compartment*kcat*km*t*X[23]);  //(10,10)  
J[302] =  (compartment*kcat*t*X[9])/(km+X[9]);  //(10,24)  
J[305] =  -(Kcat*compartment*t*X[8])/(km+X[8]);  //(10,27)  

J[318] =  -(Kcat*compartment*t*X[11])/(km+X[11]);  //(11,9)  
J[320] =  1.0/pow(km+X[10],2.0)*(X[10]*X[10]+km*X[10]*2.0+km*km+Kcat*compartment*km*t*X[26]);  //(11,11)  
J[321] =  -Kcat*compartment*km*t*1.0/pow(km+X[11],2.0)*X[8];  //(11,12)  
J[336] =  (Kcat*compartment*t*X[10])/(km+X[10]);  //(11,27)  

J[349] =  (Kcat*compartment*t*X[11])/(km+X[11]);  //(12,9)  
J[351] =  -Kcat*compartment*km*t*1.0/pow(km+X[10],2.0)*X[26];  //(12,11)  
J[352] =  1.0/pow(km+X[11],2.0)*(X[11]*X[11]+km*X[11]*2.0+km*km+Kcat*compartment*km*t*X[8]);  //(12,12)  
J[367] =  -(Kcat*compartment*t*X[10])/(km+X[10]);  //(12,27)  

J[382] =  -(Kcat*compartment*t*X[13])/(km+X[13]);  //(13,11)  
J[384] =  compartment*k1*t+1.0;  //(13,13)  
J[385] =  -Kcat*compartment*km*t*1.0/pow(km+X[13],2.0)*X[10];  //(13,14)  

J[413] =  (Kcat*compartment*t*X[13])/(km+X[13]);  //(14,11)  
J[415] =  -compartment*k1*t;  //(14,13)  
J[416] =  1.0/pow(km+X[13],2.0)*(X[13]*X[13]+km*X[13]*2.0+km*km+Kcat*compartment*km*t*X[10]);  //(14,14)  

J[434] =  -(Kcat*compartment*t*X[15])/(km+X[15]);  //(15,1)  
J[438] =  -(Kcat*compartment*t*X[15])/(km+X[15]);  //(15,5)  
J[448] =  compartment*k1*t+1.0;  //(15,15)  
J[449] =  -Kcat*compartment*km*t*1.0/pow(km+X[15],2.0)*(X[0]+X[4]);  //(15,16)  

J[465] =  (Kcat*compartment*t*X[15])/(km+X[15]);  //(16,1)  
J[469] =  (Kcat*compartment*t*X[15])/(km+X[15]);  //(16,5)  
J[479] =  -compartment*k1*t;  //(16,15)  
J[480] =  1.0/pow(km+X[15],2.0)*(X[15]*X[15]+km*X[15]*2.0+km*km+Kcat*compartment*km*t*X[0]+Kcat*compartment*km*t*X[4]);  //(16,16)  

J[510] =  -(Kcat*compartment*t*X[17])/(km+X[17]);  //(17,15)  
J[512] =  compartment*k1*t+1.0;  //(17,17)  
J[513] =  -Kcat*compartment*km*t*1.0/pow(km+X[17],2.0)*X[14];  //(17,18)  

J[541] =  (Kcat*compartment*t*X[17])/(km+X[17]);  //(18,15)  
J[543] =  -compartment*k1*t;  //(18,17)  
J[544] =  1.0/pow(km+X[17],2.0)*(X[17]*X[17]+km*X[17]*2.0+km*km+Kcat*compartment*km*t*X[14]);  //(18,18)  

J[558] =  -compartment*k1*t;  //(19,1)  
J[559] =  -compartment*k1*t;  //(19,2)  
J[576] =  1.0;  //(19,19)  

J[589] =  -(compartment*kcat*t*X[20])/(km+X[20]);  //(20,1)  
J[608] =  compartment*k1*t+1.0;  //(20,20)  
J[609] =  -compartment*kcat*km*t*1.0/pow(km+X[20],2.0)*X[0];  //(20,21)  

J[620] =  (compartment*kcat*t*X[20])/(km+X[20]);  //(21,1)  
J[639] =  -compartment*k1*t;  //(21,20)  
J[640] =  1.0/pow(km+X[20],2.0)*(X[20]*X[20]+km*X[20]*2.0+km*km+compartment*kcat*km*t*X[0]);  //(21,21)  

J[670] =  -(compartment*kcat*t*X[22])/(km+X[22]);  //(22,20)  
J[672] =  1.0/pow(km+X[21],2.0)*(X[21]*X[21]+km*X[21]*2.0+km*km+compartment*kcat*km*t*X[29]);  //(22,22)  
J[673] =  -compartment*kcat*km*t*1.0/pow(km+X[22],2.0)*X[19];  //(22,23)  
J[680] =  (compartment*kcat*t*X[21])/(km+X[21]);  //(22,30)  

J[701] =  (compartment*kcat*t*X[22])/(km+X[22]);  //(23,20)  
J[703] =  -compartment*kcat*km*t*1.0/pow(km+X[21],2.0)*X[29];  //(23,22)  
J[704] =  1.0/pow(km+X[22],2.0)*(X[22]*X[22]+km*X[22]*2.0+km*km+compartment*kcat*km*t*X[19]);  //(23,23)  
J[711] =  -(compartment*kcat*t*X[21])/(km+X[21]);  //(23,30)  

J[717] =  -(compartment*kcat*t*X[24])/(km+X[24]);  //(24,5)  
J[734] =  -(compartment*kcat*t*X[24])/(km+X[24]);  //(24,22)  
J[736] =  1.0/pow(km+X[23],2.0)*(X[23]*X[23]+km*X[23]*2.0+km*km+compartment*kcat*km*t*X[27]);  //(24,24)  
J[737] =  -compartment*kcat*km*t*1.0/pow(km+X[24],2.0)*(X[4]+X[21]);  //(24,25)  
J[740] =  (compartment*kcat*t*X[23])/(km+X[23]);  //(24,28)  

J[748] =  (compartment*kcat*t*X[24])/(km+X[24]);  //(25,5)  
J[765] =  (compartment*kcat*t*X[24])/(km+X[24]);  //(25,22)  
J[767] =  -compartment*kcat*km*t*1.0/pow(km+X[23],2.0)*X[27];  //(25,24)  
J[768] =  1.0/pow(km+X[24],2.0)*(X[24]*X[24]+km*X[24]*2.0+km*km+compartment*kcat*km*t*X[4]+compartment*kcat*km*t*X[21]);  //(25,25)  
J[771] =  -(compartment*kcat*t*X[23])/(km+X[23]);  //(25,28)  

J[775] =  -compartment*k2*t;  //(26,1)  
J[776] =  compartment*k1*t*X[25];  //(26,2)  
J[800] =  compartment*k1*t*X[1]+1.0;  //(26,26)  

J[832] =  1.0;  //(27,27)  

J[864] =  1.0;  //(28,28)  

J[896] =  1.0;  //(29,29)  

J[928] =  1.0;  //(30,30)  

J[960] =  1.0;  //(31,31)  


}