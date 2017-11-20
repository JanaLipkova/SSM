//Output: 
//        mu, sigma2 
//Input: 
//       X: species population 
//       k: reaction constant for each reaction 
void mu_sigma(double *mu, double *sigma2, double *X, double *k){ 

mu[0] = -k[0]*(X[0]*k[0] - 2*X[1]*k[2] + 2*X[0]^2*k[1]); 
mu[1] = -2*X[0]*k[1]*(X[0]*k[0] - 2*X[1]*k[2] + 2*X[0]^2*k[1]); 
mu[2] = -k[2]*(X[1]*k[2] + X[1]*k[3] - X[0]^2*k[1]); 
mu[3] = -k[3]*(X[1]*k[2] + X[1]*k[3] - X[0]^2*k[1]); 


sigma2[0] = k[0]^2*(X[0]*k[0] + 4*X[1]*k[2] + 4*X[0]^2*k[1]); 
sigma2[1] = 4*X[0]^2*k[1]^2*(X[0]*k[0] + 4*X[1]*k[2] + 4*X[0]^2*k[1]); 
sigma2[2] = k[2]^2*(X[1]*k[2] + X[1]*k[3] + X[0]^2*k[1]); 
sigma2[3] = k[3]^2*(X[1]*k[2] + X[1]*k[3] + X[0]^2*k[1]); 
} 

