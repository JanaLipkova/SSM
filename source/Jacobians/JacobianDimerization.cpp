//
//  JacobianDimerization.cpp
//  SSM_jana_Xcode
//
//  Created by Lipkova on 14/12/17.
//  Copyright (c) 2017 Lipkova. All rights reserved.
//

#include "JacobianDimerization.h"

void computeJacobianDimerization(vector<double>& Jacobian, vector<double> X, vector<double> rates, double tau)
{
    // Decaying dimerization stochiometric coeffificents:
    vector<int> v0(3,0);
    vector<int> v1(3,0);
    vector<int> v2(3,0);
    vector<int> v3(3,0);
    
    v0[0] = -1;
    
    v1[0] = -2;
    v1[1] =  1;
    
    v2[0] =  2;
    v2[1] = -1;
    
    v3[1] = -1;
    v3[2] =  1;
    
    // jacobian inputs, where dfij = @f_{i}\@x_{j}
    
    Jacobian[0] = 1. - v0[0] * rates[0] * tau - 2. * X[0] * v1[0] * rates[1] * tau + v1[0] * rates[1] * tau;
    Jacobian[1] = - v2[0] * rates[2] * tau - v3[0] * rates[3] * tau;
    Jacobian[2] = 0.;
    
    Jacobian[3] = - v0[1] * rates[0] * tau - 2. * X[0] * v1[1] * rates[1] * tau + v1[1] * rates[1] * tau;
    Jacobian[4] = 1. - v2[1] * rates[2] * tau - v3[1] * rates[3] * tau;
    Jacobian[5] = 0.;
    
    Jacobian[6] = -v0[2] * rates[0] * tau - 2. * X[0] * v1[2] * rates[1] * tau + v1[2] * rates[1] * tau;
    Jacobian[7] = -v2[2] * rates[2] * tau - v3[2] * rates[3] * tau;
    Jacobian[8] = 1.;
    
}
