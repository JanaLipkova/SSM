//
//  computeJacobian.h
//  SSM_jana_Xcode
//
//  Created by Lipkova on 14/12/17.
//  Copyright (c) 2017 Lipkova. All rights reserved.
//

/*  For given reaction system computes jacobian
 - jacobian file is created by SSM/Matlab/Derivatives/
 
 
 DESCRIPTION:
 Implicit system:
 X(t+tau) = X(t) + sum_j v_j * aj(X(t+tau))*tau + sum(..X(t)..)
 
 Newton-Rapson is used to find roots of:
 F  = X(t+tau) - sum_j v_j * aj(X(t+tau))*tau + B(X(t))
 
 where B is precmputed and is constant in Jacobian
 
 Then blow is hardcoded Jacobian for F:
 J =    @F1/@X1  @F1/@X2  @F1/@X3
 @F2/@X1  @F2/@X2  @F2/@X3
 @F3/@X1  @F3/@X2  @F3/@X3
 
 INPUT:
 Jacobian :  jacobian matrix stored as a vector
 X        :  vector systems sate X(t)
 rates    :  reaction rates
 tau      : time step
 
 OUTPUT Jacobian:  matrix stored in a vector
 */

#pragma once
#include "../HeaderFiles.h"
#include "JacobianDimerization.h"
#include "JacobianLacZLacY.h"

void computeJacobian(vector<double>& Jacobian, vector<double> X, vector<double> rates, double tau );
