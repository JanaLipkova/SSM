//
//  JacobianLacZLacy.h
//  SSM_jana_Xcode
//
//  Created by Lipkova on 14/12/17.
//  Copyright (c) 2017 Lipkova. All rights reserved.
//

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

#pragma once
#include "../HeaderFiles.h"

void computeJacobianLacZLacY(vector<double>& J, vector<double> X, vector<double> k, double t );
    