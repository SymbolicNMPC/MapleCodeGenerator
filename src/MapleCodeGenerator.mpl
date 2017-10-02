# Maple Code Generator

##
## Copyright (C) 2015-2017 Maplesoft
## Authors:    Behzad Samadi <bsamadi@maplesoft.com>
## Created:    September 2017
## Version:    0.2
## Keywords:   Maple, Code Generation, Model Predictive Control (MPC)
##
## Procedures:
##             CGMRES
##

MapleCodeGenerator := module()
description "Maple code generation for numerical optimization";
option package;

export DotProduct,
       GMRES,
       GradientDescent,
       LMGMRES,
       Minimize1D,
       NewName,
       NewtonsMethod,
       NewtonsMethodGMRES,
       nfGMRES,
       Norm2;

$include "DotProduct.mm"
$include "GMRES.mm"
$include "GradientDescent.mm"
$include "LMGMRES.mm"
$include "Minimize1D.mm"
$include "NewName.mm"
$include "NewtonsMethod.mm"
$include "NewtonsMethodGMRES.mm"
$include "nfGMRES.mm"
$include "Norm2.mm"

end module; # MCG
