##/*--------------------------------------------------------
## Linear Equation Subroutine
##	nfGMRES (Generalized Minimum Residual) 
##	Cf. C.T.Kelley: Iterative Methods for Linear and Nonlinear Equations
##						T.Ohtsuka  '00/01/24
##	axfunc(int n, double *x, double *ax): a function that gives A*x
##	n: dim x
##	b: the right-hand-side of the equation, A*x = b
##	x: initial guess, to be replaced with the solution
##	kmax: number of iteration
##	err: residual norms, |b - A*x_i| (i=1,...,n+1)
##----------------------------------------------------------*/
##
## Translated to Maple by Behzad Samadi, Aug. 2016
## Note that err is a vector of size kmax+1

nfGMRES := proc(axfunc, n, b, x, kmax, err)

	local i,j,k;
	local rho, nu; 
	local cvec, svec, gvec, tmpvec, hmat, vmat; 
	local temp;

	cvec := Vector(kmax+1,1); 
	svec := Vector(kmax+1); 
	gvec := Vector(kmax+1);
	tmpvec := Vector(n); 
	hmat := Matrix(kmax+1,kmax+1); 
	vmat := Matrix(n, kmax+1); 

	axfunc(x, tmpvec, _rest); 

	tmpvec := b-tmpvec;

	rho := sqrt(add(tmpvec[i]*tmpvec[i], i=1..n));

	gvec[1] := rho; 

    for i from 2 to kmax+1 do
      gvec[i]:=0;
    end do;

	err[1] := rho;
	
    for i to n do
      vmat[i,1]:=tmpvec[i]/rho;
    end do;

	
	for k to kmax do 
	
		axfunc(vmat[..,k], tmpvec, _rest); 
		vmat[..,k+1] := tmpvec;

		## /* Modified Gram-Schmidt */
		for j to k do
		   hmat[k,j] := add(vmat[i,j]*vmat[i,k+1], i=1..n);
		   vmat[..,k+1] := vmat[..,k+1]-hmat[k,j]*vmat[..,j];
		end do;
		
		hmat[k,k+1] := sqrt(add(vmat[i,k+1]*vmat[i,k+1], i=1..n));
		
		## /* No Breakdown? */
		if hmat[k,k+1] > 0.0 then
		   vmat[..,k+1]:=vmat[..,k+1]/hmat[k,k+1];
		else
        ##  error("gmres() : breakdown \n");
            return;
		end if;
		
		## /* Givens Rotation */
        for j from 1 to k-1 do
       	   temp := cvec[j]*hmat[k,j]-svec[j]*hmat[k,j+1];
      	   hmat[k,j+1] := svec[j]*hmat[k,j]+cvec[j]*hmat[k,j+1];
      	   hmat[k,j] := temp;
        end do;
		
		nu := sqrt(hmat[k,k] * hmat[k,k] + hmat[k,k+1] * hmat[k,k+1]);
		if nu > 0 then
			cvec[k] := hmat[k,k] / nu;
			svec[k] := - hmat[k,k+1] / nu;
			hmat[k,k] := cvec[k] * hmat[k,k] - svec[k] * hmat[k,k+1];
			hmat[k,k+1] := 0;
			
			## givrot(1, cvec+k, svec+k, gvec+k);
			temp := cvec[k] * gvec[k] - svec[k] * gvec[k+1];
		    gvec[k+1] := svec[k] * gvec[k] + cvec[k] * gvec[k+1];
			gvec[k] := temp;
						
        ## else 
        ##    WARNING("nu is zero!\n");
		end if;
		
		## /* Residual Update */
		rho := abs(gvec[k+1]);
		err[k+1] := rho;
	end do;
	
	## /* Solve hmat * y = gvec (hmat: upper triangular) */
	for i from kmax to 1 by -1 do
	   nu := gvec[i];
	   for j from i+1 to kmax do
	      nu := nu - hmat[j,i] * cvec[j];
	   end do;
	   cvec[i] := nu/hmat[i,i];
	end do;
	
	## /* Answer */
	for i to n do
	   nu := 0;
	   for j to kmax do
	      nu := nu + vmat[i,j] * cvec[j];
	   end do;
	   x[i] := x[i] + nu;
	end do;
		
end proc;