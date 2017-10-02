LMGMRES := proc (f       :: procedure, 
                 J       :: procedure, 
                 H       :: procedure, 
                 X0      :: Vector, 
                 epsilon :: numeric, 
                 imax    :: posint, 
                 smax    :: positive, 
                 kmax    :: posint, 
                 lambda0 :: positive, 
                 v       :: positive, 
                 params  :: Vector) 

		   local i, iend, X, fH, Step, f_X, f_XX, Y, n, dX, lambda, M, IdMat; 

		   if _params['params'] = NULL then
		      Y := Vector(0);
		   else
		      Y := params;
		   end if;

		   n := numelems(X0);
		   X := X0;
		   lambda:=lambda0;
		   IdMat := Matrix(n,n);

		   for i to n do
		      IdMat[i,i]:=1;
		   end do;

		   dX := Vector(n);
		   f_X := <seq(10*epsilon,i=1..n)>;
		   fH := Vector(imax); 
		   Step := Vector(imax); 

		   for i to imax while Norm2(f_X) > epsilon do 
		      fH[i] := f(seq(X),seq(Y));
		      f_X := <seq(eval(J(seq(X), seq(Y))))>;
		      f_XX := convert(eval(H(seq(X), seq(Y))),Matrix);
		      M := f_XX+lambda*IdMat;
		      GMRES(M,n,f_X,dX,kmax);
		      Step[i] := Minimize1D(t->f(seq(X-t*dX), seq(Y)), 0 .. smax, epsilon, imax); 
		      if Step[i] > smax/2 then
		         lambda:=v*lambda;
		      else
		         lambda:=lambda/v;
		      end if;
		      X := X-Step[i]*dX; 
		      iend := i;
			end do; 

			fH := fH[1..iend];
			Step := Step[1..iend];
			return X, fH, Step;

end proc: