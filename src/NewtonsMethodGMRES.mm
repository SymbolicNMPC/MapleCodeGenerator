NewtonsMethodGMRES := proc (f      :: procedure, 
                            J      :: procedure, 
							H      :: procedure, 
							X0     :: Vector, 
							epsilon:: numeric, 
							imax   :: posint, 
							smax   :: positive, 
							kmax   :: posint, 
							params :: Vector) 

   local i, iend, X, fH, Step, f_X, f_XX, Y, n, dX; 

   if _params['params'] = NULL then
      Y := Vector(0);
   else
      Y := params;
   end if;

   n := numelems(X0);
   X := X0;
   dX := Vector(n);
   f_X := <seq(10*epsilon,i=1..n)>;
   fH := Vector(imax); 
   Step := Vector(imax); 

   for i to imax while Norm2(f_X) > epsilon do 
      fH[i] := f(seq(X),seq(Y));
      f_X := <seq(eval(J(seq(X), seq(Y))))>;
      f_XX := convert(eval(H(seq(X), seq(Y))),Matrix);
      GMRES(f_XX,n,f_X,dX,kmax);
      Step[i] := Minimize1D(t->f(seq(X-t*dX), seq(Y)), 0 .. smax, epsilon, imax); 
      X := X-Step[i]*dX; 
      iend := i;
   end do; 

   fH := fH[1..iend];
   Step := Step[1..iend];

   return X, fH, Step;
   
end proc: