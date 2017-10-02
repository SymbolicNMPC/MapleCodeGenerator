GradientDescent := proc (f::procedure, 
                         J::procedure, 
						 X0::Vector, 
						 epsilon::numeric, 
						 imax::posint, 
						 smax::positive, 
						 params::Vector) 

				   local i, iend, X, fH, Step, f_X, Y, n; 

				   if _params['params'] = NULL then
				      Y := Vector(0);
				   else
				      Y := params;
				   end if;

				   n := numelems(X0);
				   X := X0;
				   f_X := <seq(10*epsilon,i=1..n)>;
				   fH := Vector(imax); 
				   Step := Vector(imax); 

				   for i to imax while Norm2(f_X) > epsilon do 
				      fH[i] := f(seq(X),seq(Y));
				      f_X := <seq(eval(J(seq(X), seq(Y))))>;
				      Step[i] := Minimize1D(t->f(seq(X-t*f_X), seq(Y)), 0 .. smax, epsilon, imax); 
				      X := X-Step[i]*f_X; 
				      iend := i;
				   end do; 

				   fH := fH[1..iend];
				   Step := Step[1..iend];

				   return X, fH, Step;
end proc:
