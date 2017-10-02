Minimize1D := proc(f::procedure, 
                   R::range, 
				   epsilon::positive, 
				   imax::posint)

		    local phi, a, b, c, i, f1, f2;

		    a:=lhs(R);
		    b:=rhs(R);
		    phi := (sqrt(5.) - 1)/2:

            i := 0;
		    while abs(b - a) > epsilon and i < imax do
		       c[1] := phi*a + (1 - phi)*b;
		       c[2] := (1 - phi)*a + phi*b;
		       f1 := f(c[1]);
		       f2 := f(c[2]);      
		       if Re(f2)=f2 then
		           if f1 < f2 then
		              b := c[2];
		              c[2] := c[1];
		              c[1] := phi*a + (1 - phi)*b;
		           else
		              a := c[1];
		              c[1] := c[2];
		              c[2] := (1 - phi)*a + phi*b;
		           end if;                 
		        else # avoiding complex numbers
		             # assumption: only large numbers 
		             #             can lead to complex numbers
		           b := c[2];
		        end if;
				i := i+1;
		    end do:  

		    return a;
end proc: