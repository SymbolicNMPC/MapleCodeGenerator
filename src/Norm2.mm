Norm2 := proc (x) 
   local i;
   sqrt(add(x[i]*x[i], i = 1 .. numelems(x))); 
end proc: