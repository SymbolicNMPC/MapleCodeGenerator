DotProduct := proc(x, y)
   local i;
   add(x[i]*y[i], i = 1 .. numelems(x));
end proc: