# Generating new variable names by appending numbers to an existing name
NewName := proc(x) 
   local params, i;
   params := [_rest];
   return cat(x, seq(seq(["_",params[i]]),i=1..numelems(params)));
end proc: