function U  = velocity_monkewitz(y,n, Lambda)
F = 1./(1 + (sinh(y./sinh(1))).^(2*n));
U = 1 - Lambda+2*Lambda*F;
end