load "Lseries6.magma";
Q:=qEigenform(ModularSymbols("37A"),1000);
V:=form_curve(37,-2,-3,-2,-1);
run_thru_primes(V,Q,251,1000);

