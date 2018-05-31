function y = h_function(x, ksai)
v = 10;
s = 5;
l = 4;
c = 3;
rho = 40;
y = - v*min(x,ksai) - s*max(x-ksai,0) + l*max(ksai-x,0) + c*x + rho;
end