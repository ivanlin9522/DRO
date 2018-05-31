function [ result ] = second_dev_cov_int( v,s,l,c,rho,psy_vector,n )
% [ result ] = second_dev_cov_int( v,s,l,c,rho,psy_vector,n )
% Note that psy_vector need to be vector with n elements. n must greater or equals to 2
%

v_str=num2str(v);
s_str=num2str(s);
l_str=num2str(l);
c_str=num2str(c);
n_str=num2str(n);
rho_str=num2str(rho);
psy_vector_str=num2str(psy_vector);

v_command=['v:=' v_str];
s_command=['s:=' s_str];
l_command=['l:=' l_str];
c_command=['c:=' c_str];
rho_command=['rho:=' rho_str];
n_command=['ndata:=' n_str];

psy_command=['etas={' psy_vector_str(1,:)];
for i=2:n
    psy_command=[psy_command ',' psy_vector_str(i,:)];
end
psy_command=[psy_command '}'];



math(v_command);
math(s_command);
math(l_command);
math(c_command);
math(n_command);
math(rho_command);
math('h[x_,eta_]:=-v*Min[x,eta]-s*Max[x-eta,0]+l*Max[eta-x,0]+c*x+rho');
math(psy_command);
math('hhat[x_]:=Sum[h[x,etas[[i]]],{i,1,ndata}]/ndata');
math('Cf[a_,b_]:=Sum[(h[a,etas[[i]]]-hhat[a])*(h[b,etas[[i]]]-hhat[b]),{i,1,ndata}] / Sqrt[Sum[(h[a,etas[[i]]]-hhat[a])^2,{i,1,ndata}]*Sum[(h[b,etas[[i]]]-hhat[b])^2,{i,1,ndata}]]');
result=math('NIntegrate[Sqrt[Derivative[1,1][Cf][c,c]],{c,0,50}]');



end

