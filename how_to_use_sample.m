%1 sided, n=20, DOF=1, discrete=1
phat = zeros(5,1);
CIup = zeros(5,1);
CIdown = zeros(5,1);
j=1;
for i =40:10:80
    [phat(j), CIup(j), CIdown(j)] = multifunction(0,0,1,i)
    j = j + 1;
end