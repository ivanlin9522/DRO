phat_one = zeros(7,1);
CIup_one = zeros(7,1);
CIdown_one = zeros(7,1);
phat_two = zeros(7,1);
CIup_two = zeros(7,1);
CIdown_two = zeros(7,1);
j=1;
for i =20:10:80
    [phat_one(j), CIup_one(j), CIdown_one(j), phat_two(j), CIup_two(j), CIdown_two(j)] = multifunction(0,1,i);
    j = j + 1;
end