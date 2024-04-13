clear all
clc
N1 = 21;
N2 = 21;
V1 = 4;
V2 = 4;
power = 0;
for i = 1:V1/2
    a = (1/N1^2)*(1/(sin(((2*i-1)*pi)/(2*N2))^2));
    for j = 1:V2/2
        b = (1/N2^2)*(1/(sin(((2*j-1)*pi)/(2*N2))^2));
        power = power + a*b;
    end
end
power = power * 4;