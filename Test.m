clear;
clc;

lim = 1000;

for k = 2 : lim
    disp(k);
    o = aks(k);
    assert( isprime(k) == o );
end
