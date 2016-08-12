% Clear workspace
clear;
clc;

% Maximum number to check
lim = 1000;

% Check from 2 to lim, whether it gives same result as matlab'd defaul isprime
for k = 2 : lim
    fprintf('%d\t',k);
    tic;
    o = aks(k);
    fprintf('./   ');
    toc
    assert( isprime(k) == o );
end
