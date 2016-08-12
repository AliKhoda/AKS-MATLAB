% Clear workspace
clear;
clc;

% Maximum number to check
lim = 1000;

% Check from 2 to lim, whether aks gives same result as MATLAB's default isprime function
for k = 2 : lim
    fprintf('%d\t',k);
    tic;
    o = aks(k);
    fprintf('./   ');
    toc
    assert( isprime(k) == o );
end
