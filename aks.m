function o = aks(n)
% AKS (primality test) : gets a positive integer and determines
% whether it is prime or composite within polynomial time.
% usage : o = aks(n)
%
% arguments: (input)
%  n : integer or vpi bigger than 1
% (variable percision integer (vpi) toolbox available here :
% https://www.mathworks.com/matlabcentral/fileexchange/22725-variable-precision-integer-arithmetic )
%
% arguments: (output)
%  o : boolean value, true if n is prime, false if n is composite
%
% AKS algorithm pseudocode:
%
% Input: integer n > 1.
% 1. If (n = a^b for a in N and b > 1), output COMPOSITE.
% 2. Find the smallest r such that o_r(n) > log2(n).
% 3. If 1 < gcd(a, n) < n for some a <= r, output COMPOSITE.
% 4. If n <= r, output PRIME
% 5. For a = 1 to floor(sqrt(phi(r) * log2 (n))) do
%   if ((X + a)^n ~= X^n + a (mod X^r - 1, n)), output COMPOSITE;
% 6. Output PRIME.
%
% where:
%   o_r(n) : the order of n modulo r
%       Given r in N, n in Z with (n, r) = 1, the order of a modulo r
%       is the smallest number k such that n^k = 1 (mod r).
%   phi(r) : Euler's totient function
%       Number of positive integers up to r that are relatively prime to r
%
% For more information on AKS algorithm:
%   [1]     https://en.wikipedia.org/wiki/AKS_primality_test
%   [2]     Agrawal, Manindra, Neeraj Kayal, and Nitin Saxena. "PRIMES is in P."
%           Annals of mathematics (2004): 781-793.
%   [3]     Aaronson, Scott. "The prime facts: From Euclid to AKS."
%           Manuscript, available from http://www. scottaaronson. com/writings/prime. pdf (2003).
%
% Author: Ali Khodabakhsh
% e-mail: ali.khodabakhsh@gmail.com

% Input: integer n > 1.
% n = vpi(n);
assert(n>1,'Input error');

% 1. If (n = a^b for a in N and b > 1), output COMPOSITE.
ds = log2(n);   % find max(b)
for k = 2 : ds
    if n == floor(nthroot(n,k))^k;
        o = false;
        return
    end
end

% 2. Find the smallest r such that o_r(n) > log2(n).
r = ceil(ds)+1;     % r >= log2(n)+1
while o_r(n,r) <= ds
    r = r + 1;
end

% 3. If 1 < gcd(a, n) < n for some a <= r, output COMPOSITE.
for k = 2 : r
    if 1 < gcd(n,k) && gcd(n,k) < n
        o = false;
        return
    end
end

% 4. If n <= r, output PRIME.
if n <= r
    o = true;
    return
end

% 5. For a = 1 to floor(sqrt(phi(r) * log2(n))) do
%   if ((X + a)^n ~= X^n + a (mod X^r - 1, n)), output COMPOSITE;

% polynomial in vector form :
%   for polynomial = a_0 + a_1*x + a_2*x^2 + ... + a_k*x^k
%   vector = [a_0, a_1, a_2, ... , a_k]

modx = -1; modx(r) = 1;     % modx is polynomial vector representing x^r - 1
s2f = pmx([0 1],n,modx,n);  % s2f is polynomial vector representing x^n mod (x^r - 1, n)
temp = zeros(1,r);          % temp is an empty polynomial vector

maxa = floor(sqrt(phi(r)*ds));  % maximum a to check
for a = 1 : maxa
    temp(1) = a;
    s1 = pmx([a 1],n,modx,n);   % s1 is polynomial vector representing (x + a)^n mod (x^r - 1, n)
    s2 = s2f + temp;            % s2 is polynomial vector representing x^n + a mod (x^r - 1, n)
    if ~all(s1 == s2)
        o = false;
        return
    end
end

% 6. Output PRIME.
o = true;
return

end

function k = o_r(a, r)
% o_r : gets two numbers a and r, and determines the order of a modulo r
% usage : k = o_r(a, r)
%
% arguments: (input)
%  a, r : integer or vpi
%
% arguments: (output)
%  k : smallest number for which a^k = 1 (mod r)
%
% o_r(a) : the order of a modulo r
%   Given r in N, a in Z with (a, r) = 1, the order of a modulo r is the smallest number k such that a^k = 1 (mod r).

a = mod(a,r);
if gcd(a,r) ~= 1,   % if gcd(a, r) ~= 1, there is no k satisfying a^k = 1 (mod r)
    k = nan;
    return
end

k = 1;
while powermod(a,k,r) ~= 1  % while a^k ~= 1 (mod r), increment k
    k = k + 1;
end

end

function t = phi(r)
% phi : gets a number r and counts the positive integers up to r that are relatively prime to r. (Euler's totient function)
% usage : t = phi(r)
%
% arguments: (input)
%  r : integer or vpi
%
% arguments: (output)
%  t : integer representing the number of positive integers up to r that are relatively prime to r
%
% phi(r) : Euler's totient function
%   Number of positive integers up to r that are relatively prime to r
%
% For more information on Euler's totient function:
%   [1]     https://en.wikipedia.org/wiki/Euler%27s_totient_function

f = unique(factor(r));  % f = unique factors of r
t = r*prod(1-1./f);

end

function v = pmx(x,k,r,n)
% pmx : gets two vector representation of polynomials x and r, and two
% integers k and n, and returns the vector v, representing remainder x^k mod (r, n)
% usage : v = pmx(x,k,r,n)
%
% arguments: (input)
%  x, r : vector representation of two polynomials
%  k, r : positive integers
%
% arguments: (output)
%  v : remainder of x^k mod(r, n)

assert(k>0);    % make sure k is positive

if k == 1   % base case
    v = x;
elseif mod(k,2) == 0    % if k is even, v = x^(k/2) * x^(k/2) (mod r, n);
    v = pmx(x,k/2,r,n);     % calculate x^(k/2) (mod r, n)
    v = multx(v,v);
else    % if k is odd, v = x^(k-1) * x (mod r, n);
    v = pmx(x,k - 1,r,n);   % calculate x^(k-1) (mod r, n)
    v = multx(v,x);
end

v = mx(v,r,n);  % compute v mod (r, n)

v = fixlen(v,r);    % make sure the length of the output vectors is fixed

end

function v = multx(x,y)
% multx : gets two vector representation of polynomials x and y, and
% computes vector representation of x * y
% usage : v = multx(x,y)
%
% arguments: (input)
%  x, y : vector representation of two polynomials
%
% arguments: (output)
%  v : vector representation of x * y

mat = x'*y;     % multiply vectors and get a matrix
vlen = length(x) + length(y) - 1;   % length of output vector
v = zeros(1, vlen);     % initialize output vector
for k = 1 : vlen    % for each element in output vector, k
    for l = 1 : k   % add all elements in matrix mat corresponding to x^k
        if k-l+1<=length(x) && l<=length(y)     % of course if elements exist
            v(k) = v(k) + mat(k-l+1,l);
        end
    end
end
v = shorten(v);     % trim the zeros of the tail of the vector if there is any

end

function v = mx(x,r,n)
% mx : gets two vector representation of polynomials x and r, and an
% integer n, and computes x mod (r, n)
% usage : v = mx(x,r,n)
%
% arguments: (input)
%  x, r : vector representation of two polynomials
%  n : integer
%
% arguments: (output)
%  v : vector representation of x mod (r, n)

len = length(r);
v = mod(x, n);  % initialize v by modulo n of all elements in x
v = shorten(v);     % remove trailing zeros
while length(v) >= len  % while v is divisible by r
    m = minv(r(end), n)*v(end);     % find an integer m, such that m * v(end)/r(end) = 1 (mod n)
    v(end-len+1:end) = v(end-len+1:end) - m*r;  % subtract v by r * m`
    v = mod(v,n);   % find modulo n of all elements in v
    v = shorten(v);     % remove trailing zeros
end

end

function r = powermod(x,k,n)
% powermod : gets three numbers x, k, and n, and computes x^k (mod n)
% efficiently
% usage : r = powermod(x,k,n)
%
% arguments: (input)
%  x, k, n : integer
%
% arguments: (output)
%  r : x^k (mod n)

assert(k>0);    % make sure k is positive
x = mod(x,n);   % put x in range of 1 and n

if k == 1   % base case
    r = x;
elseif mod(k,2) == 0    % if k is even, r = x^(k/2) * x^(k/2) (mod n);
    r = powermod(x,k/2,n);  % calculate x^(k/2) (mod r, n)
    r = mod(r^2,n);
else    % if k is odd, v = x^(k-1) * x (mod n);
    r = powermod(x,k - 1,n);    % calculate x^(k-1) (mod r, n)
    r = mod(r*x,n);
end

end

function t = minv(a, n)
% minv : gets two integers a and n, and returns multiplicative modular
% inverse of a (mod n), that is t, when a * t = 1 (mod n), using Extended Euclidean algorithm
% usage : t = minv(a, n)
%
% arguments: (input)
%  a, n : two integers
%
% arguments: (output)
%  t : integer, such that a * t = 1 (mod n)
% 
% Ref:
%   [1]     https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Computing_multiplicative_inverses_in_modular_structures

% initialize t, newt, r, and newr
t = 0;
nt = 1;
r = n;
nr = a;

while nr ~= 0   % while newr is not zero
    q = floor(r/nr);    % calculate quotient
    [t, nt] = update(t, nt, q);     % update t and nt using q
    [r, nr] = update(r, nr, q);     % update r and nr using q
end

if r>1  % if r>1, a is not invertible
    error('a is not invertible');
elseif t<0  % fix negative t
    t = t + n;
end

end

function [x, nx] = update(x, nx, q)
% update : updates x and nx, based on q for the extend euclidean algorithm
% for calculating multiplicative modular inverse
% usage : [x, nx] = update(x, nx, q)
%
% arguments: (input)
%  x, nx, q : integers
%
% arguments: (output)
%  x, nx : updated x and nx such that x = nx, and nx = x - q * nx
% 
% Ref:
%   [1]     https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Computing_multiplicative_inverses_in_modular_structures

ux = nx;
unx = x - q * nx;

x = ux;
nx = unx;

end

function v = shorten(v)
% shorten : gets a vector v and trims zeros from the end of it
% usage : v = shorten(v)
%
% arguments: (input)
%  v : a vector
%
% arguments: (output)
%  v : vector v with no trailing zeros

ind = find(v~=0,1,'last');  % find last non zero element
v = v(1:ind);

end

function v = fixlen(v,r)
% fixlen : gets two vectors v and r, and pads zeros to the end of v so that its length matches r
% usage : v = fixlen(v,r)
%
% arguments: (input)
%  v, r : two vectors
%
% arguments: (output)
%  v : vector v with a length equal to r

assert(length(v)<=length(r));   % Assert that v is shorter than r
if length(v)<length(r)
    v(length(r)) = 0;
end

end
