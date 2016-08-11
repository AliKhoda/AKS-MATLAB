function o = aks(n)

% Input: integer n > 1.
% n = vpi(n);
assert(n>1,'Input error');

% 1. If (n = a^b for a ? N and b > 1), output COMPOSITE.
ds = log2(n);
for k = 2 : ds
    if n == floor(nthroot(n,k))^k;
        o = false;
        return
    end
end

% 2. Find the smallest r such that o_r(n) > log2(n).
r = ceil(ds)+1;
while o_r(n,r) < ds
    r = r + 1;
end

% 3. If 1 < gcd(a, n) < n for some a ? r, output COMPOSITE.
for k = 2 : r
    if gcd(n,k) > 1 && gcd(n,k) < n
        o = false;
        return
    end
end

% 4. If n <= r, output PRIME.
if n<=r
    o = true;
    return
end

% 5. For a = 1 to floor(sqrt(?(r) log n)) do
% if ((X + a)^n ~= X^n + a (mod X^r - 1, n)), output COMPOSITE;
modx = -1;
modx(r) = 1;
temp = zeros(1,r);
maxa = floor(sqrt(totient(r)*ds));
s2f = pmx([0 1],n,modx,n);
for a = 1 : maxa
    temp(1) = a;
    s1 = pmx([a 1],n,modx,n);
    s2 = s2f + temp;
    if ~all(s1 == s2)
        o = false;
        return
    end
end

% 6. Output PRIME.
o = true;
return

end

function k = o_r(a,r)

a = mod(a,r);
if gcd(a,r) ~= 1,
    k = nan;
    return
end

k = 1;
while powermod(a,k,r) ~= 1
    k = k + 1;
end
    
end

function r = isdiv(n,p)

r = mod(n,p) == 0;

end

function t = totient(n)

f=unique(factor(n));
t=n*prod(1-1./f);

end

function v = pmx(x,k,r,n)

assert(k>0);

if k == 1
    v = x;
elseif mod(k,2) == 0
    v = pmx(x,k/2,r,n);
    v = multx(v,v);
else
    v = pmx(x,k - 1,r,n);
    v = multx(v,x);
end

v = mx(v,r,n);

v = fixlen(v,r);

end

function v = multx(x,y)

mat = x'*y;
v = zeros(1, 2*length(x));
for k = 1 : length(x)*2
    for l = 1 : k
        if k-l+1<=length(x) && l<=length(y)
            v(k) = v(k) + mat(k-l+1,l);
        end
    end
end
v = shorten(v);

end

function v = mx(x,r,n)

len = length(r);
v = x;
while length(v) >= len
    assert(mod(x(end),r(end)) == 0);
    m = v(end)/r(end);
    v(end-len+1:end) = v(end-len+1:end) - m*r;
    v = shorten(v);
end

v = mod(v,n);
    
end

function v = shorten(v)

ind = find(v~=0,1,'last');
v = v(1:ind);

end

function v = fixlen(v,r)

assert(length(v)<=length(r));
if length(v)<length(r)
    v(length(r)) = 0;
end

end

function r = powermod(x,k,n)

assert(k>0);
x = mod(x,n);

if k == 1
    r = x;
elseif mod(k,2) == 0
    r = powermod(x,k/2,n);
    r = mod(r^2,n);
else
    r = powermod(x,k - 1,n);
    r = mod(r*x,n);
end

end
