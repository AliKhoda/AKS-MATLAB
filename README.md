# AKS-MATLAB

An implementation of the original AKS algorithm in MATLAB as described in the Annals of mathematics

The AKS primality test (also known as Agrawal–Kayal–Saxena primality test and cyclotomic AKS test) is a deterministic primality-proving algorithm created and published by Manindra Agrawal, Neeraj Kayal, and Nitin Saxena, computer scientists at the Indian Institute of Technology Kanpur, on August 6, 2002, in a paper titled "PRIMES is in P".
The algorithm determines whether a number is prime or composite within polynomial time. The asymptotic time complexity of the algorithm is O~((log2(n))^15/2).

AKS algorithm pseudocode:

Input: integer n > 1.  
1. If (n = a^b for a in N and b > 1), output COMPOSITE.  
2. Find the smallest r such that o_r(n) > log2(n).  
3. If 1 < gcd(a, n) < n for some a <= r, output COMPOSITE.  
4. If n <= r, output PRIME.  
5. For a = 1 to floor(sqrt(phi(r) * log2 (n))) do  
      if ((X + a)^n ~= X^n + a (mod X^r - 1, n)), output COMPOSITE;  
6. Output PRIME.  

where:  
- o_r(n) (the order of n modulo r) : Given r in N, n in Z with (n, r) = 1, the order of a modulo r is the smallest number k such that n^k = 1 (mod r).  
- phi(r) (Euler's totient function) : Number of positive integers up to r that are relatively prime to r

References:  
[1] https://en.wikipedia.org/wiki/AKS_primality_test  
[2] Agrawal, Manindra, Neeraj Kayal, and Nitin Saxena. "PRIMES is in P." Annals of mathematics (2004): 781-793.
