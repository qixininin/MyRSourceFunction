# MyRSourceFunction


Example1: Simulating genotypes
```
n = 200
m = 1000
freq = runif(m, 0.2, 0.5)
Dprime = runif(m, 0.1, 0.2)
ld = DprimetoD(freq, Dprime)
X = GenerateLDGeno(freq = freq, ld = ld, N = n)
```
