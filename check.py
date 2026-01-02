from sympy.ntheory.residue_ntheory import is_primitive_root
from sympy.ntheory import isprime, primefactors
import sys

N = int(sys.argv[-1])
if not isprime(N):
    print(f"Input is not prime: {primefactors(N)}")
    exit()

for root in range(2, N):
    if is_primitive_root(root, N):
        break

print(f"Minimum primitive root: {root}")
residue = (pow(root, N - 1, N * N) - 1) // N
print(f"(r^(p-1) % p^2 - 1) /p: {residue}")
