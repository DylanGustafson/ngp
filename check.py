from sympy.ntheory.residue_ntheory import is_primitive_root
from sympy.ntheory import isprime, primefactors
import sys

def build_line(N):
    if not isprime(N):
        print(f"{N} is not prime: {primefactors(N)}")
        exit()

    for root in range(2, N):
        if is_primitive_root(root, N):
            break

    q = (pow(root, N, N * N) - root) // N
    pfacs = primefactors(N - 1)

    return f"{N} ({root} {q}) " + " ".join([str(i) for i in pfacs])

arg = sys.argv[-1]
try:
    N = int(arg)
    print(build_line(N))
    exit()

except ValueError:
    pass

with open(arg, 'r') as f:
    lines = f.read().splitlines()

line_count = len(lines)
for i, line in enumerate(lines):
    print(f"{i + 1}/{line_count}", end="\r")
    Nstr = line.split(" ", 1)[0]
    if not Nstr.isdecimal():
        continue

    true_line = build_line(int(Nstr))
    if true_line != line:
        print("Value Discrepency!")
        print("Given: " + line)
        print("True : " + true_line)
        exit()

print("\nPassed!")
