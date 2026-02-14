from sympy.ntheory.residue_ntheory import is_primitive_root
from sympy.ntheory import isprime, primefactors
from subprocess import run
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

if len(sys.argv) == 1:
    fname = sys.argv[0]
    print("Usage:")
    print(f"  Check a single number        : python3 {fname} number")
    print(f"  Check a stderr output file   : python3 {fname} file_name")
    print(f"  Run ngp-bin and check output : python3 {fname} start_value total_size")
    print("\nLine format: number (root q_value) pfac1 pfac2 ...")
    exit()

elif len(sys.argv) > 2:
    cmd = ['./ngp-bin', sys.argv[1], sys.argv[2], sys.argv[2]]
    print("Running:", " ".join(cmd))
    output = run(cmd, capture_output=True, text=True).stderr
    print("Done!", end=" ")

elif sys.argv[-1].isdecimal():
    N = int(sys.argv[-1])
    print("True output:", build_line(N))
    exit()

else:
    with open(sys.argv[-1], 'r') as f:
        output = f.read()

print("Checking...")
lines = output.splitlines()
line_count = len(lines)
for i, out_line in enumerate(lines):
    print(f"{i + 1}/{line_count}", end="\r")

    Nstr = out_line.split(" ", 1)[0]
    if not Nstr.isdecimal():
        continue

    true_line = build_line(int(Nstr))
    if true_line != out_line:
        print("Value Discrepency!")
        print(" Output: " + out_line)
        print(" True  : " + true_line)
        exit()

print("\nPassed!")
