
import sys

with open(sys.argv[1], "r") as f:
    s = f.read()
    print("wtm", sum(int(l.split()[-1]) for l in s.splitlines() if l.startswith("wtm")))
    print("btm", sum(int(l.split()[-1]) for l in s.splitlines() if l.startswith("btm")))