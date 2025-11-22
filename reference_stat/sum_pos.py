
import sys
filename = sys.argv[1]
wtm = filename.split("/")[-1].split(".tbs")[0]
btm = "k" + "k".join(wtm.split("k")[2:0:-1])
with open(filename, "r") as f:
    s = f.read()
    wtm_sum = sum(int(l.split()[-1]) for l in s.splitlines() if l.startswith("wtm"))
    wtm_broken = sum(int(l.split()[-1]) for l in s.splitlines() if l.startswith("wtm: Broken"))
    btm_sum = sum(int(l.split()[-1]) for l in s.splitlines() if l.startswith("btm"))
    btm_broken = sum(int(l.split()[-1]) for l in s.splitlines() if l.startswith("btm: Broken"))
    print("wtm", wtm, wtm_sum, "broken", wtm_broken)
    print("btm", btm, btm_sum, "broken", btm_broken)