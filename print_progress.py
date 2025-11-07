import os
from datetime import datetime
files = [f[:-9] for f in os.listdir("egtbs") if f.endswith(".zip")]
for i in range(2,7):
    print(f"{i}: {sum(1 for f in files if len(f) == i):3d}")
print("  ", len(files), "   ", datetime.now())