import os
import pathlib

os.makedirs("6men_minimal", exist_ok=True)
os.makedirs("6men_minimal/0pawns", exist_ok=True)
os.makedirs("6men_minimal/1pawns", exist_ok=True)
os.makedirs("6men_minimal/2pawns", exist_ok=True)
os.makedirs("6men_minimal/3pawns", exist_ok=True)
os.makedirs("6men_minimal/4pawns", exist_ok=True)
os.makedirs("6men_remaining", exist_ok=True)
os.makedirs("6men_remaining/0pawns", exist_ok=True)
os.makedirs("6men_remaining/1pawns", exist_ok=True)
os.makedirs("6men_remaining/2pawns", exist_ok=True)
os.makedirs("6men_remaining/3pawns", exist_ok=True)
os.makedirs("6men_remaining/4pawns", exist_ok=True)

tbs = {}
for d, _, files in list(os.walk("6men")) + list(os.walk("6men_minimal")) + list(os.walk("6men_remaining")):
    for file in files:
        if file.endswith(".bz"):
            
            tbs[file[:-len(".egtb.bz")]] = (pathlib.Path(d, file), pathlib.Path(d, file).stat().st_size)
            
            
            
size = sum(size for _, size in tbs.values())
print(size / 1024**4)

minimal_size = 0
for egtb_id, (path, size) in tbs.items():
    _, a, b = egtb_id.split("K")
    mirror_egtb_id = f"K{b}K{a}"
    if egtb_id != mirror_egtb_id:
        if egtb_id < mirror_egtb_id:
            mirror_path, mirror_size = tbs[mirror_egtb_id]
            minimal_size += min(mirror_size, size)
            
            old_path = path if size < mirror_size else mirror_path
            if "minimal" not in str(old_path):
                new_path = pathlib.Path("6men_minimal", old_path.relative_to("6men"))
                os.rename(old_path, new_path)
                print(old_path, "=>", new_path)

            old_path = mirror_path if size < mirror_size else path
            if "remaining" not in str(old_path):
                new_path = pathlib.Path("6men_remaining", old_path.relative_to("6men"))
                os.rename(old_path, new_path)
                print(old_path, "=>", new_path)

            
    else:
        minimal_size += size
        old_path = path
        if "minimal" not in str(old_path):
            new_path = pathlib.Path("6men_minimal", old_path.relative_to("6men"))
            os.rename(old_path, new_path)
            print(old_path, "=>", new_path)
print(minimal_size / 1024**4)