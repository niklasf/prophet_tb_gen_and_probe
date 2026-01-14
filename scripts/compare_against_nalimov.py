import requests
from time import sleep
with open("longest_mates.csv", "r") as f:
    count = 0
    correct = 0
    for l in f.readlines()[1:]:
        egtb, _, _, fen, dtm_str, _ = l.split(",")
        print(egtb, fen, dtm_str)
        if dtm_str == "":
            continue
        dtm = int(dtm_str)
        print(egtb, fen, "dtm:", dtm)
        if (len(egtb) <= 5): continue
        
        res = requests.get("https://www.k4it.de/egtb/fetch.php", params={
            "action": "egtb",
            "fen": fen
        })
        if (res.status_code != 200):
            print("Error:", res.status_code)
            continue
        mate_in_text = res.text.splitlines()[0].partition("|")[2]
        print("\tvs Nalimov:", mate_in_text, end=" ")
        if "Not found" in mate_in_text:
            print("... no tablebase info")
        else:
            mate_in = int(mate_in_text.split()[-1])
            is_equal = mate_in == dtm // 2 + 1
            print("... matches:", is_equal)
            count += 1
            correct += is_equal
        
        sleep(1)
    print(f"Total tested: {count}, Correct: {correct}/{count}")