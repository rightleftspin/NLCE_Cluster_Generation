import json
import glob
import os

old_file_list = sorted(glob.glob(os.path.join(os.getcwd(), "old", "*.json")))
new_file_list = sorted(glob.glob(os.path.join(os.getcwd(), "new", "*.json")))

for old, new in zip(old_file_list, new_file_list):
   # print(old[78::])
    #print(new[78::])
    old_json = json.load(open(old, mode="r", encoding="utf-8"))
    new_json = json.load(open(new, mode="r", encoding="utf-8"))
    if old_json == new_json:
        print("true")
    else:
        for key1, key2 in zip(old_json, new_json):
            old_json[key1] = [set(tuple(x) for x in old_json[key1][0]), old_json[key1][1]]
            new_json[key2] = [set(tuple(x) for x in new_json[key2][0]), new_json[key2][1]]
        if old_json == new_json:
            print("now_true")




