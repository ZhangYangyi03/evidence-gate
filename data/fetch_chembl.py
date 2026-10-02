# 抓 ChEMBL 公共活性数据（EGFR / CHEMBL203，IC50），纯 urllib，无第三方依赖。
# 用法: python fetch_chembl.py [target_id] [standard_type] [max_records]
import json, sys, time, urllib.request

def fetch(target="CHEMBL203", stype="IC50", cap=3000):
    out, off = [], 0
    while off < cap:
        url = ("https://www.ebi.ac.uk/chembl/api/data/activity.json"
               f"?target_chembl_id={target}&standard_type={stype}&limit=1000&offset={off}")
        d = json.loads(urllib.request.urlopen(url, timeout=60).read())
        keep = [a for a in d.get("activities", [])
                if a.get("pchembl_value") and a.get("canonical_smiles") and a.get("standard_relation") == "="]
        out += keep
        print(f"offset={off} kept={len(out)} total={d.get('page_meta',{}).get('total_count')}", file=sys.stderr)
        if not d.get("page_meta", {}).get("next"):
            break
        off += 1000
        time.sleep(0.2)
    return out

if __name__ == "__main__":
    t = sys.argv[1] if len(sys.argv) > 1 else "CHEMBL203"
    s = sys.argv[2] if len(sys.argv) > 2 else "IC50"
    c = int(sys.argv[3]) if len(sys.argv) > 3 else 3000
    recs = fetch(t, s, c)
    json.dump(recs, open(f"raw_{t}_{s}.json", "w", encoding="utf-8"))
    print(f"{len(recs)} records -> raw_{t}_{s}.json")
