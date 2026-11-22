"""Probe RCSB REST for candidate entries: title, method, resolution, organism, entity list."""
import sys, json, urllib.request

def get(url):
    with urllib.request.urlopen(url, timeout=30) as f:
        return json.load(f)

for pdb in sys.argv[1:]:
    p = pdb.upper()
    try:
        e = get(f"https://data.rcsb.org/rest/v1/core/entry/{p}")
    except Exception as ex:
        print(f"{p}: FETCH FAIL {ex}"); continue
    title = e.get("struct", {}).get("title", "?")
    method = [m.get("method") for m in e.get("exptl", [])]
    res = e.get("rcsb_entry_info", {}).get("resolution_combined")
    pcit = e.get("rcsb_primary_citation", {}) or {}
    cit = f"{', '.join(pcit.get('rcsb_authors',[])[:2])} et al. {pcit.get('rcsb_journal_abbrev','')} {pcit.get('year','')}; doi:{pcit.get('pdbx_database_id_doi','')}"
    npoly = e.get("rcsb_entry_info", {}).get("polymer_entity_count")
    print(f"\n=== {p}  {method} {res} A")
    print(f"  TITLE: {title}")
    print(f"  CIT:   {cit}")
    print(f"  entities: {npoly}")
    ids = e.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids", [])
    for eid in ids:
        try:
            pe = get(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{p}/{eid}")
        except Exception:
            continue
        d = pe.get("rcsb_polymer_entity", {}).get("pdbx_description", "?")
        typ = pe.get("entity_poly", {}).get("rcsb_entity_polymer_type", "?")
        ln = pe.get("entity_poly", {}).get("rcsb_sample_sequence_length")
        src = [s.get("scientific_name") for s in (pe.get("rcsb_entity_source_organism") or [])]
        chains = pe.get("rcsb_polymer_entity_container_identifiers", {}).get("auth_asym_ids", [])
        if typ != "Protein" or len(ids) < 12:
            print(f"   ent{eid:>3} {typ:<10} len={ln:<6} chains={','.join(chains):<12} {src} :: {d}")
    nprot = sum(1 for _ in ids)
