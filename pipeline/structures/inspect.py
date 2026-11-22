"""Inspect the downloaded mmCIF: chains, entity names, resolved residue counts."""
import sys, gemmi, numpy as np

path = sys.argv[1]
st = gemmi.read_structure(path)
st.setup_entities()
doc = gemmi.cif.read(path)
blk = doc.sole_block()
# entity id -> description
ent_desc = {r[0]: r[1] for r in blk.find('_entity.', ['id', 'pdbx_description'])}
# auth chain -> entity id
ch2ent = {}
for r in blk.find('_struct_asym.', ['id', 'entity_id']):
    pass
for r in blk.find('_atom_site.', ['auth_asym_id', 'label_entity_id']):
    ch2ent.setdefault(r[0], r[1])

model = st[0]
print(f"model chains: {len(model)}")
rows = []
for ch in model:
    pol = ch.get_polymer()
    kind = str(pol.check_polymer_type())
    nres = len([r for r in ch if r.name not in ('HOH',)])
    eid = ch2ent.get(ch.name, '?')
    desc = ent_desc.get(eid, '?').strip("'\"")
    first = ch[0].seqid.num if len(ch) else None
    last = ch[len(ch)-1].seqid.num if len(ch) else None
    rows.append((ch.name, eid, kind, nres, first, last, desc))
for r in sorted(rows, key=lambda x: -x[3]):
    print(f"{r[0]:>3}  ent={r[1]:>3}  {r[2]:<28} n={r[3]:<5} {r[4]}..{r[5]}  {r[6]}")
