import pymupdf, os

d = r"c:\Users\28985\.trae-cn\attachments\6a8d4ed36c15549fcd338e7d"
files = [f for f in os.listdir(d) if f.lower().endswith(".pdf")]

out = []
for f in files:
    path = os.path.join(d, f)
    out.append("="*80)
    out.append("FILE: " + f)
    out.append("="*80)
    doc = pymupdf.open(path)
    full = ""
    for page in doc:
        full += page.get_text()
    out.append("--- HEAD ---")
    out.append(full[:600])
    idx = full.lower().find("keywords")
    if idx != -1:
        out.append("--- KEYWORDS ---")
        out.append(full[idx:idx+400])
    ai = full.lower().find("abstract")
    if ai != -1:
        out.append("--- ABSTRACT ---")
        out.append(full[ai:ai+1800])
    doc.close()

with open(r"e:\尼古丁UC科研材料\555\_pdfs_extract.txt","w",encoding="utf-8") as fh:
    fh.write("\n".join(out))
print("total chars:", sum(len(x) for x in out))
print("saved")
