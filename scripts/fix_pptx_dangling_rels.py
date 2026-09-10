"""
scripts/fix_pptx_dangling_rels.py

Make a pandoc-produced .pptx openable in PowerPoint.

TWO BUGS
--------
(1) EMPTY SHAPES. With a `date:` field in the YAML, pandoc emits a bare
`<p:sp />` into the slide's shape tree: a shape with no non-visual properties
and no shape properties, which is invalid PresentationML. PowerPoint refuses
the whole file with "PowerPoint can't read..." (HRESULT 0x80070570, which reads
like a disk error and is not). LibreOffice ignores it, so `soffice --convert-to
pdf` succeeds and the defect goes unnoticed. Removing `date:` from the YAML
avoids it at source; this script also strips the empty shapes.

(2) DANGLING RELATIONSHIPS
--------------------------
Quarto/pandoc copies `ppt/_rels/presentation.xml.rels` from the reference
document but does not copy every part that file points at. Our reference doc
(MN-proxy-Ghana-presentation.pptx) carries three customXml parts and an
authors.xml (comment authors); pandoc drops the parts and keeps the
relationships. The result is four relationships whose Target does not exist in
the package.

PowerPoint validates that every relationship target resolves and refuses the
file with "PowerPoint can't read...". LibreOffice ignores dangling targets,
which is why `soffice --convert-to pdf` succeeded and the defect went unnoticed.

The dangling relationships are never referenced by r:id from the body of the
owning part, so removing them is safe and changes nothing about the slides.

WHAT THIS DOES
--------------
Walks every *.rels part, drops each internal (non-external) Relationship whose
Target is missing from the package AND whose Id is not referenced from the body
of the part that owns the .rels. Anything still referenced is left alone and
reported as an error, because that would be a different and more serious defect.

    python scripts/fix_pptx_dangling_rels.py <file.pptx> [more.pptx ...]

Edits in place, after writing <file>.bak once.
"""
import posixpath
import re
import shutil
import sys
import zipfile

REL_RE = re.compile(r"<Relationship\b[^>]*/>")
ATTR = lambda name, s: (re.search(r'%s="([^"]*)"' % name, s) or [None, None])[1]


def owning_part(rels_name):
    """ppt/_rels/presentation.xml.rels -> ppt/presentation.xml"""
    d = posixpath.dirname(posixpath.dirname(rels_name))
    base = posixpath.basename(rels_name)[: -len(".rels")]
    return posixpath.join(d, base) if base else None


def fix(path):
    with zipfile.ZipFile(path) as z:
        names = set(z.namelist())
        items = {n: z.read(n) for n in z.namelist()}
        infos = {i.filename: i for i in z.infolist()}

    removed, kept_but_broken, empty_shapes = [], [], 0

    # (1) strip bare <p:sp/> elements from slides, layouts and masters
    for n in [x for x in items if re.match(r"ppt/(slides|slideLayouts|slideMasters|notesSlides)/[^/]+\.xml$", x)]:
        text = items[n].decode("utf8")
        new = re.sub(r"<p:sp\s*/>", "", text)
        if new != text:
            empty_shapes += len(re.findall(r"<p:sp\s*/>", text))
            items[n] = new.encode("utf8")

    for rels in [n for n in list(items) if n.endswith(".rels")]:
        body_name = owning_part(rels)
        body = items.get(body_name, b"").decode("utf8", "replace") if body_name else ""
        used = set(re.findall(r'r:(?:id|embed|link)="(rId\d+)"', body))

        text = items[rels].decode("utf8")
        out, changed = text, False
        for tag in REL_RE.findall(text):
            if ATTR("TargetMode", tag) == "External":
                continue
            target = ATTR("Target", tag) or ""
            if target.startswith(("http://", "https://")):
                continue
            # A .rels at X/_rels/Y.rels resolves its targets against X.
            # For the package root (_rels/.rels) that base is "".
            base = posixpath.dirname(posixpath.dirname(rels))
            full = posixpath.normpath(posixpath.join(base, target)) if base else posixpath.normpath(target)
            if full in names:
                continue
            rid = ATTR("Id", tag)
            if rid in used:
                kept_but_broken.append((rels, rid, target))
                continue
            out = out.replace(tag, "")
            changed = True
            removed.append((rels, rid, target))
        if changed:
            items[rels] = out.encode("utf8")

    if not removed and not empty_shapes:
        print("%s: nothing to fix" % path)
    else:
        shutil.copyfile(path, path + ".bak")
        with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
            for n, data in items.items():
                zi = infos.get(n)
                z.writestr(zi if zi is not None else n, data)
        print("%s: removed %d empty shape(s), %d dangling relationship(s)"
              % (path, empty_shapes, len(removed)))
        for r, rid, t in removed:
            print("    %s  %s -> %s" % (r, rid, t))

    for r, rid, t in kept_but_broken:
        print("  !! %s %s points at missing %s AND is referenced in the body; "
              "not safe to strip, investigate" % (r, rid, t))
    return len(kept_but_broken) == 0


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    ok = all(fix(p) for p in sys.argv[1:])
    sys.exit(0 if ok else 1)
