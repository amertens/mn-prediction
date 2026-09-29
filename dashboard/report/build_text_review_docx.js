// =============================================================================
// dashboard/report/build_text_review_docx.js
//
// Turns a capture folder (capture_dashboard_text.R + slice_screens.py) into a
// Word document for reviewing the dashboard's writing: per page, the
// screenshot and then every piece of visible text in reading order, so edits
// can be suggested with Word's comments or tracked changes.
//
//   node build_text_review_docx.js <capture_dir> <out.docx> <front.json>
//
// front.json: { "title": ..., "subtitle": ..., "intro": [paragraphs],
//               "changes_title": ..., "changes": [bullets] }
// =============================================================================
const fs = require("fs");
const path = require("path");
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, ImageRun, Table, TableRow, TableCell,
  WidthType, ShadingType, BorderStyle, AlignmentType, LevelFormat, Footer, PageNumber, PageBreak,
} = require("docx");

const [capDir, outPath, frontPath] = process.argv.slice(2);
const man = JSON.parse(fs.readFileSync(path.join(capDir, "manifest.json"), "utf8"));
const slices = JSON.parse(fs.readFileSync(path.join(capDir, "slices.json"), "utf8"));
const front = JSON.parse(fs.readFileSync(frontPath, "utf8"));

const CONTENT_W = 9360;          // 6.5 inches in DXA
const IMG_W = 624;               // 6.5 inches in pixels at 96 dpi
const words = (s) => (s || "").split(/\s+/).filter((w) => /[A-Za-z0-9]/.test(w)).length;

function blockWords(b) {
  if (b.t === "table") return [b.head || [], ...(b.rows || [])].flat().reduce((a, c) => a + words(c), 0);
  if (b.t === "ctrl") return words(b.s) + (b.opts || []).reduce((a, c) => a + words(c), 0);
  if (["h", "p", "li", "note", "dt", "btn"].includes(b.t)) return words(b.s);
  return 0; // charts, maps and interactive data tables are not counted
}
const pageWords = (pg) => (pg.blocks || []).reduce((a, b) => a + blockWords(b), 0)
  + ((pg.extra && pg.extra.captions) || []).reduce((a, c) => a + words(c.caption), 0);

const run = (text, opts = {}) => new TextRun({ text, ...opts });
const para = (text, opts = {}) => new Paragraph({ children: [run(text, opts.run || {})], spacing: { after: 100 }, ...(opts.para || {}) });

function table(head, rows, small = true) {
  const n = Math.max(head.length, ...rows.map((r) => r.length), 1);
  const pad = (r) => [...r, ...Array(n - r.length).fill("")].slice(0, n);
  const base = Math.floor(CONTENT_W / n);
  const widths = Array(n).fill(base); widths[n - 1] += CONTENT_W - base * n;
  const mkRow = (cells, header) => new TableRow({
    tableHeader: header,
    children: pad(cells).map((c, i) => new TableCell({
      width: { size: widths[i], type: WidthType.DXA },
      shading: header ? { fill: "E8EEF1", type: ShadingType.CLEAR, color: "auto" } : undefined,
      margins: { top: 40, bottom: 40, left: 80, right: 80 },
      children: [new Paragraph({ children: [run(String(c || ""), { size: small ? 16 : 18, bold: header })] })],
    })),
  });
  const all = [];
  if (head && head.length) all.push(mkRow(head, true));
  rows.forEach((r) => all.push(mkRow(r, false)));
  return new Table({ width: { size: CONTENT_W, type: WidthType.DXA }, columnWidths: widths, rows: all });
}

function blockToDocx(b) {
  const out = [];
  switch (b.t) {
    case "h": out.push(new Paragraph({ heading: HeadingLevel.HEADING_3, children: [run(b.s)] })); break;
    case "p": out.push(para(b.s)); break;
    case "dt": out.push(para(b.s, { run: { bold: true } })); break;
    case "li": out.push(new Paragraph({ numbering: { reference: "bullets", level: 0 }, children: [run(b.s)], spacing: { after: 60 } })); break;
    case "note": out.push(new Paragraph({
      children: [run(b.s, { size: 20 })],
      shading: { fill: "F3F4F6", type: ShadingType.CLEAR, color: "auto" },
      border: { left: { style: BorderStyle.SINGLE, size: 18, color: "2C7BB6", space: 8 } },
      indent: { left: 160 }, spacing: { before: 80, after: 140 },
    })); break;
    case "btn": out.push(para(`[Button] ${b.s}`, { run: { italics: true, color: "555555", size: 18 } })); break;
    case "ctrl": {
      const bits = [`[Control] ${b.s || "(no label)"}`];
      if (b.opts && b.opts.length) bits.push(`Options: ${b.opts.join(" | ")}`);
      if (b.sel) bits.push(`Showing: ${b.sel}`);
      out.push(para(bits.join(". "), { run: { italics: true, color: "555555", size: 18 } }));
      break;
    }
    case "widget": out.push(para(b.s, { run: { italics: true, color: "888888", size: 18 } })); break;
    case "table": out.push(table(b.head || [], b.rows || [])); out.push(para("")); break;
    case "rtable":
      out.push(para("[Interactive table, first rows shown]", { run: { italics: true, color: "888888", size: 18 } }));
      out.push(table(b.head || [], b.rows || [])); out.push(para("")); break;
    default: break;
  }
  return out;
}

// ── front matter ─────────────────────────────────────────────────────────────
const children = [];
children.push(new Paragraph({ heading: HeadingLevel.TITLE, children: [run(front.title)] }));
if (front.subtitle) children.push(para(front.subtitle, { run: { size: 26, color: "444444" } }));
children.push(para(`Captured ${man.captured} from the dashboard at ${man.url}.`, { run: { size: 18, color: "666666" } }));
(front.intro || []).forEach((t) => children.push(para(t)));

const totals = man.pages.map((pg) => [pg.title, String(pageWords(pg))]);
const grand = man.pages.reduce((a, pg) => a + pageWords(pg), 0) + words(man.banner);
children.push(new Paragraph({ heading: HeadingLevel.HEADING_2, children: [run("Words per page")] }));
children.push(para("Words a reader sees on each page, including labels and notes but not the text inside charts, maps and interactive data tables.", { run: { size: 18, color: "555555" } }));
children.push(table(["Page", "Words"], [...totals, ["Banner shown on every page", String(words(man.banner))], ["Total", String(grand)]], false));

if (front.changes && front.changes.length) {
  children.push(new Paragraph({ heading: HeadingLevel.HEADING_2, children: [run(front.changes_title || "What changed")] }));
  front.changes.forEach((c) => children.push(new Paragraph({ numbering: { reference: "bullets", level: 0 }, children: [run(c)], spacing: { after: 60 } })));
}
children.push(new Paragraph({ heading: HeadingLevel.HEADING_2, children: [run("Shown on every page")] }));
children.push(para(`Banner: ${man.banner}`));
if (man.footer) children.push(para(`Footer: ${man.footer}`));

// ── pages ────────────────────────────────────────────────────────────────────
man.pages.forEach((pg) => {
  children.push(new Paragraph({ children: [new PageBreak()] }));
  children.push(new Paragraph({ heading: HeadingLevel.HEADING_1, children: [run(pg.title)] }));
  children.push(para(`${pageWords(pg)} words on this page.`, { run: { size: 18, color: "666666" } }));
  pg.shots.forEach((png, si) => {
    if (si > 0) children.push(para("The same page with cross-hatching switched on (rank range view):", { run: { italics: true, size: 18, color: "555555" } }));
    (slices[png] || []).forEach(([jpg, w, h]) => {
      children.push(new Paragraph({
        alignment: AlignmentType.CENTER, spacing: { after: 80 },
        children: [new ImageRun({ type: "jpg", data: fs.readFileSync(path.join(capDir, "slices", jpg)),
                                  transformation: { width: IMG_W, height: Math.round(IMG_W * h / w) } })],
      }));
    });
  });
  children.push(new Paragraph({ heading: HeadingLevel.HEADING_2, children: [run("Text on this page")] }));
  (pg.blocks || []).forEach((b) => blockToDocx(b).forEach((x) => children.push(x)));
  const caps = pg.extra && pg.extra.captions;
  if (caps && caps.length) {
    children.push(new Paragraph({ heading: HeadingLevel.HEADING_3, children: [run("Caption under the map, for each view")] }));
    caps.forEach((c) => children.push(new Paragraph({ numbering: { reference: "bullets", level: 0 }, spacing: { after: 60 },
      children: [run(`${c.layer}: `, { bold: true }), run(c.caption)] })));
  }
});

const doc = new Document({
  creator: "Micronutrient Burden dashboard",
  title: front.title,
  styles: {
    default: { document: { run: { font: "Calibri", size: 21 } } },
    paragraphStyles: [
      { id: "Title", name: "Title", basedOn: "Normal", next: "Normal", run: { size: 40, bold: true, color: "0F7B8A" }, paragraph: { spacing: { after: 120 } } },
      { id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true, run: { size: 32, bold: true, color: "0F7B8A" }, paragraph: { spacing: { before: 120, after: 120 }, outlineLevel: 0 } },
      { id: "Heading2", name: "Heading 2", basedOn: "Normal", next: "Normal", quickFormat: true, run: { size: 26, bold: true, color: "333333" }, paragraph: { spacing: { before: 200, after: 100 }, outlineLevel: 1 } },
      { id: "Heading3", name: "Heading 3", basedOn: "Normal", next: "Normal", quickFormat: true, run: { size: 22, bold: true, color: "C8641E" }, paragraph: { spacing: { before: 160, after: 60 }, outlineLevel: 2 } },
    ],
  },
  numbering: { config: [{ reference: "bullets", levels: [{ level: 0, format: LevelFormat.BULLET, text: "•", alignment: AlignmentType.LEFT,
    style: { paragraph: { indent: { left: 540, hanging: 270 } } } }] }] },
  sections: [{
    properties: { page: { size: { width: 12240, height: 15840 }, margin: { top: 1080, bottom: 1080, left: 1440, right: 1440 } } },
    footers: { default: new Footer({ children: [new Paragraph({ alignment: AlignmentType.CENTER,
      children: [run("Page ", { size: 16, color: "888888" }), new TextRun({ children: [PageNumber.CURRENT], size: 16, color: "888888" })] })] }) },
    children,
  }],
});

Packer.toBuffer(doc).then((buf) => {
  fs.writeFileSync(outPath, buf);
  console.log(`wrote ${outPath}: ${man.pages.length} pages, ${grand} words of dashboard text`);
});
