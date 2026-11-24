// Shared by the atlas pages. data/atlas.json is written by `whale build-atlas`
// (src/whale/atlas.py), in the Pages workflow, from the studies merged to main.

export const $ = (id) => document.getElementById(id);
export const esc = (s) => String(s ?? "").replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
export const fmt = (n) => (n == null ? "–" : Number(n).toLocaleString("en-US"));
export const pct = (x, d = 0) => (x == null ? "–" : `${(100 * x).toFixed(d)}%`);
export const num = (x, d = 2) => (x == null ? "–" : Number(x).toFixed(d));
export const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

export async function load() {
  const res = await fetch("data/atlas.json", { cache: "no-store" });
  if (!res.ok) throw new Error(`data/atlas.json: ${res.status}`);
  return res.json();
}

// The header's provenance line and the downloads list, common to both pages.
export function frame(data) {
  $("updated").innerHTML = `Built ${esc(data.generated_at.replace("T", " ").replace("Z", " UTC"))}` +
    (data.commit ? ` from <a href="${esc(data.repo_url)}/commit/${esc(data.commit)}">${esc(data.commit)}</a>` : "") +
    ` · ${esc(data.genome_assembly || "")} · rebuilt whenever a study merges`;
  const dl = $("downloads");
  if (dl) {
    dl.innerHTML = data.tables.map((t) => `<a href="data/${esc(t)}">${esc(t.replace("tables/", ""))}</a>`).join(" · ");
  }
  const ex = $("excluded");
  if (ex) {
    ex.innerHTML = data.excluded.length
      ? "Not pooled: " + data.excluded.map((e) => `${esc(e.study_id)} (${esc(e.reason)})`).join("; ") + "."
      : "";
  }
}

// Sequential blue, interpolated in sRGB between the ramp's two ends; t in [0, 1].
export function seq(t) {
  const a = hex(css("--seq-0")), b = hex(css("--seq-1"));
  t = Math.max(0, Math.min(1, t));
  return `rgb(${a.map((v, i) => Math.round(v + (b[i] - v) * t)).join(",")})`;
}
function hex(h) {
  h = h.replace("#", "");
  return [0, 2, 4].map((i) => parseInt(h.slice(i, i + 2), 16));
}

const tipEl = () => $("tip");
export function showTip(e, html) {
  const t = tipEl();
  t.innerHTML = html;
  t.style.display = "block";
  const x = Math.min(e.clientX + 14, window.innerWidth - t.offsetWidth - 8);
  const y = Math.min(e.clientY + 14, window.innerHeight - t.offsetHeight - 8);
  t.style.left = `${x}px`;
  t.style.top = `${y}px`;
}
export const hideTip = () => { tipEl().style.display = "none"; };

export function fail(err) {
  document.querySelector("main").insertAdjacentHTML("beforeend",
    `<div class="card empty">The atlas data could not be loaded (${esc(err.message)}). It is built by <code>whale build-atlas</code>.</div>`);
}
