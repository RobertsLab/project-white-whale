"""fetch-fastq with the gannet mirror: mirror copies, ENA fallback, mixed sources.

The local server from test_fastq.py serves ENA's files by name and the mirror
under /raw/<ACCESSION>/. Mirror copies are different bytes from ENA's, as
files converted from NCBI's .sra are. No network access beyond localhost.
"""
from __future__ import annotations

import csv
import hashlib

from click.testing import CliRunner
from test_fastq import PAYLOAD, _files, server  # noqa: F401 - server is a fixture

from whale.cli import main
from whale.fastq import apply_mirror, fetch_all, load_ledger, record_staged
from whale.mirror import LEDGER_COLUMNS

ACC = "PRJNA1"
MIRROR = {name: b"converted-" + data[::-1] for name, data in PAYLOAD.items()}


def _md5(b: bytes) -> str:
    return hashlib.md5(b).hexdigest()


def _ledger(rows):
    lines = ["\t".join(LEDGER_COLUMNS)]
    for r in rows:
        lines.append("\t".join(str(r.get(c, "")) for c in LEDGER_COLUMNS))
    return ("\n".join(lines) + "\n").encode()


def _mirror_row(name, *, md5=None, ena_md5=None):
    return {"run_accession": "SRR1", "mate": name[-10], "file": name, "url": "https://s3/sra/SRR1",
            "md5": md5 or _md5(MIRROR[name]), "bytes": len(MIRROR[name]), "source": "ncbi_sra",
            "ena_md5": ena_md5 or _md5(PAYLOAD[name]), "mirrored_on": "2026-10-23"}


def _serve_mirror(handler, base, rows, files=MIRROR):
    root = base.split("/vol1")[0]
    handler.routes[f"/raw/{ACC}/mirror_manifest.tsv"] = _ledger(rows)
    for name, data in files.items():
        handler.routes[f"/raw/{ACC}/{name}"] = data
    return f"{root}/raw/{ACC}"


def test_apply_mirror_counts():
    files = _files("http://ena")
    ledger = {files[0].file: _mirror_row(files[0].file),
              files[1].file: _mirror_row(files[1].file, ena_md5="0" * 32)}
    out, counts = apply_mirror(files, ledger, "http://gannet/raw/X/")
    assert counts == {"mirror": 1, "not_mirrored": 0, "stale": 1}
    assert out[0].source == "mirror" and out[0].url == f"http://gannet/raw/X/{files[0].file}"
    assert out[0].fallback is files[0] and out[0].ena_md5 == files[0].md5
    assert out[1] is files[1]
    assert apply_mirror(files, None, "http://g")[1] == {"mirror": 0, "not_mirrored": 2, "stale": 0}


def test_copies_from_mirror_and_records_source(server, tmp_path):  # noqa: F811
    base, handler = server
    url = _serve_mirror(handler, base, [_mirror_row(n) for n in PAYLOAD])
    files, counts = apply_mirror(_files(base), load_ledger(f"{url}/mirror_manifest.tsv"), url)
    assert counts["mirror"] == 2
    assert fetch_all(files, tmp_path / "raw")["downloaded"] == 2
    for name in PAYLOAD:
        assert (tmp_path / "raw" / name).read_bytes() == MIRROR[name]
    assert all(r.startswith(f"/raw/{ACC}/") for r in handler.requests if "manifest" not in r)

    record = tmp_path / "fastq_staged.tsv"
    record_staged(files, tmp_path / "raw", record, machine="raven", today="2026-10-23")
    rows = list(csv.DictReader(record.open(), delimiter="\t"))
    assert {r["source"] for r in rows} == {"mirror"}
    assert rows[0]["md5"] == _md5(MIRROR[rows[0]["file"]]) and rows[0]["ena_md5"] == _md5(PAYLOAD[rows[0]["file"]])


def test_broken_mirror_copy_falls_back_to_ena(server, tmp_path):  # noqa: F811
    base, handler = server
    name = "SRR1_1.fastq.gz"
    url = _serve_mirror(handler, base, [_mirror_row(name, md5="f" * 32)], {name: MIRROR[name]})
    files, _ = apply_mirror(_files(base), load_ledger(f"{url}/mirror_manifest.tsv"), url)
    res = fetch_all(files, tmp_path, log=lambda m: None)
    assert res == {"downloaded": 2, "skipped": 0, "failed": {}}
    assert (tmp_path / name).read_bytes() == PAYLOAD[name]
    assert (tmp_path / f"{name}.md5").read_text().split()[2] == "ena"


def test_file_staged_from_ena_is_not_refetched_once_mirrored(server, tmp_path):  # noqa: F811
    base, handler = server
    ena = _files(base)
    fetch_all(ena, tmp_path)
    url = _serve_mirror(handler, base, [_mirror_row(n) for n in PAYLOAD])
    files, _ = apply_mirror(ena, load_ledger(f"{url}/mirror_manifest.tsv"), url)
    handler.requests.clear()
    assert fetch_all(files, tmp_path)["skipped"] == 2
    assert handler.requests == []
    record = tmp_path / "rec.tsv"
    record_staged(files, tmp_path, record)
    assert {r["source"] for r in csv.DictReader(record.open(), delimiter="\t")} == {"ena"}


def test_partial_from_ena_is_not_resumed_from_mirror(server, tmp_path):  # noqa: F811
    base, handler = server
    name = "SRR1_1.fastq.gz"
    ena = _files(base)[0]
    (tmp_path / f"{name}.part").write_bytes(PAYLOAD[name][:500])
    (tmp_path / f"{name}.part.src").write_text(ena.https_url + "\n")
    url = _serve_mirror(handler, base, [_mirror_row(name)], {name: MIRROR[name]})
    files, _ = apply_mirror([ena], load_ledger(f"{url}/mirror_manifest.tsv"), url)
    handler.requests.clear()
    assert fetch_all(files, tmp_path)["downloaded"] == 1
    assert handler.requests == [f"/raw/{ACC}/{name}"]   # from the start, no Range
    assert (tmp_path / name).read_bytes() == MIRROR[name]


def test_no_ledger_is_none(server):  # noqa: F811
    base, _ = server
    assert load_ledger(f"{base.split('/vol1')[0]}/raw/NOPE/mirror_manifest.tsv") is None


def _cli_repo(tmp_path, base):
    repo = tmp_path / "repo"
    (repo / "registry" / "controlled_vocabularies").mkdir(parents=True)
    study = repo / "data" / "studies" / "X2020_Y"
    study.mkdir(parents=True)
    with open(study / "fastq_manifest.tsv", "w") as fh:
        fh.write("run_accession\tmate\tfile\turl\tmd5\tbytes\tsource\n")
        for f in _files(base):
            fh.write(f"SRR1\t1\t{f.file}\t{f.url}\t{f.md5}\t{f.bytes}\tfastq_ftp\n")
    return repo


def test_cli_sources(server, tmp_path, monkeypatch):  # noqa: F811
    base, handler = server
    repo = _cli_repo(tmp_path, base)
    monkeypatch.setenv("WHALE_REPO", str(repo))
    name = "SRR1_1.fastq.gz"
    url = _serve_mirror(handler, base, [_mirror_row(name)], {name: MIRROR[name]})
    args = ["fetch-fastq", "--study", "X2020_Y", "--mirror-url", url]

    res = CliRunner().invoke(main, [*args, "--source", "mirror", "--dest", str(tmp_path / "m")])
    assert res.exit_code != 0 and "not usable from the mirror" in res.output

    res = CliRunner().invoke(main, [*args, "--dest", str(tmp_path / "auto")])
    assert res.exit_code == 0, res.output
    assert "1 file(s) from the mirror, 1 not mirrored" in res.output
    assert (tmp_path / "auto" / name).read_bytes() == MIRROR[name]
    assert (tmp_path / "auto" / "SRR1_2.fastq.gz").read_bytes() == PAYLOAD["SRR1_2.fastq.gz"]

    res = CliRunner().invoke(main, [*args, "--source", "ena", "--dest", str(tmp_path / "ena")])
    assert res.exit_code == 0 and "Mirror" not in res.output
    assert (tmp_path / "ena" / name).read_bytes() == PAYLOAD[name]
