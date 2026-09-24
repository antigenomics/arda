"""``arda qc report`` — the self-contained dashboard.

The guarantee this file exists for is **no external reference**: the page must open on an
air-gapped login node, off a USB stick, or as an email attachment, months after the results
directory is gone. That breaks silently — a CDN `<script>` renders fine on the machine that made
it and blank everywhere else — so it is asserted directly rather than trusted.
"""

from __future__ import annotations

import json
import re
from html.parser import HTMLParser

import pytest

from arda.qcreport import build, render

_DOC = {
    "samples": ["S0", "S1"],
    "labels": {"S0": {"project": "TRIAL9", "batch": "RUN1"},
               "S1": {"project": "TRIAL9", "batch": "RUN1"}},
    "stats": {
        "S0": {"sample": {"": {"reads": 210, "clonotypes": 19}},
               "run": {"map": {"mapped_fraction": 0.6364}},
               "junction_aa_len": {"IGH:15": {"reads": 3}}},
        "S1": {"sample": {"": {"reads": 103, "clonotypes": 16}},
               "run": {"map": {"mapped_fraction": 0.3121}},
               "junction_aa_len": {"IGH:15": {"reads": 2}}},
    },
    "outliers": [], "z_flag": 3.5, "min_group": 5,
}


@pytest.fixture
def page():
    return build(_DOC, title="test cohort")


def test_the_page_fetches_nothing(page):
    """No script src, no stylesheet link, no CSS url(), no fetch. The SVG namespace URI is not a
    network reference — browsers never retrieve it — so it is excluded explicitly rather than by
    a blanket 'no http' rule that would be a lie."""
    assert not re.search(r"<script[^>]*\ssrc=", page)
    assert not re.search(r"<link[^>]*\shref=", page)
    assert "@import" not in page
    assert not re.search(r"url\(\s*['\"]?https?:", page)
    assert not re.search(r"\bfetch\s*\(|XMLHttpRequest|importScripts", page)
    external = [m for m in re.findall(r"https?://[^\s'\"<>]+", page)
                if m != "http://www.w3.org/2000/svg"]
    assert external == []


def test_the_data_round_trips_out_of_the_inlined_block(page):
    blob = re.search(r'<script type="application/json" id="qc-data">(.*?)</script>',
                     page, re.S).group(1)
    assert json.loads(blob.replace("<\\/", "</")) == _DOC


def test_a_closing_script_tag_inside_the_data_cannot_end_the_block_early():
    """A sample id or a file path containing `</script` would truncate the page at that byte."""
    doc = dict(_DOC, samples=["</script><b>x"], labels={}, stats={"</script><b>x": {}})
    page = build(doc)
    assert "</script><b>x" not in page.split('id="qc-data"')[1].split("</script>")[0]
    assert page.count('<script type="application/json"') == 1


def test_the_page_is_well_formed_html_and_carries_every_panel(page):
    class Parser(HTMLParser):
        def __init__(self):
            super().__init__()
            self.ids, self.stack = set(), []

        def handle_starttag(self, tag, attrs):
            d = dict(attrs)
            if "id" in d:
                self.ids.add(d["id"])
            if tag not in ("meta", "br", "input", "img", "link", "hr"):
                self.stack.append(tag)

        def handle_endtag(self, tag):
            if self.stack and self.stack[-1] == tag:
                self.stack.pop()

    p = Parser()
    p.feed(page)
    assert p.stack == [], f"unclosed: {p.stack}"
    assert {"qc-data", "table", "metric", "dist", "prov", "f-project", "f-batch",
            "m-pick", "d-scope"} <= p.ids


def test_a_single_sample_stats_json_renders_as_a_cohort_of_one(tmp_path):
    """One run is a cohort of one — the same page, not a second renderer."""
    stats = tmp_path / "PT01.stats.json"
    stats.write_text(json.dumps({"arda_version": "2.23.0", "rows": 2,
                                 "stats": {"sample": {"": {"reads": 7}}}}))
    out = render(stats, tmp_path / "PT01.qc.html")
    page = out.read_text()
    blob = json.loads(re.search(r'id="qc-data">(.*?)</script>', page, re.S).group(1))
    assert blob["samples"] == ["PT01"]
    assert blob["stats"]["PT01"]["sample"][""]["reads"] == 7


def test_a_value_naming_a_later_placeholder_is_not_substituted_into():
    """Chained `.replace()` substitutes into what it has already substituted, so a sample id of
    `__JS__` would have the whole script spliced in where its name belongs."""
    doc = dict(_DOC, samples=["__JS__"], labels={"__JS__": {"project": "", "batch": ""}},
               stats={"__JS__": {"sample": {"": {"reads": 3}}}})
    page = build(doc, title="__PAYLOAD__")
    blob = page.split('id="qc-data">')[1].split("</script>")[0]
    assert "const DOC" not in blob
    assert page.count("const DOC = JSON.parse") == 1
    assert '"__JS__"' in blob


def test_the_title_is_escaped():
    """It comes from a filename, and a filename is not trusted markup."""
    page = build(_DOC, title="<script>alert(1)</script>")
    assert "<h1>&lt;script&gt;alert(1)&lt;/script&gt;</h1>" in page
    assert "<title>&lt;script&gt;" in page
