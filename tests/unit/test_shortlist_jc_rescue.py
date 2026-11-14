"""A contradictory shortlist is not a substitute for full-reference competition."""
from types import SimpleNamespace
import pytest
from arda.annotate import mapper


def test_jc_winner_over_shortlisted_v_requires_full_reference(tmp_path, monkeypatch):
    monkeypatch.setattr(mapper.mmseqs, "search", lambda *a, **k: None)
    monkeypatch.setattr(mapper.mmseqs, "convertalis", lambda *a, **k: None)
    monkeypatch.setattr(mapper, "_db_keys", lambda *a: {"q": 0})
    monkeypatch.setattr(mapper, "_segment_rows", lambda *a: [
        dict(query="q", target="V|V", bits=100), dict(query="q", target="J|J", bits=80)])
    monkeypatch.setattr(mapper, "_align_implied", lambda *a, **k:
                        ({"q": dict(target="JC", bits=139)} if k.get("side")=="J"
                         else {"q": dict(target="short_VJ", bits=63)}, set()))
    rescued=[]
    def full(*a, **k):
        rescued.extend(a[2]);return {"q": dict(target="full_VJ", bits=200)}
    monkeypatch.setattr(mapper, "_full_rescue", full)
    ref=SimpleNamespace(jc_combinations=lambda:{("J", "C"):"JC"})
    best,_=mapper._segment_best_hits(tmp_path/"q",tmp_path/"s",tmp_path/"f",tmp_path,
        ref,threads=1,sensitivity=7,mm_strand=2,max_seqs=300,kmer=12,
        search_type=3,seqs={"q":"ACGT"},combos={("V","J"):"short_VJ"})
    assert rescued==["q"]
    assert best["q"]["target"]=="full_VJ"


@pytest.mark.parametrize("bounds", [dict(cdr3_start=4, fwr4_start=10),
                                   dict(cdr3_start=4), dict(fwr4_start=10)])
def test_rejected_extended_anchor_gets_original_extension_batch(tmp_path, monkeypatch, bounds):
    searches=[]
    monkeypatch.setattr(mapper.mmseqs,"createdb",lambda *a,**k:None)
    def search(*a,**k): searches.append(k.get("zdrop",200))
    monkeypatch.setattr(mapper.mmseqs,"search",search)
    monkeypatch.setattr(mapper.mmseqs,"top_hit",lambda *a:a[1])
    monkeypatch.setattr(mapper.mmseqs,"convertalis",lambda *a,**k:None)
    monkeypatch.setattr(mapper,"_best_hits",lambda *a:{"q":dict(target="VJ",qstart=1,qend=12)})
    def transfer(qid,*a,**k):
        return dict(sequence_id=qid,**bounds,
                    junction="TGTAAATGG" if searches[-1]==40 else "")
    monkeypatch.setattr(mapper,"transfer_hit",transfer)
    ref=SimpleNamespace(get=lambda *a:SimpleNamespace(locus="TRB"),d_germlines={},anchors={})
    out=mapper._annotate_chunk([("q","TGTAAATGGAAA")],ref,tmp_path/"db","nt",
                              threads=1,sensitivity=7,mm_strand=2)
    assert out[0]["junction"]=="TGTAAATGG"
    assert searches==[200,40]


def test_no_segment_hit_still_competes_against_full_reference(tmp_path, monkeypatch):
    """Weak separate V/J evidence must not discard a valid whole-scaffold alignment."""
    monkeypatch.setattr(mapper.mmseqs, 'search', lambda *a, **k: None)
    monkeypatch.setattr(mapper.mmseqs, 'convertalis', lambda *a, **k: None)
    monkeypatch.setattr(mapper, '_db_keys', lambda *a: {'q': 0})
    monkeypatch.setattr(mapper, '_segment_rows', lambda *a: [])
    rescued = []

    def full(*a, **k):
        rescued.extend(a[2])
        return {'q': dict(target='full_VJ', bits=100)}

    monkeypatch.setattr(mapper, '_full_rescue', full)
    ref = SimpleNamespace(jc_combinations=lambda: {})
    best, report = mapper._segment_best_hits(
        tmp_path/'q', tmp_path/'s', tmp_path/'f', tmp_path, ref,
        threads=1, sensitivity=7, mm_strand=2, max_seqs=300, kmer=None,
        search_type=3, seqs={'q': 'ACGT'}, combos={})
    assert rescued == ['q']
    assert best['q']['target'] == 'full_VJ'
    assert report['rescued'] == 1
    assert report['no_segment_hit'] == 0  # no unsearched residual
    assert report['reasons']['no_segment_hit'] == 1


def test_rejected_shortlist_gets_standard_full_reference_before_old_extension(tmp_path, monkeypatch):
    """A rejected nominee has not yet competed under the ordinary full-search settings."""
    monkeypatch.setattr(mapper.mmseqs, 'createdb', lambda *a, **k: None)
    monkeypatch.setattr(mapper, '_segment_best_hits', lambda *a, **k:
                        ({'q': dict(target='shortlist', qstart=1, qend=12)}, {}))
    searches = []
    monkeypatch.setattr(mapper.mmseqs, 'search', lambda *a, **k: searches.append(k['zdrop']))
    monkeypatch.setattr(mapper.mmseqs, 'top_hit', lambda *a: a[1])
    monkeypatch.setattr(mapper.mmseqs, 'convertalis', lambda *a, **k: None)
    monkeypatch.setattr(mapper, '_best_hits', lambda *a:
                        {'q': dict(target='full', qstart=1, qend=12)})
    def transfer(qid, seq, hit, *a, **k):
        return dict(sequence_id=qid, cdr3_start=4,
                    junction='TGTAAATGG' if hit['target']=='full' and searches[-1]==200 else '')
    monkeypatch.setattr(mapper, 'transfer_hit', transfer)
    ref = SimpleNamespace(get=lambda *a: SimpleNamespace(locus='TRB'), d_germlines={}, anchors={})
    out = mapper._annotate_chunk([('q','TGTAAATGGAAA')], ref, tmp_path/'db', 'nt',
                                threads=1, sensitivity=7, mm_strand=2, segment_db=tmp_path/'segments')
    assert out[0]['junction'] == 'TGTAAATGG'
    assert searches == [200]
