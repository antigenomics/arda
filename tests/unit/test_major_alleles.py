"""Major-only selection preserves gene ambiguity and searches actual templates."""
from arda import select_major_alleles


def test_major_order_numeric_fallback_and_starless_d():
    names = ['V1*10', 'V1*02', 'V2*02', 'V2*01', 'D1', 'V3*09']
    assert select_major_alleles(names) == {'V1': 'V1*02', 'V2': 'V2*01', 'D1': 'D1', 'V3': 'V3*09'}
    assert select_major_alleles(reversed(names)) == select_major_alleles(names)


def test_major_map_uses_only_eligible_available_anchors(tmp_path, monkeypatch):
    from arda import major_allele_map
    from arda import paths
    monkeypatch.setattr(paths, 'vdj_dir', lambda _: tmp_path)
    (tmp_path/'markup.tsv').write_text('locus\tv_call\tj_call\nTRB\tV1*01,V1*02,V2*10,V2*02\tJ1*01\n')
    (tmp_path/'cdr3_anchors.tsv').write_text('locus\tallele\tstatus\tfunctionality\nTRB\tV1*01\tpartial\tF\nTRB\tV1*02\tok\tF\nTRB\tV2*10\tok\tF\nTRB\tV2*02\tok\tF\nTRB\tJ1*01\tok\tF\nTRB\tV3*01\tok\tF\n')
    (tmp_path/'d_germlines.fasta').write_text('>TRB|D1\nACGT\n>TRB|D2*02\nACGT\n>TRB|D2*10\nAAAA\n')
    assert major_allele_map('human','TRB') == {'V1':'V1*02','V2':'V2*02','J1':'J1*01','D1':'D1','D2':'D2*02'}


def test_major_map_keeps_genes_with_only_partial_alignment_targets(tmp_path, monkeypatch):
    from arda import major_allele_map, paths
    monkeypatch.setattr(paths, 'vdj_dir', lambda _: tmp_path)
    (tmp_path/'markup.tsv').write_text('locus\tv_call\tj_call\nTRB\tV1*01,V1*02,V2*01\tJ1*01\n')
    (tmp_path/'cdr3_anchors.tsv').write_text('locus\tallele\tstatus\tfunctionality\nTRB\tV1*01\tpartial\tF\nTRB\tV1*02\tok\tF\nTRB\tV2*01\tpartial\tF\nTRB\tJ1*01\tok\tF\n')
    (tmp_path/'d_germlines.fasta').write_text('')
    assert major_allele_map('human') == {'V1':'V1*02','V2':'V2*01','J1':'J1*01'}


def test_major_reference_searches_exact_selected_templates(tmp_path, monkeypatch):
    import csv
    from arda import paths
    from arda.annotate import major
    from arda.annotate.reference import Reference, REGIONS

    base=tmp_path/'base'
    base.mkdir()
    monkeypatch.setattr(paths, 'vdj_dir', lambda _: base)
    monkeypatch.setattr(major, 'data_dir', lambda: tmp_path/'cache')
    rows=[]
    for i, allele in enumerate(('V1*01','V1*02','V1*10','V2*01')):
        rows.append(dict(scaffold_id=f's{i}',locus='TRB',v_call=allele,j_call='J1*01',c_call='',
                         v_sequence_end=12,j_sequence_start=13,vj_end=24,
                         **{f'{r}_{end}': -1 for r in REGIONS for end in ('start','end')}))
    for name in ('markup.tsv','markup.aa.tsv'):
        with (base/name).open('w') as handle:
            writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
    for name in ('alleles.fasta','alleles.aa.fasta'):
        (base/name).write_text(''.join(f'>s{i}\n{seq}\n' for i,seq in enumerate(('A'*24,'C'*12+'T'*12,'G'*24,'T'*24))))
    with (base/'combinations.tsv').open('w') as handle:
        handle.write('scaffold_id\tlocus\tv_calls\tj_calls\tn_pad\n')
        for row in rows:
            handle.write(f"{row['scaffold_id']}\tTRB\t{row['v_call']}\tJ1*01\t0\n")
    (base/'cdr3_anchors.tsv').write_text('locus\tallele\tstatus\tfunctionality\nTRB\tV1*01\tpartial\tF\nTRB\tV1*02\tok\tF\nTRB\tV1*10\tok\tF\nTRB\tV2*01\tok\tF\nTRB\tJ1*01\tok\tF\n')
    (base/'d_germlines.fasta').write_text('>TRB|D1*02\nACGT\n>TRB|D1*10\nAAAA\n>TRB|D2\nGGGG\n')
    ref=Reference('human','nt',base/'alleles.fasta',{}, {'TRB':[('D1*02','ACGT'),('D1*10','AAAA'),('D2','GGGG')]})
    result=major.major_reference(ref)
    assert {key for key in result.entries if '|' not in key} == {'s1','s3'}
    assert '>s1\nCCCCCCCCCCCCTTTTTTTTTTTT\n' in result.target_fasta.read_text()
    assert result.d_germlines['TRB'] == [('D1*02','ACGT'),('D2','GGGG')]
    assert {e.v_call for key,e in result.entries.items() if '|' not in key} == {'V1*02','V2*01'}
    assert (base/'alleles.fasta').read_text().startswith('>s0\n'+'A'*24)
    assert major.major_reference(ref).target_fasta == result.target_fasta


def test_public_search_defaults_are_major_only():
    import inspect
    from arda import alignment_candidates, annotate_sequences
    from arda.annotate.mapper import annotate_records, annotate_file
    from arda.annotate.reference import load_reference
    from arda.rnaseq.map import map_rnaseq
    for function in (alignment_candidates, annotate_sequences, annotate_records, annotate_file, load_reference, map_rnaseq):
        assert inspect.signature(function).parameters['major_only'].default is True


def test_name_selector_import_stays_lightweight():
    import subprocess
    import sys
    result = subprocess.run([sys.executable, '-c', "import sys; from arda import select_major_alleles; assert select_major_alleles(['V*10','V*02']) == {'V':'V*02'}; assert 'polars' not in sys.modules"], capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stderr
