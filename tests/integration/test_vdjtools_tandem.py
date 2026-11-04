"""Optional cross-package contract: run with vdjtools installed; no model refitting."""
import pytest

pytest.importorskip('vdjtools')

from arda.annotate.dmap import map_d_junction
from vdjtools.model import load_bundled
from vdjtools.model.generate import generate


@pytest.mark.parametrize('locus', ['TRB', 'TRD', 'IGH'])
def test_vdjtools_tandem_emissions_and_single_d_control(locus):
    rows = generate(load_bundled(locus), 1000, seed=157, engine='native', threads=1)
    tandem = rows.filter(rows['d2_call'].fill_null('') != '')
    assert tandem.height > 0
    recovered = 0
    for row in tandem.iter_rows(named=True):
        args = row['junction_nt'], row['v_call'], row['j_call']
        call = map_d_junction(*args)
        single = map_d_junction(*args, map_dd=False)
        assert not single.is_dd and not single.d2_call
        if call.is_dd:
            recovered += 1
            assert call.d_sequence_end < call.d2_sequence_start
            junction = row['junction_nt']
            d1 = junction[call.d_sequence_start-1:call.d_sequence_end]
            d2 = junction[call.d2_sequence_start-1:call.d2_sequence_end]
            interior = junction[call.v_sequence_end:call.j_sequence_start-1]
            assert call.np1 + d1 + call.np2 + d2 + call.np3 == interior
    # Short, trimmed TRB D segments can fail the established evidence gate.
    if locus in ('TRD', 'IGH'):
        assert recovered > 0
    print(f'{locus}: generated tandem={tandem.height}, aligned tandem={recovered}, draws=1000')
