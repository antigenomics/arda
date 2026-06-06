// arda markup-transfer hot path.
//
// Given a reference (target) sequence with known region coordinates (FR1-4 /
// CDR1-3) and an MMseqs2 alignment of a query to that target, project the
// reference region boundaries onto the query by walking the gapped aligned
// strings (qaln/taln) once.
//
// Coordinate convention: all public inputs/outputs are 1-based, closed — the
// AIRR convention. q_start/t_start are the 1-based alignment start positions in
// the query / target (mmseqs qstart/tstart). Regions not covered by the
// alignment (e.g. a 5'-truncated query) come back as (-1, -1).

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <string>
#include <utility>
#include <vector>

namespace py = pybind11;

using Interval = std::pair<int, int>;

// Project a single reference interval onto query coordinates. Kept as a simple
// primitive (used by unit tests); transfer_regions is the batched workhorse.
// ref_aln_offset / qry_aln_offset are 0-based (tstart-1 / qstart-1); ref_start /
// ref_end are 0-based inclusive. Returns 0-based inclusive (query_start, query_end).
static Interval project_region(const std::string &qaln, const std::string &taln,
                               int ref_aln_offset, int qry_aln_offset,
                               int ref_start, int ref_end) {
    int ref_pos = ref_aln_offset;
    int qry_pos = qry_aln_offset;
    int q_start = -1, q_end = -1;
    const size_t n = std::min(qaln.size(), taln.size());
    for (size_t i = 0; i < n; ++i) {
        const bool ref_gap = (taln[i] == '-');
        const bool qry_gap = (qaln[i] == '-');
        if (!ref_gap && ref_pos >= ref_start && ref_pos <= ref_end && !qry_gap) {
            if (q_start < 0) q_start = qry_pos;
            q_end = qry_pos;
        }
        if (!ref_gap) ++ref_pos;
        if (!qry_gap) ++qry_pos;
        if (ref_pos > ref_end + 1) break;
    }
    return {q_start, q_end};
}

// Batched projection: walk the alignment once and project every region.
// region_starts/region_ends are 1-based closed coordinates on the TARGET
// (reference scaffold). Returns one 1-based closed (q_start, q_end) per region,
// or (-1, -1) where a region is not covered by the alignment.
//
// Indel semantics: query residues inserted *within* a region (gap in target)
// fall between the region's first and last aligned query positions, so they are
// included by the [q_start, q_end] span. Deleted reference positions (gap in
// query) simply contribute no query base. A region boundary landing on a gap
// clamps to the nearest aligned query base.
static std::vector<Interval> transfer_regions(
    const std::string &qaln, const std::string &taln,
    int q_start, int t_start,
    const std::vector<int> &region_starts,
    const std::vector<int> &region_ends) {

    const size_t R = region_starts.size();
    std::vector<int> rs(R), re(R), qs(R, -1), qe(R, -1);
    int max_re = -1;
    for (size_t k = 0; k < R; ++k) {
        rs[k] = region_starts[k] - 1;  // -> 0-based inclusive
        re[k] = region_ends[k] - 1;
        if (re[k] > max_re) max_re = re[k];
    }

    int ref_pos = t_start - 1;
    int qry_pos = q_start - 1;
    const size_t n = std::min(qaln.size(), taln.size());
    for (size_t i = 0; i < n; ++i) {
        const bool ref_gap = (taln[i] == '-');
        const bool qry_gap = (qaln[i] == '-');
        if (!ref_gap && !qry_gap) {
            for (size_t k = 0; k < R; ++k) {
                if (ref_pos >= rs[k] && ref_pos <= re[k]) {
                    if (qs[k] < 0) qs[k] = qry_pos;
                    qe[k] = qry_pos;
                }
            }
        }
        if (!ref_gap) ++ref_pos;
        if (!qry_gap) ++qry_pos;
        if (ref_pos > max_re + 1) break;  // past every region
    }

    std::vector<Interval> out(R);
    for (size_t k = 0; k < R; ++k) {
        out[k] = (qs[k] < 0) ? Interval{-1, -1}
                             : Interval{qs[k] + 1, qe[k] + 1};  // -> 1-based closed
    }
    return out;
}

PYBIND11_MODULE(_markup, m) {
    m.doc() = "arda markup-transfer hot path (C++/pybind11)";
    m.attr("__version__") = "0.2.0";
    m.def("project_region", &project_region,
          py::arg("qaln"), py::arg("taln"), py::arg("ref_aln_offset"),
          py::arg("qry_aln_offset"), py::arg("ref_start"), py::arg("ref_end"),
          "Project a single 0-based inclusive reference interval onto 0-based "
          "inclusive query coordinates. Returns (-1,-1) if no overlap.");
    m.def("transfer_regions", &transfer_regions,
          py::arg("qaln"), py::arg("taln"), py::arg("q_start"), py::arg("t_start"),
          py::arg("region_starts"), py::arg("region_ends"),
          "Project multiple 1-based closed reference (target) intervals onto "
          "1-based closed query coordinates in a single alignment walk. Returns "
          "one (q_start, q_end) per region; (-1,-1) where uncovered.");
}
