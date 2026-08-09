// arda denoise: the hot paths of quality-aware clonotype error correction.
//
// Why this exists
// ---------------
// Stage 2 decides which clonotypes are real. Two of its steps scale with the READ count or with
// the SQUARE of the clonotype count, and both were Python loops:
//
//   * per-read mean junction Phred -- one pass over every base of every junction-bearing read
//     (310,559 reads x ~48 nt on a MIGEC library, 88,697 x ~45 on an IGH repertoire);
//   * the wide-radius parent search -- for each clonotype, the most abundant neighbour within
//     `max_subs` substitutions. seqtree's index handles the k <= 3 case the abundance model needs,
//     but the quality-directed rescue below searches k up to ~12, where an edit-bounded index
//     degenerates and a direct length-bucketed Hamming scan with early exit is both simpler and
//     faster. On IGH_repertoire that is 31,943 clonotypes -- 5.1e8 ordered pairs.
//
// What the measurements say, because they set every rule here
// -----------------------------------------------------------
// (benchmark repo: results/round22/junction_key_and_ladder.md, results/round22/jurkat_errors/)
//
//   * **The abundance model is only justified for k <= 2, with 3 as headroom.** On Jurkat TRB,
//     every one of the 82 two-substitution variants has an OBSERVED one-substitution intermediate
//     on its path to the parent, and the observed count (82) is within an order of magnitude of
//     the binomial prediction (11.2). That is a ladder, and chain collapse walks it.
//   * **At k >= 4 there is no ladder -- 0 of 13 (k=4), 0 of 18 (k=5), 0 of 14 (k>=6) have any
//     observed intermediate**, and the independent-error model predicts 0.0019 clonotypes at k=4
//     where 13 are observed. Those are not accumulated substitutions. They are single bad reads:
//     median mean junction Phred falls monotonically 31.4 (k=1) -> 16.5 (k=8), and the k >= 5
//     class is 100 % sub-Q30 while the dominant clone is 5.9 %.
//   * So a far neighbour may only be collapsed on QUALITY evidence, never on abundance alone --
//     which is why `nearest_more_abundant` takes a candidate mask rather than being turned loose.
//
// Read conservation
// -----------------
// Nothing here deletes anything. Both functions return indices; routing reads to a parent is the
// caller's job, and the caller's invariant is that the sum of duplicate_count does not fall.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace py = pybind11;

namespace {

// Per-read mean Phred over a quality string, Phred+33. Returns -1.0 for a read whose quality is
// missing or whose length disagrees with its junction -- ABSENT EVIDENCE, which the caller must
// treat as "no information", never as "bad". A length disagreement is the one corruption nothing
// downstream can detect (a right-length string from the wrong strand or offset), so it is refused
// here rather than averaged.
std::vector<double> mean_phred(const std::vector<std::string>& junctions,
                               const std::vector<std::string>& quals) {
    if (junctions.size() != quals.size())
        throw std::invalid_argument("junctions and quals must be the same length");
    std::vector<double> out(junctions.size(), -1.0);
    for (size_t i = 0; i < junctions.size(); ++i) {
        const std::string& q = quals[i];
        if (q.empty() || q.size() != junctions[i].size()) continue;
        uint64_t sum = 0;
        for (unsigned char c : q) sum += static_cast<uint64_t>(c) - 33u;
        out[i] = static_cast<double>(sum) / static_cast<double>(q.size());
    }
    return out;
}

// Fraction of a read's junction bases below `cut` (Phred+33 assumed). -1.0 on missing/mismatched
// quality, same contract as mean_phred.
std::vector<double> frac_below(const std::vector<std::string>& junctions,
                               const std::vector<std::string>& quals, int cut) {
    if (junctions.size() != quals.size())
        throw std::invalid_argument("junctions and quals must be the same length");
    const int thr = cut + 33;
    std::vector<double> out(junctions.size(), -1.0);
    for (size_t i = 0; i < junctions.size(); ++i) {
        const std::string& q = quals[i];
        if (q.empty() || q.size() != junctions[i].size()) continue;
        size_t n = 0;
        for (unsigned char c : q) n += (static_cast<int>(c) < thr);
        out[i] = static_cast<double>(n) / static_cast<double>(q.size());
    }
    return out;
}

// For each clonotype flagged in `candidates`, the index of the most abundant clonotype within
// `max_subs` substitutions whose count is at least `min_ratio` times its own; -1 if none.
//
// Only EQUAL-LENGTH sequences are compared: a length change is an indel, a different event with a
// different rate, and pretending otherwise silently aligns two junctions that do not correspond.
// Sequences are bucketed by length so the scan never crosses lengths, and the per-pair loop exits
// as soon as the mismatch budget is blown -- at k = 3 over a 48 nt junction that is typically 4
// base comparisons, not 48.
//
// Ties are broken by count, then by the sequence itself, so the result does not depend on input
// order. That is not cosmetic: `correct` was nondeterministic for exactly this class of reason.
std::vector<int64_t> nearest_more_abundant(const std::vector<std::string>& seqs,
                                           const std::vector<int64_t>& counts,
                                           const std::vector<bool>& candidates,
                                           int max_subs, double min_ratio) {
    const size_t n = seqs.size();
    if (counts.size() != n || candidates.size() != n)
        throw std::invalid_argument("seqs, counts and candidates must be the same length");
    std::vector<int64_t> parent(n, -1);

    std::unordered_map<size_t, std::vector<size_t>> by_len;
    for (size_t i = 0; i < n; ++i) by_len[seqs[i].size()].push_back(i);

    for (size_t ci = 0; ci < n; ++ci) {
        if (!candidates[ci]) continue;
        const std::string& a = seqs[ci];
        if (a.empty()) continue;
        const auto it = by_len.find(a.size());
        if (it == by_len.end()) continue;
        // A parent must be strictly more abundant, and by at least `min_ratio`. Without the
        // strictness two clonotypes of equal count can each be the other's parent and the caller's
        // root walk never terminates.
        const double need = static_cast<double>(counts[ci]) * min_ratio;
        int64_t best = -1, best_count = -1;
        for (size_t nj : it->second) {
            if (nj == ci) continue;
            const int64_t cj = counts[nj];
            if (cj <= counts[ci] || static_cast<double>(cj) < need) continue;
            if (cj < best_count) continue;                 // cannot win; skip the compare entirely
            const std::string& b = seqs[nj];
            int mm = 0;
            bool ok = true;
            for (size_t p = 0; p < a.size(); ++p) {
                if (a[p] != b[p] && ++mm > max_subs) { ok = false; break; }
            }
            if (!ok || mm == 0) continue;
            if (cj > best_count || (cj == best_count && best >= 0 &&
                                    b < seqs[static_cast<size_t>(best)])) {
                best = static_cast<int64_t>(nj);
                best_count = cj;
            }
        }
        parent[ci] = best;
    }
    return parent;
}

// Indices of candidates that CONTAIN `segment` verbatim.
//
// The tie-list primitive. A read aligned to germline `G` over `[gstart, gend)` is explained exactly
// as well by any other germline containing that same stretch -- the read carries no base that
// distinguishes them, so naming only `G` is a claim the data does not support. Measured on Ramos:
// 59 of 60 disagreements with IgBLAST are pairs where BOTH germlines match at identity 1.0000 over
// 63-70 nt, and arda emitted 0 tie lists in 504 calls against IgBLAST's 11.68 %.
//
// Substring search rather than a position-wise compare, because ungapped germline sequences are
// not in a common coordinate system -- the same stretch sits at different offsets in two alleles,
// and comparing `other[gstart:gend]` would miss it. `std::search` is fine here: the needle is
// ~60-200 nt, the haystacks a few hundred, and the caller memoises on (allele, gstart, gend), which
// on an amplicon collapses hundreds of thousands of reads onto a few thousand distinct spans.
//
// ⛔ Returns candidates in the order given, so the caller's ordering (and therefore the emitted
// call string) is deterministic.
std::vector<int64_t> containing(const std::string& segment,
                                const std::vector<std::string>& candidates) {
    std::vector<int64_t> out;
    if (segment.empty()) return out;
    for (size_t i = 0; i < candidates.size(); ++i) {
        const std::string& c = candidates[i];
        if (c.size() < segment.size()) continue;
        if (std::search(c.begin(), c.end(), segment.begin(), segment.end()) != c.end())
            out.push_back(static_cast<int64_t>(i));
    }
    return out;
}

// Substitution distance to a reference sequence, -1 when the lengths differ. Used to build the
// ladder/cliff diagnostic the framework reports, so the evidence for a collapse is inspectable
// rather than asserted.
std::vector<int> subs_to(const std::vector<std::string>& seqs, const std::string& ref) {
    std::vector<int> out(seqs.size(), -1);
    for (size_t i = 0; i < seqs.size(); ++i) {
        if (seqs[i].size() != ref.size()) continue;
        int mm = 0;
        for (size_t p = 0; p < ref.size(); ++p) mm += (seqs[i][p] != ref[p]);
        out[i] = mm;
    }
    return out;
}

}  // namespace

PYBIND11_MODULE(_denoise, m) {
    m.doc() = "Hot paths of quality-aware clonotype error correction (see denoise.cpp header).";
    m.attr("__version__") = "0.2.0";
    m.def("mean_phred", &mean_phred, py::arg("junctions"), py::arg("quals"),
          "Per-read mean Phred over the junction; -1.0 when the quality is absent or its length "
          "disagrees with the junction (absent evidence, not bad evidence).");
    m.def("frac_below", &frac_below, py::arg("junctions"), py::arg("quals"), py::arg("cut"),
          "Per-read fraction of junction bases below `cut`; -1.0 on absent/mismatched quality.");
    m.def("nearest_more_abundant", &nearest_more_abundant,
          py::arg("seqs"), py::arg("counts"), py::arg("candidates"),
          py::arg("max_subs"), py::arg("min_ratio") = 1.0,
          "For each flagged clonotype, the most abundant equal-length neighbour within `max_subs` "
          "substitutions and at least `min_ratio` times its count; -1 if none. Deterministic.");
    m.def("containing", &containing, py::arg("segment"), py::arg("candidates"),
          "Indices of candidates containing `segment` verbatim -- the germlines a read aligned "
          "over that stretch cannot distinguish from the one it was called against.");
    m.def("subs_to", &subs_to, py::arg("seqs"), py::arg("ref"),
          "Substitution distance from each sequence to `ref`; -1 when the lengths differ.");
}
