// arda segment mapper: answer "which V and which J" structurally, instead of asking a general
// homology search engine.
//
// Why this exists
// ---------------
// The two-pass path runs `mmseqs search` against a 924-target segment reference purely to learn,
// per read, its best V allele and its best J allele with coordinates -- `_SEGMENT_FORMAT` is
// `query,target,bits,qstart,qend,tstart`, no cigar, no backtrace. Measured (100,000 amplicon reads,
// 8 threads, results/round10/aligner_options.md):
//
//     mmseqs search vs the segment DB      2770 ms    0.036 M reads/s
//     arda's _prefilter seed scan             4 ms    25.1  M reads/s     ~700x
//
// Both walk every k-mer of every read against the same reference. The prefilter simply records a
// boolean and throws away WHICH target the seed came from and at what offset -- which is exactly
// what the segment pass exists to recover. This module keeps it.
//
// What makes that sound here, and would not be sound for a general aligner:
//
//   * matches are to GERMLINE, at ~99 % identity, not remote homology. Measured: mmseqs'
//     `--exact-kmer-matching 1` loses zero hits on this reference, so its similarity-based k-mer
//     lookup is buying sensitivity nobody needs;
//   * segment assignment needs UNGAPPED extension only. Germline V and J carry no indels relative
//     to a read except sequencing error and IG SHM, and the gapped alignment is still done once,
//     later, against the full V+pad+J scaffold;
//   * the junction is NON-TEMPLATED and aligns to nothing. V and J are extended independently and
//     the junction is whatever lies between them -- there is no N-pad to gap through;
//   * the reference is tiny (~236 kb) and fixed, so the index is ~2 MB and cache-resident.
//
// ⛔ This does NOT replace the second pass. The final call is still made by mmseqs against the
//    full scaffold, so this only has to NOMINATE the right candidates. That is the regression
//    contract: the AIRR output must not move, not that the scores match mmseqs'.
//
// Index layout
// ------------
// CSR over sorted unique k-mer codes: `codes_[i]` owns `postings_[starts_[i] .. starts_[i+1])`.
// A posting packs (target, offset) into one uint32. Only the FORWARD strand of each target is
// indexed; the read is scanned forwards and as its reverse complement, which keeps the diagonal
// arithmetic uniform and halves the index against the prefilter's both-strands layout.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace py = pybind11;

namespace {

std::array<int8_t, 256> make_base_idx() {
    std::array<int8_t, 256> t{};
    t.fill(-1);
    t['A'] = t['a'] = 0; t['C'] = t['c'] = 1; t['G'] = t['g'] = 2; t['T'] = t['t'] = 3;
    return t;
}
const std::array<int8_t, 256> BASE_IDX = make_base_idx();

constexpr int PRE_K = 12;                                  // first-level bitset width
constexpr size_t PRE_BITS = size_t(1) << (2 * PRE_K);      // 2 MB, L2-resident

// Ungapped scoring. Deliberately simple and NOT an attempt to reproduce mmseqs' bit scores: this
// stage only ranks candidates, and the winner is re-scored by mmseqs against the full scaffold.
constexpr int MATCH = 2;
constexpr int MISMATCH = -3;
constexpr int XDROP = 20;      // stop extending once the score falls this far below its peak

// Default significance floor, calibrated against mmseqs rather than chosen.
//
// mmseqs applies `-e 1e-3`; this scheme has no e-value, so without a floor every seeded diagonal
// with a positive extension is reported. Measured on 100,000 amplicon reads against the 924-target
// reference, reads with a hit per side:
//
//              V        J       C
//   no floor  53,600  47,043  43,010
//   >= 40     51,874  46,997     471
//   mmseqs    52,004  47,041     473
//
// 40 reproduces mmseqs to within 0.25 % on V and 0.09 % on J. The C column is the tell: half the
// spurious constant-region hits score exactly **38**, which is a bare 16-mer seed (16 x MATCH = 32)
// plus a couple of flanking matches -- a seed and nothing more. 40 is the first value above that.
constexpr int MIN_SCORE = 40;

// Indel detection. An ungapped extension follows ONE diagonal, so a read carrying an insertion or
// deletion relative to germline scores only up to the indel, and its two halves land on two
// diagonals of the SAME target offset by the indel length. That is a signature, not a nuisance --
// it is present in the vote list before any extension runs and costs one pass to read.
//
// Measured on 341,294 real IGH mates (two IGH RepSeq amplicons, IgBLAST truth, benchmark round 12):
// 3.18 % of reads carry a V indel, and the rate tracks SHM load, because AID produces indels and
// not only substitutions --
//
//   V identity   >=98     95-98    90-95    <90
//   indel reads  0.74 %   1.63 %   3.56 %   8.00 %
//
// On hypermutated IGH that is a 1-in-12 event, in exactly the population every fast path already
// handles worst. A flagged read is sent for GAPPED realignment, never dropped, so a false positive
// costs a little speed and cannot cost a read. That asymmetry is why the support threshold is set
// low rather than tuned tight.
constexpr uint32_t MIN_DIAG_SEEDS = 3;   // seeds a diagonal needs before it counts as real support

struct Hit {
    int32_t target;
    int32_t score;
    int32_t qstart, qend;      // 1-based inclusive, in ORIGINAL read orientation
    int32_t tstart, tend;      // 1-based inclusive, on the forward target
    int32_t rc;                // 1 if the read matched as its reverse complement
    int32_t split;             // 1 if this target carries the two-diagonal signature of an indel
};

std::string revcomp(const std::string& s) {
    std::string out(s.size(), 'N');
    for (size_t i = 0; i < s.size(); ++i) {
        const char c = s[s.size() - 1 - i];
        switch (c) {
            case 'A': case 'a': out[i] = 'T'; break;
            case 'C': case 'c': out[i] = 'G'; break;
            case 'G': case 'g': out[i] = 'C'; break;
            case 'T': case 't': out[i] = 'A'; break;
            default:            out[i] = 'N'; break;
        }
    }
    return out;
}

class SegmentMapper {
public:
    SegmentMapper(const std::vector<std::string>& seqs, const std::vector<int>& groups, int k)
        : k_(k), seqs_(seqs), groups_(groups) {
        if (k_ < PRE_K || k_ > 31) throw std::invalid_argument("k must be between 12 and 31");
        if (seqs_.size() != groups_.size())
            throw std::invalid_argument("seqs and groups must be the same length");
        if (seqs_.size() > 0xFFFFu) throw std::invalid_argument("at most 65535 targets");
        mask_ = (uint64_t(1) << (2 * k_)) - 1;
        shift_ = 2 * (k_ - PRE_K);
        n_groups_ = 0;
        for (int g : groups_) n_groups_ = std::max(n_groups_, g + 1);

        // (code, target, offset) triples, then collapse to CSR.
        std::vector<std::pair<uint64_t, uint32_t>> raw;
        for (size_t t = 0; t < seqs_.size(); ++t) {
            const std::string& s = seqs_[t];
            if (s.size() > 0xFFFFu) throw std::invalid_argument("target longer than 65535 nt");
            int run = 0;
            uint64_t code = 0;
            for (size_t i = 0; i < s.size(); ++i) {
                const int8_t b = BASE_IDX[static_cast<unsigned char>(s[i])];
                if (b < 0) { run = 0; code = 0; continue; }
                code = ((code << 2) | uint64_t(b)) & mask_;
                if (++run < k_) continue;
                const uint32_t off = uint32_t(i + 1 - k_);           // 0-based start of the window
                raw.emplace_back(code, (uint32_t(t) << 16) | off);
            }
        }
        std::sort(raw.begin(), raw.end());

        codes_.reserve(raw.size());
        starts_.reserve(raw.size() + 1);
        postings_.reserve(raw.size());
        for (size_t i = 0; i < raw.size();) {
            const uint64_t c = raw[i].first;
            codes_.push_back(c);
            starts_.push_back(uint32_t(postings_.size()));
            while (i < raw.size() && raw[i].first == c) postings_.push_back(raw[i++].second);
        }
        starts_.push_back(uint32_t(postings_.size()));
        codes_.shrink_to_fit(); starts_.shrink_to_fit(); postings_.shrink_to_fit();

        pre_.assign(PRE_BITS / 64, 0);
        for (uint64_t c : codes_) {
            const uint64_t p = c >> shift_;
            pre_[p >> 6] |= (uint64_t(1) << (p & 63));
        }
    }

    size_t size() const { return codes_.size(); }
    size_t postings() const { return postings_.size(); }
    int k() const { return k_; }
    int n_groups() const { return n_groups_; }

    // Best hit per (read, group), plus exactly-tied hits up to `max_tied` for groups that allow
    // them. Returns a flat list of tuples in `_SEGMENT_FORMAT` order plus the strand and indel
    // flags:
    //   (query_index, target_index, score, qstart, qend, tstart, rc, split)
    // `max_indel` = 0 disables indel detection and `split` is then always 0.
    py::list map(py::sequence queries, int max_tied, int min_score, int threads,
                 int max_indel) const {
        const size_t n = size_t(py::len(queries));
        py::list out;
        if (n == 0) return out;
        if (max_tied < 1) max_tied = 1;

        // Same pattern as Prefilter::filter -- take every pointer under the GIL, then drop it.
        // `PyUnicode_AsUTF8AndSize` hands back a view into the str's own buffer, and `queries`
        // keeps the objects alive for the call.
        std::vector<std::pair<const char*, size_t>> views;
        views.reserve(n);
        for (auto item : queries) {
            Py_ssize_t len = 0;
            const char* p = PyUnicode_AsUTF8AndSize(item.ptr(), &len);
            if (p == nullptr) throw py::error_already_set();
            views.emplace_back(p, size_t(len));
        }

        std::vector<std::vector<Hit>> per_read(n);
        size_t nthread = threads > 0 ? size_t(threads) : 1;
        nthread = std::min(nthread, n);
        auto worker = [&](size_t lo, size_t hi) {
            std::vector<uint64_t> votes;                 // (target<<32) | biased diagonal
            std::vector<Hit> hits;
            for (size_t i = lo; i < hi; ++i)
                per_read[i] = best_hits(views[i].first, views[i].second, max_tied, min_score,
                                        max_indel, votes, hits);
        };
        {
            py::gil_scoped_release rel;
            if (nthread <= 1) {
                worker(0, n);
            } else {
                std::vector<std::thread> pool;
                pool.reserve(nthread);
                const size_t step = (n + nthread - 1) / nthread;
                for (size_t t = 0; t < nthread; ++t) {
                    const size_t lo = t * step, hi = std::min(n, lo + step);
                    if (lo >= hi) break;
                    pool.emplace_back(worker, lo, hi);
                }
                for (auto& th : pool) th.join();
            }
        }
        for (size_t i = 0; i < n; ++i)
            for (const Hit& h : per_read[i])
                out.append(py::make_tuple(i, h.target, h.score, h.qstart, h.qend, h.tstart, h.rc,
                                          h.split));
        return out;
    }

private:
    // Seed, vote by diagonal, extend the best diagonals ungapped, keep the best per group.
    std::vector<Hit> best_hits(const char* s, size_t n, int max_tied, int min_score, int max_indel,
                               std::vector<uint64_t>& votes, std::vector<Hit>& hits) const {
        hits.clear();
        if (n < size_t(k_)) return {};
        const std::string fwd(s, n);
        const std::string rev = revcomp(fwd);
        for (int rc = 0; rc < 2; ++rc)
            collect(rc ? rev : fwd, rc, min_score, max_indel, votes, hits);

        // Best per group, then exact ties on score. Ties break on target index, which is the
        // reference's own order -- the same shape of rule `_segment_rows` applies, so the choice
        // is deterministic rather than dependent on scan order.
        std::sort(hits.begin(), hits.end(), [](const Hit& a, const Hit& b) {
            if (a.score != b.score) return a.score > b.score;
            return a.target < b.target;
        });
        std::vector<Hit> out;
        std::vector<int> kept(size_t(n_groups_), 0);
        std::vector<int> top(size_t(n_groups_), INT32_MIN);
        for (const Hit& h : hits) {
            const int g = groups_[size_t(h.target)];
            if (kept[size_t(g)] == 0) {
                top[size_t(g)] = h.score;
                kept[size_t(g)] = 1;
                out.push_back(h);
            } else if (h.score == top[size_t(g)] && kept[size_t(g)] < max_tied) {
                ++kept[size_t(g)];
                out.push_back(h);
            }
        }
        return out;
    }

    // ⛔ `min_score` is a PARAMETER, not a member. `map` is const and every worker thread shares
    // one SegmentMapper, so stashing per-call state on the object would be a data race that only
    // shows up under threads -- the class of bug this project has spent a lot of time removing.
    void collect(const std::string& q, int rc, int min_score, int max_indel,
                 std::vector<uint64_t>& votes, std::vector<Hit>& hits) const {
        votes.clear();
        const size_t n = q.size();
        int run = 0;
        uint64_t code = 0;
        for (size_t i = 0; i < n; ++i) {
            const int8_t b = BASE_IDX[static_cast<unsigned char>(q[i])];
            if (b < 0) { run = 0; code = 0; continue; }
            code = ((code << 2) | uint64_t(b)) & mask_;
            if (++run < k_) continue;
            const uint64_t p = code >> shift_;
            if (!((pre_[p >> 6] >> (p & 63)) & 1)) continue;      // level 1: one cached access
            const auto it = std::lower_bound(codes_.begin(), codes_.end(), code);
            if (it == codes_.end() || *it != code) continue;
            const size_t ci = size_t(it - codes_.begin());
            const uint32_t qpos = uint32_t(i + 1 - size_t(k_));   // 0-based window start in q
            for (uint32_t pi = starts_[ci]; pi < starts_[ci + 1]; ++pi) {
                const uint32_t post = postings_[pi];
                const uint32_t t = post >> 16, tpos = post & 0xFFFFu;
                // Bias the diagonal so it is non-negative and packs into the low half.
                const int64_t diag = int64_t(qpos) - int64_t(tpos) + 0x8000;
                votes.push_back((uint64_t(t) << 32) | uint64_t(diag & 0xFFFFFFFF));
            }
        }
        if (votes.empty()) return;
        std::sort(votes.begin(), votes.end());
        // Votes sort by (target << 32 | biased diagonal), so every diagonal of a target is
        // contiguous and ascending. One walk per target-run therefore yields both the per-diagonal
        // seed support and the two-diagonal indel signature, without a second data structure.
        for (size_t s = 0; s < votes.size();) {
            const uint32_t t = uint32_t(votes[s] >> 32);
            size_t e = s;
            while (e < votes.size() && uint32_t(votes[e] >> 32) == t) ++e;

            // Does this target carry two well-supported diagonals a plausible indel apart?
            int32_t split = 0;
            if (max_indel > 0) {
                int32_t prev = 0;
                bool have_prev = false;
                for (size_t i = s; i < e;) {
                    size_t j = i;
                    while (j < e && votes[j] == votes[i]) ++j;
                    if (uint32_t(j - i) >= MIN_DIAG_SEEDS) {
                        const int32_t d = int32_t(uint32_t(votes[i] & 0xFFFFFFFF)) - 0x8000;
                        // Ascending, so the gap is positive. Bound it: two diagonals far apart are
                        // a repeat or a spurious seed, not one indel, and routing those to gapped
                        // realignment would buy nothing at real cost.
                        if (have_prev && d - prev <= max_indel) { split = 1; break; }
                        prev = d;
                        have_prev = true;
                    }
                    i = j;
                }
            }

            // One extension per distinct (target, diagonal). Reads are short, so this list is small.
            for (size_t i = s; i < e;) {
                size_t j = i;
                while (j < e && votes[j] == votes[i]) ++j;
                const int32_t diag = int32_t(uint32_t(votes[i] & 0xFFFFFFFF)) - 0x8000;
                Hit h;
                if (extend(q, rc, t, diag, h) && h.score >= min_score) {
                    h.split = split;
                    hits.push_back(h);
                }
                i = j;
            }
            s = e;
        }
    }

    // Ungapped X-drop extension from the whole overlap the diagonal implies.
    bool extend(const std::string& q, int rc, uint32_t t, int32_t diag, Hit& out) const {
        const std::string& tgt = seqs_[t];
        const int64_t qn = int64_t(q.size()), tn = int64_t(tgt.size());
        // q[i] pairs with tgt[i - diag]; the overlap is where both are in range.
        const int64_t lo = std::max<int64_t>(0, diag);
        const int64_t hi = std::min<int64_t>(qn, tn + diag);
        if (hi - lo < k_) return false;

        int best = 0, cur = 0;
        int64_t bs = lo, be = lo - 1, run_start = lo;
        for (int64_t i = lo; i < hi; ++i) {
            const char a = q[size_t(i)], b = tgt[size_t(i - diag)];
            const bool ok = (a == b) && BASE_IDX[static_cast<unsigned char>(a)] >= 0;
            cur += ok ? MATCH : MISMATCH;
            if (cur > best) { best = cur; bs = run_start; be = i; }
            if (cur < 0) { cur = 0; run_start = i + 1; }
            else if (best - cur > XDROP) break;
        }
        if (be < bs) return false;
        out.target = int32_t(t);
        out.score = best;
        out.tstart = int32_t(bs - diag) + 1;             // 1-based on the forward target
        out.tend = int32_t(be - diag) + 1;
        out.rc = rc;
        if (rc == 0) {
            out.qstart = int32_t(bs) + 1;
            out.qend = int32_t(be) + 1;
        } else {
            // Back to the ORIGINAL read's coordinates. mmseqs signals a minus-strand hit with
            // qstart > qend and arda's `_align_implied` relies on that, so the convention is kept.
            out.qstart = int32_t(qn - bs);
            out.qend = int32_t(qn - be);
        }
        return true;
    }

    int k_, shift_, n_groups_;
    uint64_t mask_;
    std::vector<std::string> seqs_;
    std::vector<int> groups_;
    std::vector<uint64_t> codes_;
    std::vector<uint32_t> starts_;
    std::vector<uint32_t> postings_;
    std::vector<uint64_t> pre_;
};

}  // namespace

PYBIND11_MODULE(_segmap, m) {
    m.doc() = "Structure-aware segment mapper: best V and best J per read, without a homology search.";
    py::class_<SegmentMapper>(m, "SegmentMapper")
        .def(py::init<const std::vector<std::string>&, const std::vector<int>&, int>(),
             py::arg("sequences"), py::arg("groups"), py::arg("k") = 16,
             "Index `sequences` (forward strand only). `groups` assigns each target a side, and "
             "the mapper returns the best hit per (read, group).")
        .def("map", &SegmentMapper::map,
             py::arg("queries"), py::arg("max_tied") = 8, py::arg("min_score") = MIN_SCORE,
             py::arg("threads") = 1, py::arg("max_indel") = 0,
             "(query_index, target_index, score, qstart, qend, tstart, rc, split) per best/tied "
             "hit. `max_indel` > 0 flags targets carrying two well-supported diagonals that far "
             "apart -- the signature of an indel, which one ungapped extension cannot score.")
        .def_property_readonly("size", &SegmentMapper::size, "Distinct indexed k-mers.")
        .def_property_readonly("postings", &SegmentMapper::postings, "Total (k-mer, position) entries.")
        .def_property_readonly("k", &SegmentMapper::k, "k-mer length.")
        .def_property_readonly("n_groups", &SegmentMapper::n_groups, "Number of target groups.");
}
