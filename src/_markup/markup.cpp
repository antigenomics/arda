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

#include <array>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace py = pybind11;

using Interval = std::pair<int, int>;

// ---------------------------------------------------------------------------
// Fast sequence primitives (translation / back-translation / reverse-complement).
//
// These are deliberately API-compatible with mirpy's mir.basic.mirseq so that
// mirpy can later `import arda` and reuse them. Standard genetic code; codons
// with any non-ACGT base translate to 'X'; stop codons are '*'; a trailing
// partial codon is dropped (arda) — translate_bidi mirrors mirpy's '_' padding.
// ---------------------------------------------------------------------------

// 2-bit base index: A=0, C=1, G=2, T=3; anything else = -1.
static std::array<int8_t, 256> make_base_idx() {
    std::array<int8_t, 256> t{};
    t.fill(-1);
    t['A'] = t['a'] = 0; t['C'] = t['c'] = 1; t['G'] = t['g'] = 2; t['T'] = t['t'] = 3;
    return t;
}
static const std::array<int8_t, 256> BASE_IDX = make_base_idx();

// Codon table indexed by (b0*16 + b1*4 + b2) with A=0,C=1,G=2,T=3.
static const char AA_TABLE[64] = {
    'K','N','K','N', 'T','T','T','T', 'R','S','R','S', 'I','I','M','I',  // A__
    'Q','H','Q','H', 'P','P','P','P', 'R','R','R','R', 'L','L','L','L',  // C__
    'E','D','E','D', 'A','A','A','A', 'G','G','G','G', 'V','V','V','V',  // G__
    '*','Y','*','Y', 'S','S','S','S', '*','C','W','C', 'L','F','L','F',  // T__
};

// Complement LUT.
static std::array<char, 256> make_comp() {
    std::array<char, 256> t{};
    for (int i = 0; i < 256; ++i) t[i] = 'N';
    t['A'] = t['a'] = 'T'; t['T'] = t['t'] = 'A';
    t['G'] = t['g'] = 'C'; t['C'] = t['c'] = 'G'; t['N'] = t['n'] = 'N';
    return t;
}
static const std::array<char, 256> COMP = make_comp();

// Human (Kazusa) most-frequent codon per amino acid, for mock back-translation.
static char BT_TABLE[128][4];
static bool init_bt() {
    for (auto &c : BT_TABLE) { c[0] = 'N'; c[1] = 'N'; c[2] = 'N'; c[3] = '\0'; }
    const char *aa = "ARNDCQEGHILKMFPSTWYV";
    const char *co[] = {"GCC","AGG","AAC","GAC","TGC","CAG","GAG","GGC","CAC","ATC",
                        "CTG","AAG","ATG","TTC","CCC","AGC","ACC","TGG","TAC","GTG"};
    for (int i = 0; aa[i]; ++i) {
        unsigned u = static_cast<unsigned char>(aa[i]);
        BT_TABLE[u][0] = co[i][0]; BT_TABLE[u][1] = co[i][1]; BT_TABLE[u][2] = co[i][2];
    }
    return true;
}
static const bool BT_INIT = init_bt();

static std::string translate(const std::string &nt, int frame) {
    std::string out;
    const int n = static_cast<int>(nt.size());
    if (frame < 0) frame = 0;
    out.reserve((n - frame) / 3 + 1);
    for (int i = frame; i + 3 <= n; i += 3) {
        const int8_t a = BASE_IDX[(unsigned char)nt[i]];
        const int8_t b = BASE_IDX[(unsigned char)nt[i + 1]];
        const int8_t c = BASE_IDX[(unsigned char)nt[i + 2]];
        out.push_back((a < 0 || b < 0 || c < 0) ? 'X' : AA_TABLE[a * 16 + b * 4 + c]);
    }
    return out;
}

static int detect_coding_frame(const std::string &nt) {
    int best_frame = 0, best_stops = -1;
    for (int f = 0; f < 3; ++f) {
        int stops = 0;
        const int n = static_cast<int>(nt.size());
        for (int i = f; i + 3 <= n; i += 3) {
            const int8_t a = BASE_IDX[(unsigned char)nt[i]];
            const int8_t b = BASE_IDX[(unsigned char)nt[i + 1]];
            const int8_t c = BASE_IDX[(unsigned char)nt[i + 2]];
            if (a >= 0 && b >= 0 && c >= 0 && AA_TABLE[a * 16 + b * 4 + c] == '*') ++stops;
        }
        if (best_stops < 0 || stops < best_stops) {
            best_stops = stops; best_frame = f;
            if (stops == 0) break;
        }
    }
    return best_frame;
}

static std::string reverse_complement(const std::string &nt) {
    std::string out(nt.size(), 'N');
    for (size_t i = 0; i < nt.size(); ++i)
        out[nt.size() - 1 - i] = COMP[(unsigned char)nt[i]];
    return out;
}

static std::string back_translate(const std::string &aa, const std::string &unknown) {
    std::string out;
    out.reserve(aa.size() * 3);
    for (char ch : aa) {
        unsigned u = static_cast<unsigned char>(ch);
        if (u < 128 && BT_TABLE[u][0] != 'N')
            out.append(BT_TABLE[u], 3);
        else
            out.append(unknown);
    }
    return out;
}

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

// Gapless local alignment of a (short) D germline against a query interior.
//
// D segments are tiny (~8-31 nt) and exonuclease-trimmed on both ends, so the
// useful signal is the best-scoring contiguous ungapped match between any
// substring of the interior and any substring of D. mmseqs' k-mer prefilter is
// unreliable at this length, so we brute-force every diagonal (Kadane's
// maximum-subarray per diagonal) with match=+1 / mismatch=-1 scoring. The D set
// per locus is small (≤ ~40 alleles) and the interior is short, so this is cheap.
//
// Returns (score, start, end): best score and the 0-based inclusive offsets of
// the matched segment within `interior`. (0, -1, -1) if no positive-scoring
// segment exists. Comparison is case-insensitive; N (or any non-matching base)
// counts as a mismatch.
static std::tuple<int, int, int> d_local_align(const std::string &interior,
                                                const std::string &d) {
    const int n = static_cast<int>(interior.size());
    const int m = static_cast<int>(d.size());
    int best = 0, bs = -1, be = -1;
    auto up = [](char c) -> char {
        return (c >= 'a' && c <= 'z') ? static_cast<char>(c - 32) : c;
    };
    // Each diagonal is a fixed offset = i - j (i over interior, j over D).
    for (int off = -(m - 1); off <= n - 1; ++off) {
        int cur = 0, cur_start = -1;
        const int i_lo = off > 0 ? off : 0;
        const int i_hi = (n - 1 < off + m - 1) ? n - 1 : off + m - 1;
        for (int i = i_lo; i <= i_hi; ++i) {
            const int j = i - off;
            const int sc = (up(interior[i]) == up(d[j])) ? 1 : -1;
            if (cur <= 0) { cur = sc; cur_start = i; }
            else { cur += sc; }
            if (cur > best) { best = cur; bs = cur_start; be = i; }
        }
    }
    return {best, bs, be};
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

    // Fast sequence primitives (also consumable by mirpy).
    m.def("translate", &translate, py::arg("nt"), py::arg("frame") = 0,
          "Translate a nucleotide string from `frame` (0/1/2). Non-ACGT codons "
          "-> 'X', stops -> '*', trailing partial codon dropped.");
    m.def("detect_coding_frame", &detect_coding_frame, py::arg("nt"),
          "Return the reading frame (0/1/2) with the fewest stop codons.");
    m.def("reverse_complement", &reverse_complement, py::arg("nt"),
          "Reverse-complement a nucleotide string (non-ACGT -> 'N').");
    m.def("back_translate", &back_translate, py::arg("aa"), py::arg("unknown") = "NNN",
          "Mock back-translation using the most-frequent human (Kazusa) codon per "
          "amino acid; unknown residues -> `unknown` (default 'NNN').");
    m.def("d_local_align", &d_local_align, py::arg("interior"), py::arg("d"),
          "Gapless local alignment (match=+1, mismatch=-1) of a short D germline "
          "against a query interior. Returns (score, start, end) with 0-based "
          "inclusive offsets of the best segment in `interior`; (0,-1,-1) if none.");
}
