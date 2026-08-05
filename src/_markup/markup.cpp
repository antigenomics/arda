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

#include <algorithm>
#include <array>
#include <cctype>
#include <climits>
#include <map>
#include <stdexcept>
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
// Returns (score, i_start, i_end, d_start, d_end): best score, the 0-based
// inclusive offsets of the matched segment within `interior`, and the 0-based
// inclusive offsets of the aligned span within the D germline `d`. Because the
// alignment is gapless, the two spans have equal length and lie on one diagonal.
// (0, -1, -1, -1, -1) if no positive-scoring segment exists. Comparison is
// case-insensitive; N (or any non-matching base) counts as a mismatch.
static std::tuple<int, int, int, int, int> d_local_align(const std::string &interior,
                                                          const std::string &d) {
    const int n = static_cast<int>(interior.size());
    const int m = static_cast<int>(d.size());
    int best = 0, bs = -1, be = -1, best_off = 0;
    auto up = [](char c) -> char {
        return (c >= 'a' && c <= 'z') ? static_cast<char>(c - 32) : c;
    };
    // Each diagonal is a fixed offset = i - j (i over interior, j over D), so the
    // D-side coordinate of any interior position i on that diagonal is j = i - off.
    for (int off = -(m - 1); off <= n - 1; ++off) {
        int cur = 0, cur_start = -1;
        const int i_lo = off > 0 ? off : 0;
        const int i_hi = (n - 1 < off + m - 1) ? n - 1 : off + m - 1;
        for (int i = i_lo; i <= i_hi; ++i) {
            const int j = i - off;
            const int sc = (up(interior[i]) == up(d[j])) ? 1 : -1;
            if (cur <= 0) { cur = sc; cur_start = i; }
            else { cur += sc; }
            if (cur > best) { best = cur; bs = cur_start; be = i; best_off = off; }
        }
    }
    const int ds = (bs < 0) ? -1 : bs - best_off;
    const int de = (be < 0) ? -1 : be - best_off;
    return {best, bs, be, ds, de};
}

// Stitch a contig's alignment to the scaffold from its constituent reads'
// alignments -- the alternative to re-aligning the assembled contig with mmseqs.
//
// A contig is a consensus of reads that already carry per-read alignments to the
// same scaffold (Stage 1 emits qaln/taln/qstart/tstart per read). Given the
// assembly layout (each read's 0-based offset within the contig, in contig
// orientation), the contig->scaffold alignment is the per-scaffold-column
// consensus of the reads' alignments. This is the O(reads x columns) reduction
// worth doing in C++ when a scRNA-seq sample has ~10^5 contigs.
//
// Inputs (one contig): the reads' aligned strings qalns[i]/talns[i] (coding
// strand, '-' for gaps), 1-based read start qstarts[i] and scaffold start
// tstarts[i], 0-based contig offset offsets[i] of the read's first base, and the
// contig nucleotide string (the authoritative consensus base at each position).
//
// Returns (qaln, taln, qstart, qend, tstart, tend): the contig->scaffold gapped
// alignment with contig on the query side, 1-based contig qstart/qend and 1-based
// scaffold tstart/tend. Downstream (segment_cigars / region transfer) is identical
// to a freshly aligned contig, so both paths share transfer_hit.
//
// Ceiling (ponytail: first-seen read wins a column; assumes a contiguous ungapped
// consensus layout -- the normal output of reference-anchored assembly). The first
// read to cover a scaffold column fixes its base and any insertion preceding it;
// later reads do not override it (no majority vote). A scaffold column inside
// [tstart,tend] that no read covers is a coverage gap -- not representable as one
// contig -- and throws. Upgrade path: majority vote per column + explicit gap
// handling if real assemblies produce disagreeing overlaps.
static std::tuple<std::string, std::string, int, int, int, int> merge_alignment(
    const std::vector<std::string> &qalns, const std::vector<std::string> &talns,
    const std::vector<int> &qstarts, const std::vector<int> &tstarts,
    const std::vector<int> &offsets, const std::string &contig) {

    const size_t R = qalns.size();
    if (R == 0)
        throw std::invalid_argument("merge_alignment: no reads");

    // Scaffold column span [t_lo, t_hi] (1-based) = union of the reads' aligned
    // spans. tend_i = tstart_i + (non-gap chars in talns[i]) - 1.
    int t_lo = INT_MAX, t_hi = INT_MIN;
    for (size_t i = 0; i < R; ++i) {
        int tref = 0;
        for (char c : talns[i]) if (c != '-') ++tref;
        const int te = tstarts[i] + tref - 1;
        if (tstarts[i] < t_lo) t_lo = tstarts[i];
        if (te > t_hi) t_hi = te;
    }
    const int W = t_hi - t_lo + 1;
    std::vector<char> qbase(W, 0);   // consensus contig base per column; '-' = deletion; 0 = uncovered
    std::vector<char> tbase(W, 0);   // scaffold base per column
    std::vector<int> qpos_at(W, 0);  // 1-based contig position of the M at this column
    std::map<int, std::string> ins;  // inserted contig bases keyed by the scaffold column they precede

    for (size_t i = 0; i < R; ++i) {
        const std::string &qa = qalns[i], &ta = talns[i];
        const size_t n = std::min(qa.size(), ta.size());
        int qpos = qstarts[i], tpos = tstarts[i];
        std::string pending;  // inserted bases since the last aligned column, in this read
        for (size_t k = 0; k < n; ++k) {
            const bool qg = (qa[k] == '-'), tg = (ta[k] == '-');
            if (!qg && !tg) {                        // M
                const int cp = offsets[i] + qpos;    // 1-based contig position
                const int idx = tpos - t_lo;
                if (qbase[idx] == 0) {
                    if (cp < 1 || cp > (int)contig.size())
                        throw std::invalid_argument("merge_alignment: read offset outside contig");
                    qbase[idx] = contig[cp - 1];
                    tbase[idx] = ta[k];
                    qpos_at[idx] = cp;
                    if (!pending.empty()) ins[tpos] = pending;  // insertion precedes this column
                }
                pending.clear();
                ++qpos; ++tpos;
            } else if (!qg && tg) {                  // I: contig base absent from scaffold
                const int cp = offsets[i] + qpos;
                if (cp >= 1 && cp <= (int)contig.size()) pending.push_back(contig[cp - 1]);
                ++qpos;
            } else if (qg && !tg) {                  // D: scaffold base absent from contig
                const int idx = tpos - t_lo;
                if (qbase[idx] == 0) { qbase[idx] = '-'; tbase[idx] = ta[k]; }
                pending.clear();
                ++tpos;
            }
        }
    }

    std::string out_q, out_t;
    out_q.reserve(W + 8); out_t.reserve(W + 8);
    int qstart = -1, qend = -1;
    for (int t = t_lo; t <= t_hi; ++t) {
        auto it = ins.find(t);
        if (it != ins.end()) { out_q += it->second; out_t.append(it->second.size(), '-'); }
        const int idx = t - t_lo;
        if (qbase[idx] == 0)
            throw std::invalid_argument("merge_alignment: uncovered scaffold column (coverage gap)");
        out_q.push_back(qbase[idx]);
        out_t.push_back(tbase[idx]);
        if (qbase[idx] != '-') {                     // an M: it carries a contig position
            if (qstart < 0) qstart = qpos_at[idx];
            qend = qpos_at[idx];
        }
    }
    return {out_q, out_t, qstart, qend, t_lo, t_hi};
}


// ---------------------------------------------------------------------------
// Per-segment AIRR CIGARs.
//
// The pure-Python `arda.annotate.cigar.segment_cigars` walks the alignment column by column and
// makes two function calls per column (`_classify`, `_germline_pos`). Profiled over the real-read
// fixture that was 47,309 + 46,943 calls for 524 mapped reads, and 63 % of `transfer_hit`'s
// 130 us -- which is the term that caps arda on receptor-rich libraries, where the prefilter has
// nothing left to remove and every surviving read builds a record.
//
// Straight port, column for column, so the two implementations can be diffed on real data.
// V is [1, t_vend], J is [t_jstart, t_vjend], C is > t_vjend; the N-pad between V and J belongs to
// no germline (it is the np region). A zero boundary means the segment is absent -- a `J + C`
// scaffold has t_vend == 0.
// ---------------------------------------------------------------------------

static inline int seg_key(int tpos, int t_vend, int t_jstart, int t_vjend) {
    if (t_vend && tpos <= t_vend) return 0;                              // v
    if (t_jstart && t_vjend && tpos >= t_jstart && tpos <= t_vjend) return 1;  // j
    if (t_vjend && tpos > t_vjend) return 2;                             // c
    return -1;
}

static inline int germline_pos(int key, int t, int t_jstart, int t_vjend) {
    if (key == 0) return t;                       // V: scaffold position IS the germline position
    if (key == 1) return t - t_jstart + 1;
    return t - t_vjend;                           // C: germline starts one past the V-J end
}

static std::string build_cigar_cpp(int q_lead, int g_lead, const std::string& ops, int q_trail) {
    std::string body;
    for (size_t i = 0; i < ops.size();) {
        size_t j = i;
        while (j < ops.size() && ops[j] == ops[i]) ++j;
        body += std::to_string(j - i);
        body += ops[i];
        i = j;
    }
    std::string out;
    if (q_lead > 0) out += std::to_string(q_lead) + "S";
    if (g_lead > 0) out += std::to_string(g_lead) + "N";
    out += body;
    // Trailing germline N is intentionally omitted (it is optional in AIRR).
    if (q_trail > 0) out += std::to_string(q_trail) + "S";
    return out;
}

std::map<std::string, std::string> segment_cigars(
        const std::string& qaln, const std::string& taln,
        int qstart, int tstart, int qlen,
        int t_vend, int t_jstart, int t_vjend) {
    static const char* NAMES[3] = {"v_cigar", "j_cigar", "c_cigar"};
    std::string ops[3];
    int q_first[3] = {0, 0, 0}, q_last[3] = {0, 0, 0}, g_first[3] = {0, 0, 0};
    bool seen[3] = {false, false, false};

    int q = qstart, t = tstart;
    const size_t n = std::min(qaln.size(), taln.size());
    for (size_t i = 0; i < n; ++i) {
        const bool cq = qaln[i] != '-';
        const bool ct = taln[i] != '-';
        const char op = (cq && ct) ? 'M' : (cq ? 'I' : 'D');
        const int key = seg_key(t, t_vend, t_jstart, t_vjend);  // on an insertion, t is the next base
        if (key >= 0) {
            ops[key] += op;
            if (cq) {
                if (!seen[key]) { q_first[key] = q; seen[key] = true; }
                q_last[key] = q;
            }
            if (ct && g_first[key] == 0) g_first[key] = germline_pos(key, t, t_jstart, t_vjend);
        }
        if (cq) ++q;
        if (ct) ++t;
    }

    std::map<std::string, std::string> out;
    for (int k = 0; k < 3; ++k) {
        if (!seen[k]) continue;                   // no query base aligned to this segment
        const int g_lead = (g_first[k] ? g_first[k] : 1) - 1;
        std::string cig = build_cigar_cpp(q_first[k] - 1, g_lead, ops[k], qlen - q_last[k]);
        if (!cig.empty()) out[NAMES[k]] = cig;
    }
    return out;
}


// Fractional identity over germline positions in target range [t_lo, t_hi].
//
// Second per-column loop in `transfer_hit`, and once `segment_cigars` moved to C++ it was 44 % of
// what remained. Same walk: count each target-consuming column that falls in range, and how many
// of those are an exact base match. Returns -1.0 when nothing is covered, which the caller maps to
// the empty AIRR field -- 0.0 would be a real identity of zero and is a different statement.
double aln_identity(const std::string& qaln, const std::string& taln,
                    int tstart, int t_lo, int t_hi) {
    int t = tstart, covered = 0, ident = 0;
    const size_t n = std::min(qaln.size(), taln.size());
    for (size_t i = 0; i < n; ++i) {
        const char ta = taln[i];
        if (ta == '-') continue;                     // query-only column: no germline position
        if (t >= t_lo && t <= t_hi) {
            ++covered;
            const char qa = qaln[i];
            if (qa != '-' && std::toupper((unsigned char)qa) == std::toupper((unsigned char)ta))
                ++ident;
        }
        ++t;
    }
    return covered ? double(ident) / double(covered) : -1.0;
}

PYBIND11_MODULE(_markup, m) {
    m.doc() = "arda markup-transfer hot path (C++/pybind11)";
    m.attr("__version__") = "0.4.0";
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
    m.def("segment_cigars", &segment_cigars,
          py::arg("qaln"), py::arg("taln"), py::arg("qstart"), py::arg("tstart"),
          py::arg("qlen"), py::arg("t_vend"), py::arg("t_jstart"), py::arg("t_vjend"),
          "Per-segment AIRR CIGARs {v_cigar, j_cigar, c_cigar} for the segments with a body. "
          "Boundaries are 1-based scaffold positions; 0 means the segment is absent.");
    m.def("aln_identity", &aln_identity,
          py::arg("qaln"), py::arg("taln"), py::arg("tstart"), py::arg("t_lo"), py::arg("t_hi"),
          "Fractional identity over germline positions in target range [t_lo, t_hi]; "
          "-1.0 when no germline position is covered.");
    m.def("d_local_align", &d_local_align, py::arg("interior"), py::arg("d"),
          "Gapless local alignment (match=+1, mismatch=-1) of a short D germline "
          "against a query interior. Returns (score, i_start, i_end, d_start, d_end) "
          "with 0-based inclusive offsets of the best segment in `interior` and in "
          "the D germline; (0,-1,-1,-1,-1) if none.");
    m.def("merge_alignment", &merge_alignment,
          py::arg("qalns"), py::arg("talns"), py::arg("qstarts"), py::arg("tstarts"),
          py::arg("offsets"), py::arg("contig"),
          "Stitch a contig's scaffold alignment from its reads' alignments (the "
          "alternative to re-aligning the contig). Per-read qaln/taln, 1-based read "
          "qstart and scaffold tstart, 0-based contig offset, plus the contig string. "
          "Returns (qaln, taln, qstart, qend, tstart, tend). First read wins a "
          "column; a coverage gap or out-of-range offset throws.");
}
