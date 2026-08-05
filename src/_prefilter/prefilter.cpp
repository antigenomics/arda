// arda k-mer prefilter: reject reads that cannot possibly align, before MMseqs2 sees them.
//
// Why this exists
// ---------------
// On bulk RNA-seq the receptor fraction is 0.02-3 %, so `mmseqs search` spends essentially all
// of its time proving that reads are not receptor reads. Measured on 4 M reads of SRR10611239
// (947 map, 0.024 %): the search alone is 48.9 s. Every read costs the same whether it hits or
// not, and the fitted cost model says so directly -- `wall ~ reads/46,353 + hits/350`, i.e. the
// scan term is set by the READ COUNT, not by the answer.
//
// A read can only align to a V(D)J scaffold if it shares an exact 16-mer with one, up to the
// mismatches SHM introduces. Testing that is a hash lookup; proving it with Smith-Waterman is
// not. This module does the cheap test so MMseqs2 only sees the survivors.
//
// The design is OPTIMIZATION.md section 3-4, and its measurements pick every constant here:
//
//   * **k = 16, exact, both strands, >= 1 hit.** k=12 passes 62-65 % of reads (useless); k>=18
//     adds nothing over 16. Requiring >=2 or >=3 hits barely moves the pass rate (4.64 % ->
//     3.91 %) while FN climbs -- a real read sharing one exact 16-mer usually shares many, so the
//     extra hits buy nothing. `min_hits` is still a parameter because SHM-heavy libraries are the
//     one place that might not hold, and it must be measurable without a rebuild.
//   * **Index the constant region, not just V+J.** This is the finding that decides the whole
//     component: against a V+J-only reference, FN is 16.29 % and **69.27 % of what it loses are
//     J->C reads**. Indexing C takes FN to 0.53 %. The caller passes the target sequences, so
//     this is the caller's job to get right -- see arda/prefilter.py, which indexes the same
//     `alleles.fasta` MMseqs2 searches, so the two cannot drift apart.
//   * **No Bloom filter.** 171,918 seeds is 1.4 MB as a sorted uint64 array. A Bloom filter
//     trades exactness for memory that is not scarce here.
//
// Two-level lookup: a 4^12 bitset (2 MB, L2-resident) rejects most windows in one memory access;
// only survivors pay the binary search into the sorted seed array.
//
// Sequences containing N: windows covering an N are dropped, never guessed. One substitution
// destroys k consecutive windows, which is why the read must be scanned at every offset.

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

// 2-bit base index: A=0, C=1, G=2, T=3; anything else (N, IUPAC codes) = -1 and breaks the run.
std::array<int8_t, 256> make_base_idx() {
    std::array<int8_t, 256> t{};
    t.fill(-1);
    t['A'] = t['a'] = 0; t['C'] = t['c'] = 1; t['G'] = t['g'] = 2; t['T'] = t['t'] = 3;
    return t;
}
const std::array<int8_t, 256> BASE_IDX = make_base_idx();

// Bits of the k-mer code used by the first-level bitset. 12 bases -> 4^12 bits -> 2 MB.
constexpr int PRE_K = 12;
constexpr size_t PRE_BITS = size_t(1) << (2 * PRE_K);

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

class Prefilter {
public:
    Prefilter(const std::vector<std::string>& targets, int k) : k_(k) {
        if (k_ < PRE_K || k_ > 32)
            throw std::invalid_argument("k must be between 12 and 32");
        mask_ = (k_ == 32) ? ~uint64_t(0) : ((uint64_t(1) << (2 * k_)) - 1);
        shift_ = 2 * (k_ - PRE_K);

        // Index BOTH strands of the reference rather than canonicalising k-mers. Canonical form
        // would halve the index but forces the same computation on every query window; the index
        // is 1-3 MB either way, and the query side is where 4 M reads x ~85 windows are spent.
        for (const auto& t : targets) {
            add_all(t);
            add_all(revcomp(t));
        }
        std::sort(seeds_.begin(), seeds_.end());
        seeds_.erase(std::unique(seeds_.begin(), seeds_.end()), seeds_.end());
        seeds_.shrink_to_fit();

        pre_.assign(PRE_BITS / 64, 0);
        for (uint64_t s : seeds_) {
            const uint64_t p = s >> shift_;
            pre_[p >> 6] |= (uint64_t(1) << (p & 63));
        }
    }

    size_t size() const { return seeds_.size(); }
    int k() const { return k_; }

    // Number of distinct query windows found in the index, capped: the caller only ever compares
    // against `min_hits`, so counting past it is wasted work on the reads that pass.
    int hits(const std::string& s, int min_hits) const {
        return scan(s.data(), s.size(), min_hits);
    }

    int scan(const char* s, size_t n, int min_hits) const {
        int run = 0, found = 0;
        uint64_t code = 0;
        for (size_t i = 0; i < n; ++i) {
            const int8_t b = BASE_IDX[static_cast<unsigned char>(s[i])];
            if (b < 0) { run = 0; code = 0; continue; }  // N: drop every window covering it
            code = ((code << 2) | uint64_t(b)) & mask_;
            if (++run < k_) continue;
            const uint64_t p = code >> shift_;
            if (!((pre_[p >> 6] >> (p & 63)) & 1)) continue;   // level 1: one cache-resident access
            if (std::binary_search(seeds_.begin(), seeds_.end(), code)) {
                if (++found >= min_hits) return found;         // early exit: the common pass case
            }
        }
        return found;
    }

    // Per-sequence pass/fail. Returns uint8 rather than bool because std::vector<bool> is a
    // bitfield and cannot be written from several threads without a race.
    //
    // Takes a py::sequence, NOT a std::vector<std::string>. pybind11's automatic conversion would
    // allocate and copy every read before a single thread starts -- 4 M allocations under the GIL,
    // measured at 2.34 s on 32 cluster threads against 5.39 M reads/s on ONE laptop core, i.e. the
    // copy was several times the scan it was feeding. `PyUnicode_AsUTF8AndSize` hands back a
    // pointer into the str object's own buffer instead; the list keeps every object alive for the
    // duration of the call, so the views stay valid after the GIL is dropped.
    std::vector<uint8_t> mask(py::sequence queries, int min_hits, int threads) const {
        const size_t n = size_t(py::len(queries));
        std::vector<uint8_t> out(n, 0);
        if (n == 0) return out;
        if (min_hits < 1) min_hits = 1;

        std::vector<std::pair<const char*, size_t>> views;
        views.reserve(n);
        for (auto item : queries) {
            Py_ssize_t len = 0;
            const char* p = PyUnicode_AsUTF8AndSize(item.ptr(), &len);
            if (p == nullptr) throw py::error_already_set();
            views.emplace_back(p, size_t(len));
        }

        size_t nthread = threads > 0 ? size_t(threads) : 1;
        nthread = std::min(nthread, n);

        auto worker = [&](size_t lo, size_t hi) {
            for (size_t i = lo; i < hi; ++i)
                out[i] = scan(views[i].first, views[i].second, min_hits) >= min_hits ? 1 : 0;
        };
        // The GIL is dropped only now, once every pointer has been taken -- `PyUnicode_AsUTF8AndSize`
        // is a CPython call and must not run without it. `views` and `out` are then touched through
        // disjoint index ranges, so the threads need no synchronisation.
        py::gil_scoped_release rel;
        if (nthread <= 1) {
            worker(0, n);
            return out;
        }
        std::vector<std::thread> pool;
        pool.reserve(nthread);
        const size_t step = (n + nthread - 1) / nthread;
        for (size_t t = 0; t < nthread; ++t) {
            const size_t lo = t * step, hi = std::min(n, lo + step);
            if (lo >= hi) break;
            pool.emplace_back(worker, lo, hi);
        }
        for (auto& th : pool) th.join();
        return out;
    }

    // Filter `records` (a sequence of (id, sequence) tuples) down to the survivors, in one call.
    //
    // `mask` alone is not enough: around it the caller had to build a 4 M-element list of
    // sequences to pass in, and unpack a 4 M-element list of ints coming back -- two full Python
    // passes over every read, which is why the 2.66x the scan gained showed up as 1.16x end to
    // end. Here the only Python objects created are the survivors, and on bulk that is 0.5-2 % of
    // the input. The returned tuples are the caller's own, not copies.
    py::list filter(py::sequence records, int seq_index, int min_hits, int threads) const {
        const size_t n = size_t(py::len(records));
        py::list out;
        if (n == 0) return out;
        if (min_hits < 1) min_hits = 1;

        std::vector<PyObject*> items;
        std::vector<std::pair<const char*, size_t>> views;
        items.reserve(n);
        views.reserve(n);
        for (auto rec : records) {
            PyObject* item = rec.ptr();
            PyObject* s = PySequence_GetItem(item, seq_index);   // new reference
            if (s == nullptr) throw py::error_already_set();
            Py_ssize_t len = 0;
            const char* p = PyUnicode_AsUTF8AndSize(s, &len);
            Py_DECREF(s);   // `records` still holds the tuple, which holds the str
            if (p == nullptr) throw py::error_already_set();
            items.push_back(item);
            views.emplace_back(p, size_t(len));
        }

        std::vector<uint8_t> keep(n, 0);
        size_t nthread = threads > 0 ? size_t(threads) : 1;
        nthread = std::min(nthread, n);
        auto worker = [&](size_t lo, size_t hi) {
            for (size_t i = lo; i < hi; ++i)
                keep[i] = scan(views[i].first, views[i].second, min_hits) >= min_hits ? 1 : 0;
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
            if (keep[i]) out.append(py::reinterpret_borrow<py::object>(items[i]));
        return out;
    }

private:
    void add_all(const std::string& s) {
        int run = 0;
        uint64_t code = 0;
        for (size_t i = 0; i < s.size(); ++i) {
            const int8_t b = BASE_IDX[static_cast<unsigned char>(s[i])];
            if (b < 0) { run = 0; code = 0; continue; }
            code = ((code << 2) | uint64_t(b)) & mask_;
            if (++run >= k_) seeds_.push_back(code);
        }
    }

    int k_;
    int shift_;
    uint64_t mask_;
    std::vector<uint64_t> seeds_;
    std::vector<uint64_t> pre_;
};

}  // namespace

PYBIND11_MODULE(_prefilter, m) {
    m.doc() = "Exact k-mer prefilter: reject reads that cannot align, before MMseqs2 sees them.";
    py::class_<Prefilter>(m, "Prefilter")
        .def(py::init<const std::vector<std::string>&, int>(),
             py::arg("targets"), py::arg("k") = 16,
             "Index every k-mer of `targets` and of their reverse complements.")
        .def("hits", &Prefilter::hits, py::arg("sequence"), py::arg("min_hits") = 1,
             "Indexed k-mers found in `sequence`, counted no further than `min_hits`.")
        .def("mask", &Prefilter::mask,
             py::arg("sequences"), py::arg("min_hits") = 1, py::arg("threads") = 1,
             "1 for each sequence with >= min_hits indexed k-mers, 0 otherwise.")
        .def("filter", &Prefilter::filter,
             py::arg("records"), py::arg("seq_index") = 1, py::arg("min_hits") = 1,
             py::arg("threads") = 1,
             "The subset of `records` whose element `seq_index` has >= min_hits indexed k-mers.")
        .def_property_readonly("size", &Prefilter::size, "Distinct indexed k-mers.")
        .def_property_readonly("k", &Prefilter::k, "k-mer length.");
}
