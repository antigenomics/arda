/*
 * Copyright 2014-2015 Mikhail Shugay
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package com.antigenomics.migmap.blast

import org.junit.Test

import static com.antigenomics.migmap.blast.BlastTestUtil.*

class RefPointSearcherTest {
    @Test
    void cdr3EndTest() {
        def seq = "CAGGTGCAGCTGCAGGAGTCGGGCCCAGGACTGGTAAAGCCTTCGGAGACCCTGTCCCTCACCTGCGGTGTCTCTGGTTACTCCATAAGCAGTGGTTACTACTGGGCCTGGATCCGGCAGCCCCCAGGGAAGGGGCTGGAGTGGGTTGCGACTATCTATCATGATGGAAGATCCTACTACAACCCGTCCCTGGAAAGTCGAGTCACCATATCAGTAGACACGTCCAAGAACCAGTTCTCCCTGAGGCTGACTTCTGTGACCGCCGCAGACACGGCCGTATATTTCTGTGCGAGGGATCAGGGGCCACGAGACCACCCCAGTTCATCGTTTTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCA"
        def read = toRead(seq)
        def chunk = toChunk(read)

        def mapping = toMapping(chunk)
        assert mapping.hasCdr3
        assert mapping.complete
        assert mapping.regionMarkup.cdr3End == 333

        def readMapping = toReadMapping(mapping, read)
        assert readMapping.canonical
        assert readMapping.cdr3nt == "TGTGCGAGGGATCAGGGGCCACGAGACCACCCCAGTTCATCGTTTTGG"
    }

    @Test
    void cdr3StartTest() {
        def seq = "CAGCTGCGGCTGCAGGAGTCGGGCCCAGGACAGGTGAAGCCTTCGGAGACCCTGTCCCTCACTTGCACGGTCTCTGCTGGCGCCATCACCAGTGATCATTACTACTGGGGCTGGGTCCGCCAGCGCCCAGGGAAGGGACTAGAGTGGATTGGGAGTGTCCATAATAGTGGGAGCACCTCCTACAACCCGTCCCTTCAGAGTCGAATCACCATGTCCATAGACATGTCGAAGAACCACTTCTCCCTGCGGCTGACTTCCGTGACCGCCGCAGACACGGGTGTTTATTACTGCACCACCTATGGTTATTATTACTATTACGGCATGGACGTCTGGGGCCAAGGGACCACGGTCACCGTCTCCTCA"
        def read = toRead(seq)
        def chunk = toChunk(read)

        def mapping = toMapping(chunk)

        def readMapping = toReadMapping(mapping, read)
        assert readMapping.cdr3nt == "TGCACCACCTATGGTTATTATTACTATTACGGCATGGACGTCTGG"
    }

    @Test
    void cdr3StartTest2() {
        def seq = "CAGGTGCAGCTGATGGAGTCTGGGGGAGGCGTTGTCCAGCCTGGGAGGTCCCTGAGACTCTCCTGTGAAGCCTCTGGATTCACCTTCAGTCGTTCTTCTATGCACTGGGTCCGCCAGGCTCCAGGCAAGGGGCTGGAGTGGTTGGCACTTATTTCACGTGATGGGAACTATGAGAGTTACGCAGACTCCGTGAAGGGCCGATTCTCGATCTCCAGAGACAACTCCAAGAACACTCTCTTTCTGCAAGTGGACGGCCTGAGAGCTGAGGACACGGCTGTGTATTATTGTCTTGGTGACAACGGCCATTGGGGCCAGGGAACCCTGCTCAGCGTCTCCTCA"
        def read = toRead(seq)
        def chunk = toChunk(read)

        def mapping = toMapping(chunk)

        def readMapping = toReadMapping(mapping, read)
        assert readMapping.cdr3nt == "TGTCTTGGTGACAACGGCCATTGG"
    }

    @Test
    void borderlineCase() {
        def chunk = "# IGBLASTN 2.2.29+\n" +
                "# Query: @M01691:112:000000000-AFF3F:1:1106:14463:20458 2:N:0:1 UMI:ATGTGTACCGGG:GGGGGGGGGGGG\n" +
                "# Database: /Users/mikesh/Programming/migmap/data/database-f1894158-8d79-41c0-8288-0bd8217f7c4f/v /Users/mikesh/Programming/migmap/data/database-f1894158-8d79-41c0-8288-0bd8217f7c4f/d /Users/mikesh/Programming/migmap/data/database-f1894158-8d79-41c0-8288-0bd8217f7c4f/j\n" +
                "# Domain classification requested: imgt\n" +
                "\n" +
                "# V-(D)-J rearrangement summary for query sequence (Top V gene match, Top D gene match, Top J gene match, Chain type, stop codon, V-J frame, Productive, Strand).  Multiple equivalent top matches having the same score and percent identity, if present, are separated by a comma.\n" +
                "IGHV3-21*01\tIGHD5-5*01\tIGHJ6*04\tVH\tYes\tN/A\tNo\t+\n" +
                "\n" +
                "# V-(D)-J junction details based on top germline gene matches (V end, V-D junction, D region, D-J junction, J start).  Note that possible overlapping nucleotides at VDJ junction (i.e, nucleotides that could be assigned to either rearranging gene) are indicated in parentheses (i.e., (TACT)) but are not included under the V, D, or J gene itself\n" +
                "GAGAG\tTCCCAGAGGC\tGGATACAGCTATGGTT\tTTGAG\tTTACT\t\n" +
                "\n" +
                "# Alignment summary between query and top germline V gene hit (from, to, length, matches, mismatches, gaps, percent identity)\n" +
                "FR1-IMGT\t93\t167\t75\t75\t0\t0\t100\n" +
                "CDR1-IMGT\t168\t191\t24\t24\t0\t0\t100\n" +
                "FR2-IMGT\t192\t242\t51\t51\t0\t0\t100\n" +
                "CDR2-IMGT\t243\t266\t24\t24\t0\t0\t100\n" +
                "FR3-IMGT\t267\t376\t114\t110\t0\t4\t96.5\n" +
                "CDR3-IMGT (germline)\t377\t383\t7\t7\t0\t0\t100\n" +
                "Total\tN/A\tN/A\t295\t291\t0\t4\t98.6\n" +
                "\n" +
                "# Hit table (the first field indicates the chain type of the hit)\n" +
                "# Fields: subject id, q. start, query seq, s. start, subject seq\n" +
                "# 3 hits found\n" +
                "V\tIGHV3-21*01\t93\tGAGGTGCAGCTGGTGGAGTCTGGGGGAGGCCTGGTCAAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTATAGCATGAACTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTCTCATCCATTAGTAGTAGTAGTAGTTACATATACTAC----ACTCAGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTATTACTGTGCGAGAG\t1\tGAGGTGCAGCTGGTGGAGTCTGGGGGAGGCCTGGTCAAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTATAGCATGAACTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTCTCATCCATTAGTAGTAGTAGTAGTTACATATACTACGCAGACTCAGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTATTACTGTGCGAGAG\n" +
                "D\tIGHD5-5*01\t394\tGGATACAGCTATGGTT\t3\tGGATACAGCTATGGTT\n" +
                "J\tIGHJ6*04\t415\tTTACTACTACTACTACGGTATGGACGTCTGG\t2\tTTACTACTACTACTACGGTATGGACGTCTGG\n" +
                "# BLAST processed 1 queries"

        def mapping = toMapping(chunk)

        assert mapping.complete
    }
}
