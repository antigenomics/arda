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

import com.antigenomics.migmap.genomic.SegmentDatabase
import com.antigenomics.migmap.mutation.MutationType
import com.antigenomics.migmap.mutation.SubRegion
import org.junit.Test

import static com.antigenomics.migmap.blast.BlastTestUtil.*

class BlastParserTest {
    @Test
    void exactTest() {
        def chunk = "# IGBLASTN 2.2.29+\n" +
                "# Query: @MIG UMI:GGATATGCCGCTC:8|TAAGAGGGCAGTGGTATCAACGCAGAGTACGGATATTCTGAGGTCCGCTCTCTTGGGGGGCTTTCTGAGAGTCGTGGATCTCATGTGCAAGAAAATGAAGCACCTGTGGTTCTTCCTCCTGCTGGTGGCGGCTCCCAGATGGGTCCTGTCCCAGCTGCAGCTGCAGGAGTCGGGCCCAGGACTGGTGAAGCCCTCGGAGACCCTGTCCCTCAGGTGCACTGTCTCTGGGGGCTCCATGACAAGAACTACTGACTACTGGGGCTGGGTCCGCCAGCCCCCAGGGAAGGGACTGGAGTGGATTGCAAGTGTCTCTTATAGTGGGAGCACCACCTACAACCCGTCCCGGAAGAGTCGAGTCACAATCTCCCTAGACCCGTCCAGGAACGAACTCTCCCTGGAACTGAGGTCCATGACCGCCGCAGACACGGCTGTGTATTTCTGTGCGAGGTGGCTTGGGGAAGACATTCGGACCTTTGACTCCTGGGGCCAGGGAACCCTGGTCACCGTCTCTAA|#'.IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIBIBIIIIIIIIIIIIIIIIIIIIIBIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIBIIBIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIBIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIBIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIBIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIB;IIIIIIIIIIII.''\n" +
                "# Database: data//database-66ede6f0-a68a-4931-bd32-5cabd404d955/v data//database-66ede6f0-a68a-4931-bd32-5cabd404d955/d data//database-66ede6f0-a68a-4931-bd32-5cabd404d955/j\n" +
                "# Domain classification requested: imgt\n" +
                "\n" +
                "# V-(D)-J rearrangement summary for query sequence (Top V gene match, Top D gene match, Top J gene match, Chain type, stop codon, V-J frame, Productive, Strand).  Multiple equivalent top matches having the same score and percent identity, if present, are separated by a comma.\n" +
                "IGHV4-39*01\tIGHD6-19*01\tIGHJ4*02\tVH\tNo\tN/A\tN/A\t+\n" +
                "\n" +
                "# V-(D)-J junction details based on top germline gene matches (V end, V-D junction, D region, D-J junction, J start).  Note that possible overlapping nucleotides at VDJ junction (i.e, nucleotides that could be assigned to either rearranging gene) are indicated in parentheses (i.e., (TACT)) but are not included under the V, D, or J gene itself\n" +
                "GCGAG\tN/A\tGTGGCT\tTGGGGAAGACATTCGGAC\tCTTTG\t\n" +
                "\n" +
                "# Alignment summary between query and top germline V gene hit (from, to, length, matches, mismatches, gaps, percent identity)\n" +
                "FR1-IMGT\t152\t226\t75\t72\t3\t0\t96\n" +
                "CDR1-IMGT\t227\t256\t30\t22\t8\t0\t73.3\n" +
                "FR2-IMGT\t257\t307\t51\t47\t4\t0\t92.2\n" +
                "CDR2-IMGT\t308\t328\t21\t19\t2\t0\t90.5\n" +
                "FR3-IMGT\t329\t442\t114\t96\t18\t0\t84.2\n" +
                "CDR3-IMGT (germline)\t443\t447\t5\t5\t0\t0\t100\n" +
                "Total\tN/A\tN/A\t296\t261\t35\t0\t88.2\n" +
                "\n" +
                "# Hit table (the first field indicates the chain type of the hit)\n" +
                "# Fields: subject id, q. start, query seq, s. start, subject seq\n" +
                "# 3 hits found\n" +
                "V\tIGHV4-39*01\t152\tCAGCTGCAGCTGCAGGAGTCGGGCCCAGGACTGGTGAAGCCCTCGGAGACCCTGTCCCTCAGGTGCACTGTCTCTGGGGGCTCCATGACAAGAACTACTGACTACTGGGGCTGGGTCCGCCAGCCCCCAGGGAAGGGACTGGAGTGGATTGCAAGTGTCTCTTATAGTGGGAGCACCACCTACAACCCGTCCCGGAAGAGTCGAGTCACAATCTCCCTAGACCCGTCCAGGAACGAACTCTCCCTGGAACTGAGGTCCATGACCGCCGCAGACACGGCTGTGTATTTCTGTGCGAG\t1\tCAGCTGCAGCTGCAGGAGTCGGGCCCAGGACTGGTGAAGCCTTCGGAGACCCTGTCCCTCACCTGCACTGTCTCTGGTGGCTCCATCAGCAGTAGTAGTTACTACTGGGGCTGGATCCGCCAGCCCCCAGGGAAGGGGCTGGAGTGGATTGGGAGTATCTATTATAGTGGGAGCACCTACTACAACCCGTCCCTCAAGAGTCGAGTCACCATATCCGTAGACACGTCCAAGAACCAGTTCTCCCTGAAGCTGAGCTCTGTGACCGCCGCAGACACGGCTGTGTATTACTGTGCGAG\n" +
                "D\tIGHD6-19*01\t448\tGTGGCT\t11\tGTGGCT\n" +
                "J\tIGHJ4*02\t472\tCTTTGACTCCTGGGGCCAGGGAACCCTGGTCACCGTCTC\t5\tCTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTC"

        def mapping = toMapping(chunk)

        assert mapping.vSegment.name == "IGHV4-39*01"
        assert mapping.dSegment.name == "IGHD6-19*01"
        assert mapping.jSegment.name == "IGHJ4*02"

        assert mapping.cdr3Markup.vEnd == 8
        assert mapping.cdr3Markup.dStart == 8
        assert mapping.cdr3Markup.dEnd == 14
        assert mapping.cdr3Markup.jStart == 32

        assert mapping.regionMarkup.cdr1Start == 226
        assert mapping.regionMarkup.cdr1End == 256
        assert mapping.regionMarkup.cdr2Start == 307
        assert mapping.regionMarkup.cdr2End == 328
    }

    @Test
    void mutationCoordinatesTest() {
        // def seq = "GACACGGCTGTGTATTTCTGTGCGAGATATAGTGACTACGATTACCGGACCTTTGACTCCTGGGGCCAGGGAACCCTGGTCACCGTCTC"
        //            ................A.........--------------------------------------------------------------- IGHV3-33*01
        //            --------------------------........T..........-------------------------------------------- IGHD5/OR15-5a*01
        //            --------------------------------------------------........A.............................. IGHJ4*02
        //            78901234567890123456789012345678901234567890123456789012345678901234567890123456789012345
        //            66677777777778888888888999999999900000000001111111111222222222233333333334444444444555555
        //            22222222222222222222222222222222233333333333333333333333333333333333333333333333333333333

        // def factory = new BlastInstanceFactory("data/", "human", ["IGH"], true, false)
        // def instance = factory.create()

        // def read = new Read(
        //        "@",
        //        seq,
        //        "I" * seq.length()
        // )

        // instance.put(read)
        // instance.put(null)
        // def readMapping = instance.take()

        def chunk = "# IGBLASTN 2.2.29+\n" +
                "# Query: @|GACACGGCTGTGTATTTCTGTGCGAGATATAGTGACTACGATTACCGGACCTTTGACTCCTGGGGCCAGGGAACCCTGGTCACCGTCTC|IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n" +
                "# Database: data//database-a52bccf8-da3b-473f-bacc-f9d34e532516/v data//database-a52bccf8-da3b-473f-bacc-f9d34e532516/d data//database-a52bccf8-da3b-473f-bacc-f9d34e532516/j\n" +
                "# Domain classification requested: imgt\n" +
                "\n" +
                "# V-(D)-J rearrangement summary for query sequence (Top V gene match, Top D gene match, Top J gene match, Chain type, stop codon, V-J frame, Productive, Strand).  Multiple equivalent top matches having the same score and percent identity, if present, are separated by a comma.\n" +
                "IGHV3-33*01\tIGHD5/OR15-5a*01\tIGHJ4*02\tVH\tNo\tN/A\tN/A\t+\n" +
                "\n" +
                "# V-(D)-J junction details based on top germline gene matches (V end, V-D junction, D region, D-J junction, J start).  Note that possible overlapping nucleotides at VDJ junction (i.e, nucleotides that could be assigned to either rearranging gene) are indicated in parentheses (i.e., (TACT)) but are not included under the V, D, or J gene itself\n" +
                "TGCGA\t(GA)\tTATAGTGACTACGATTAC\tCGGAC\tCTTTG\t\n" +
                "\n" +
                "# Alignment summary between query and top germline V gene hit (from, to, length, matches, mismatches, gaps, percent identity)\n" +
                "FR3-IMGT\t1\t21\t21\t20\t1\t0\t95.2\n" +
                "CDR3-IMGT (germline)\t22\t27\t6\t6\t0\t0\t100\n" +
                "Total\tN/A\tN/A\t27\t26\t1\t0\t96.3\n" +
                "\n" +
                "# Hit table (the first field indicates the chain type of the hit)\n" +
                "# Fields: subject id, q. start, query seq, s. start, subject seq\n" +
                "# 3 hits found\n" +
                "V\tIGHV3-33*01\t1\tGACACGGCTGTGTATTTCTGTGCGAGA\t268\tGACACGGCTGTGTATTACTGTGCGAGA\n" +
                "D\tIGHD5/OR15-5a*01\t26\tGATATAGTGACTACGATTAC\t4\tGATATAGTGTCTACGATTAC\n" +
                "J\tIGHJ4*02\t51\tCTTTGACTCCTGGGGCCAGGGAACCCTGGTCACCGTCTC\t5\tCTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTC\n" +
                "# BLAST processed 1 queries"

        def mapping = toMapping(chunk)
        def mutations = mapping.mutations

        assert mutations.size() == 3
        assert mutations[0].pos == 283
        assert mutations[1].pos == 301
        assert mutations[2].pos == 325
    }

    @Test
    void indelCoordinatesTest() {
        // def seq = "ATCCACTTGGTGATCAGCACTGAGCACCGAGGATTCACCATGGAACTGGGGCTCCGCTGGGTTTTCCTTGTTGCTATTTTAGAAGGTGTCCAGTGTGAGGTGCAGCTGGTGGAGTCTGGGGGAGGCCTGGTCAAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTATAGCATGAACTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTCTCATCCATTAGTAGTACTAGTTACATATACTACGCAGACTCAGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTATTACTGTGCGAGCGATCGGAACGGTATGGACGTCTGGGGCCAAGGGACCACGGTCACCGTCTCCTCAGGGAGTGCATCCGCCCCAACCCTTTTCCCCCTCTCTGCGTTGATACCACTG"
        // def factory = new BlastInstanceFactory("data/", "human", ["IGH"], true, false)
        // def instance = factory.create()

        // def read = new Read(
        //        "@",
        //        seq,
        //        "I" * seq.length()
        // )

        // instance.put(read)
        // instance.put(null)
        // println instance.nextChunk()

        def chunk = "# IGBLASTN 2.2.29+\n" +
                "# Query: @\n" +
                "# Database: /Users/mikesh/Programming/migmap/data/database-1d3a6a74-3b2d-413f-b443-dd0f7e9e058f/v /Users/mikesh/Programming/migmap/data/database-1d3a6a74-3b2d-413f-b443-dd0f7e9e058f/d /Users/mikesh/Programming/migmap/data/database-1d3a6a74-3b2d-413f-b443-dd0f7e9e058f/j\n" +
                "# Domain classification requested: imgt\n" +
                "\n" +
                "# V-(D)-J rearrangement summary for query sequence (Top V gene match, Top D gene match, Top J gene match, Chain type, stop codon, V-J frame, Productive, Strand).  Multiple equivalent top matches having the same score and percent identity, if present, are separated by a comma.\n" +
                "IGHV3-21*01\tIGHD1-20*01\tIGHJ6*02\tVH\tNo\tN/A\tN/A\t+\n" +
                "\n" +
                "# V-(D)-J junction details based on top germline gene matches (V end, V-D junction, D region, D-J junction, J start).  Note that possible overlapping nucleotides at VDJ junction (i.e, nucleotides that could be assigned to either rearranging gene) are indicated in parentheses (i.e., (TACT)) but are not included under the V, D, or J gene itself\n" +
                "AGCGA\tTC\tGGA\t(ACG)\tGTATG\t\n" +
                "\n" +
                "# Alignment summary between query and top germline V gene hit (from, to, length, matches, mismatches, gaps, percent identity)\n" +
                "FR1-IMGT\t97\t171\t75\t75\t0\t0\t100\n" +
                "CDR1-IMGT\t172\t195\t24\t24\t0\t0\t100\n" +
                "FR2-IMGT\t196\t246\t51\t51\t0\t0\t100\n" +
                "CDR2-IMGT\t247\t267\t24\t20\t1\t3\t83.3\n" +
                "FR3-IMGT\t268\t381\t114\t114\t0\t0\t100\n" +
                "CDR3-IMGT (germline)\t382\t389\t8\t7\t1\t0\t87.5\n" +
                "Total\tN/A\tN/A\t296\t291\t2\t3\t98.3\n" +
                "\n" +
                "# Hit table (the first field indicates the chain type of the hit)\n" +
                "# Fields: subject id, q. start, query seq, s. start, subject seq\n" +
                "# 3 hits found\n" +
                "V\tIGHV3-21*01\t97\tGAGGTGCAGCTGGTGGAGTCTGGGGGAGGCCTGGTCAAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTATAGCATGAACTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTCTCATCCAT---TAGTAGTACTAGTTACATATACTACGCAGACTCAGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTATTACTGTGCGAGCGA\t1\tGAGGTGCAGCTGGTGGAGTCTGGGGGAGGCCTGGTCAAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTATAGCATGAACTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTCTCATCCATTAGTAGTAGTAGTAGTTACATATACTACGCAGACTCAGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTATTACTGTGCGAGAGA\n" +
                "D\tIGHD1-20*01\t392\tGGAACG\t10\tGGAACG\n" +
                "J\tIGHJ6*02\t395\tACGGTATGGACGTCTGGGGCCAAGGGACCACGGTCACCGTCTCCTCA\t16\tACGGTATGGACGTCTGGGGCCAAGGGACCACGGTCACCGTCTCCTCA\n" +
                "# BLAST processed 1 queries"

        def mapping = toMapping(chunk)
        def mutations = mapping.mutations

        assert mutations[0].subRegion == SubRegion.CDR2
    }

    @Test
    void subRegionTest() {
        /*def seq = "GAGGTGCAATTGGTGGAGTCTGGGGGAACCTTGGTGCAGCCGGGGGGGTCCCTGACACTCTCCTGTGTAGTCTCTGGATTCACCTTTGACACTTATGCCATAAGCTGGGTCCGCCTGGCTCCAGGGAAGGGGCTGGAATGGGTCTCAAGTACTGGTGATAGAACCTTTGCAAACTCCGTGAAGGGCCGCTTCACCATCTCCAAAGACAAGTCCAAGAACACCGTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCCATTTATTATTGTGCGAAATGTGACTTTGGAGTCAGTGGCTGGTGTAACTGGCTCGACCCCTGGGGCCAGGGAACCCTGGTCACTGTCTCCTCA"
        def factory = new BlastInstanceFactory("data/", "human", ["IGH"], true, false)
        def instance = factory.create()

        def read = new Read(
               "@",
               seq,
               "I" * seq.length()
        )

        instance.put(read)
        instance.put(null)
        println instance.nextChunk()*/

        def chunk = "# IGBLASTN 2.2.29+\n" +
                "# Query: @\n" +
                "# Database: /Users/mikesh/Programming/migmap/data/database-e94a78dd-2068-4f1a-b8d9-088c1b06dc2d/v /Users/mikesh/Programming/migmap/data/database-e94a78dd-2068-4f1a-b8d9-088c1b06dc2d/d /Users/mikesh/Programming/migmap/data/database-e94a78dd-2068-4f1a-b8d9-088c1b06dc2d/j\n" +
                "# Domain classification requested: imgt\n" +
                "\n" +
                "# V-(D)-J rearrangement summary for query sequence (Top V gene match, Top D gene match, Top J gene match, Chain type, stop codon, V-J frame, Productive, Strand).  Multiple equivalent top matches having the same score and percent identity, if present, are separated by a comma.\n" +
                "IGHV3-23*04\tIGHD6-19*01\tIGHJ5*02\tVH\tNo\tN/A\tN/A\t+\n" +
                "\n" +
                "# V-(D)-J junction details based on top germline gene matches (V end, V-D junction, D region, D-J junction, J start).  Note that possible overlapping nucleotides at VDJ junction (i.e, nucleotides that could be assigned to either rearranging gene) are indicated in parentheses (i.e., (TACT)) but are not included under the V, D, or J gene itself\n" +
                "CGAAA\tTGTGACTTTGGAGT\tCAGTGGCTGGT\tGT\tAACTG\t\n" +
                "\n" +
                "# Alignment summary between query and top germline V gene hit (from, to, length, matches, mismatches, gaps, percent identity)\n" +
                "FR1-IMGT\t1\t75\t75\t66\t9\t0\t88\n" +
                "CDR1-IMGT\t76\t99\t24\t20\t4\t0\t83.3\n" +
                "FR2-IMGT\t100\t147\t51\t45\t3\t3\t88.2\n" +
                "CDR2-IMGT\t148\t164\t24\t13\t4\t7\t54.2\n" +
                "FR3-IMGT\t165\t276\t114\t101\t11\t2\t88.6\n" +
                "CDR3-IMGT (germline)\t277\t282\t6\t6\t0\t0\t100\n" +
                "Total\tN/A\tN/A\t294\t251\t31\t12\t85.4\n" +
                "\n" +
                "# Hit table (the first field indicates the chain type of the hit)\n" +
                "# Fields: subject id, q. start, query seq, s. start, subject seq\n" +
                "# 3 hits found\n" +
                "V\tIGHV3-23*04\t1\tGAGGTGCAATTGGTGGAGTCTGGGGGAACCTTGGTGCAGCCGGGGGGGTCCCTGACACTCTCCTGTGTAGTCTCTGGATTCACCTTTGACACTTATGCCATAAGCTGGGTCCGCCTGGCTCCAGGGAAGGGGCTGGAATGGGTCTCA---------AGTACTGGTGATAGAAC---CTTTGCAAACTCCGTGAAGGGCCGCTTCACCATCTCCAAAGACAAGTCCAAGAACACCGTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCCATTTATTATTGTGCGAAA\t1\tGAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTTAGCAGCTATGCCATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTCTCAGCTATTAGTGGTAGTGGTGGTAGCACATACTACGCAGACTCCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATTCCAAGAACACGCTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCCGTATATTACTGTGCGAAA\n" +
                "D\tIGHD6-19*01\t297\tCAGTGGCTGGT\t9\tCAGTGGCTGGT\n" +
                "J\tIGHJ5*02\t310\tAACTGGCTCGACCCCTGGGGCCAGGGAACCCTGGTCACTGTCTCCTCA\t3\tAACTGGTTCGACCCCTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCA\n" +
                "# BLAST processed 1 queries"

        def mapping = toMapping(chunk)
        def mutations = mapping.mutations

        def dels = mutations.findAll { it.type == MutationType.Deletion }

        assert dels[0].subRegion == SubRegion.FR2
        assert dels[1].subRegion == SubRegion.CDR2
    }

    @Test
    void truncationTest() {
        def chunk = "# IGBLASTN 2.2.29+\n" +
                "# Query: @MIG UMI:GGATATGCCGCTC:8|TAAGAGGGCAGTGGTATCAACGCAGAGTACGGATATTCTGAGGTCCGCTCTCTTGGGGGGCTTTCTGAGAGTCGTGGATCTCATGTGCAAGAAAATGAAGCACCTGTGGTTCTTCCTCCTGCTGGTGGCGGCTCCCAGATGGGTCCTGTCCCAGCTGCAGCTGCAGGAGTCGGGCCCAGGACTGGTGAAGCCCTCGGAGACCCTGTCCCTCAGGTGCACTGTCTCTGGGGGCTCCATGACAAGAACTACTGACTACTGGGGCTGGGTCCGCCAGCCCCCAGGGAAGGGACTGGAGTGGATTGCAAGTGTCTCTTATAGTGGGAGCACCACCTACAACCCGTCCCGGAAGAGTCGAGTCACAATCTCCCTAGACCCGTCCAGGAACGAACTCTCCCTGGAACTGAGGTCCATGACCGCCGCAGACACGGCTGTGTATTTCTGTGCGAGGTGGCTTGGGGAAGACATTCGGACCTTTGACTCCTGGGGCCAGGGAACCCTGGTCACCGTCTCTAA|#'.IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIBIBIIIIIIIIIIIIIIIIIIIIIBIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIBIIBIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIBIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIBIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIBIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIB;IIIIIIIIIIII.''\n" +
                "# Database: data//database-66ede6f0-a68a-4931-bd32-5cabd404d955/v data//database-66ede6f0-a68a-4931-bd32-5cabd404d955/d data//database-66ede6f0-a68a-4931-bd32-5cabd404d955/j\n" +
                "# Domain classification requested: imgt\n" +
                "\n" +
                "# V-(D)-J rearrangement summary for query sequence (Top V gene match, Top D gene match, Top J gene match, Chain type, stop codon, V-J frame, Productive, Strand).  Multiple equivalent top matches having the same score and percent identity, if present, are separated by a comma.\n" +
                "IGHV4-39*01\tIGHD6-19*01\tIGHJ4*02\tVH\tNo\tN/A\tN/A\t+\n" +
                "\n" +
                "# V-(D)-J junction details based on top germline gene matches (V end, V-D junction, D region, D-J junction, J start).  Note that possible overlapping nucleotides at VDJ junction (i.e, nucleotides that could be assigned to either rearranging gene) are indicated in parentheses (i.e., (TACT)) but are not included under the V, D, or J gene itself\n" +
                "GCGAG\tN/A\tGTGGCT\tTGGGGAAGACATTCGGAC\tCTTTG\t\n" +
                "\n" +
                "# Alignment summary between query and top germline V gene hit (from, to, length, matches, mismatches, gaps, percent identity)\n" +
                "FR1-IMGT\t152\t226\t75\t72\t3\t0\t96\n" +
                "CDR1-IMGT\t227\t256\t30\t22\t8\t0\t73.3\n" +
                "FR2-IMGT\t257\t307\t51\t47\t4\t0\t92.2\n" +
                "CDR2-IMGT\t308\t328\t21\t19\t2\t0\t90.5\n" +
                "FR3-IMGT\t329\t442\t114\t96\t18\t0\t84.2\n" +
                "CDR3-IMGT (germline)\t443\t447\t5\t5\t0\t0\t100\n" +
                "Total\tN/A\tN/A\t296\t261\t35\t0\t88.2\n" +
                "\n" +
                "# Hit table (the first field indicates the chain type of the hit)\n" +
                "# Fields: subject id, q. start, query seq, s. start, subject seq\n" +
                "# 3 hits found\n" +
                "V\tIGHV4-39*01\t152\tCAGCTGCAGCTGCAGGAGTCGGGCCCAGGACTGGTGAAGCCCTCGGAGACCCTGTCCCTCAGGTGCACTGTCTCTGGGGGCTCCATGACAAGAACTACTGACTACTGGGGCTGGGTCCGCCAGCCCCCAGGGAAGGGACTGGAGTGGATTGCAAGTGTCTCTTATAGTGGGAGCACCACCTACAACCCGTCCCGGAAGAGTCGAGTCACAATCTCCCTAGACCCGTCCAGGAACGAACTCTCCCTGGAACTGAGGTCCATGACCGCCGCAGACACGGCTGTGTATTTCTGTGCGAG\t1\tCAGCTGCAGCTGCAGGAGTCGGGCCCAGGACTGGTGAAGCCTTCGGAGACCCTGTCCCTCACCTGCACTGTCTCTGGTGGCTCCATCAGCAGTAGTAGTTACTACTGGGGCTGGATCCGCCAGCCCCCAGGGAAGGGGCTGGAGTGGATTGGGAGTATCTATTATAGTGGGAGCACCTACTACAACCCGTCCCTCAAGAGTCGAGTCACCATATCCGTAGACACGTCCAAGAACCAGTTCTCCCTGAAGCTGAGCTCTGTGACCGCCGCAGACACGGCTGTGTATTACTGTGCGAG\n" +
                "D\tIGHD6-19*01\t448\tGTGGCT\t11\tGTGGCT\n" +
                "J\tIGHJ4*02\t472\tCTTTGACTCCTGGGGCCAGGGAACCCTGGTCACCGTCTC\t5\tCTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTC"

        def mapping = toMapping(chunk)

        assert mapping.truncations.vDel == 3
        assert mapping.truncations.dDel5 == 10
        assert mapping.truncations.dDel3 == 5
        assert mapping.truncations.jDel == 4
    }

    @Test
    void tooMuchTruncationsTest() {
        def seq = "AAGTAGTCCTTGACCAGGCACGTGATGGTGGCCGACTCCCGCAGGTTCAGCTGCTCCC" +
                "GGGCTGGTGGCAGCAAGTAGACATCGGGCCTGTGCAGGGCCACCCCCTTGGGCCGGGA" +
                "GATGGTCTGCTTCAGTGGCGAGGGCAGGTCTGTGTGGGTCACGGTGCACGTGAACCTC" +
                "TCCCCGGAATTCCAGTCATCCTCGCAGATGCTGGCCTCACCCACGGCGCTGAAAGTGG" +
                "CATTGGGGTGGCTCTCGGAGATGTTGGTGTGGGTTTTCACAGCTTCGCCATTCTTGGC" +
                "GTTGTCTCTGGAGATGGTGAATCGGCCCTTCACTGAGTCTGCGTAGTATATGTAACTA" +
                "GTACTACTAATGGATGAGACCCACTCCAGAGACAACGCCAAGAACTCACTGTATCTGC" +
                "AAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTATTACTGTGCGAGCGATCGGAA" +
                "CGGTATGGACGTCTGGGGCCAAGGGACCACGGTCACCGTCTCCTCAGGGAGTGCATCC" +
                "GCCCCAACCCTTTTCCCCCTCTCTGCGTTGATACCACTG"
        def read = toRead(seq)
        def chunk = toChunk(read)

        def mapping = toMapping(chunk)

        assert !mapping.hasCdr3
        assert !mapping.complete
    }

    @Test
    void parserCase1Test() {
        def segmentDatabase = new SegmentDatabase("data/", "human", ["IGK"])
        def parser = new BlastParser(segmentDatabase)

        def chunk = "# IGBLASTN 2.2.29+\n" +
                "# Query: @MIG UMI:TACCGCCGCTTGT:5\n" +
                "# Database: /Users/mikesh/Programming/higblast/data/database-46f68e99-c1e5-45cc-b3cd-0934b850d3d4/v /Users/mikesh/Programming/higblast/data/database-46f68e99-c1e5-45cc-b3cd-0934b850d3d4/d /Users/mikesh/Programming/higblast/data/database-46f68e99-c1e5-45cc-b3cd-0934b850d3d4/j\n" +
                "# Domain classification requested: imgt\n" +
                "\n" +
                "# V-(D)-J rearrangement summary for query sequence (Top V gene match, Top D gene match, Top J gene match, Chain type, stop codon, V-J frame, Productive, Strand).  Multiple equivalent top matches having the same score and percent identity, if present, are separated by a comma.\n" +
                "IGKV4-1*01\tN/A\tN/A\tVH\tNo\tN/A\tN/A\t+\n" +
                "\n" +
                "# V-(D)-J junction details based on top germline gene matches (V end, V-D junction, D region, D-J junction, J start).  Note that possible overlapping nucleotides at VDJ junction (i.e, nucleotides that could be assigned to either rearranging gene) are indicated in parentheses (i.e., (TACT)) but are not included under the V, D, or J gene itself\n" +
                "ACTGT\tN/A\tN/A\tN/A\tN/A\t\n" +
                "\n" +
                "# Alignment summary between query and top germline V gene hit (from, to, length, matches, mismatches, gaps, percent identity)\n" +
                "FR3-IMGT\t382\t424\t43\t33\t10\t0\t76.7\n" +
                "Total\tN/A\tN/A\t43\t33\t10\t0\t76.7\n" +
                "\n" +
                "# Hit table (the first field indicates the chain type of the hit)\n" +
                "# Fields: subject id, q. start, query seq, s. start, subject seq\n" +
                "# 3 hits found\n" +
                "V\tIGKV4-1*01\t382\tCATGAGCAGCCTGAGAGCCGAAGACACGGCCGTATATTACTGT\t240\tCATCAGCAGCCTGCAGGCTGAAGATGTGGCAGTTTATTACTGT\n" +
                "V\tIGKV6D-41*01\t382\tCATGAGCAGCCTGAGAGCCGAAGACACGGCCGTATATTACTGT\t222\tCATCAGTAGCCTGGAAGCTGAAGATGCTGCAACATATTACTGT\n" +
                "V\tIGKV3D-7*01\t382\tCATGAGCAGCCTGAGAGCCGAAGACACGGCCGTATATTACTGT\t225\tCATCAGCAGCCTGCAGCCTGAAGATTTTGCAGTTTATTACTGT"

        def mapping = parser.parse(chunk)

        assert mapping.vSegment.name == "IGKV4-1*01"
        assert mapping.dSegment.name == "."
        assert mapping.jSegment.name == "."
    }

    @Test
    void parserCase2Test() {
        def segmentDatabase = new SegmentDatabase("data/", "human", ["TRA"])
        def parser = new BlastParser(segmentDatabase)

        def chunk = "# Query: @MIG UMI:TAACAATCTGAAC:11\n" +
                "# Database: /Users/mikesh/Programming/higblast/data/database-ba5de9d1-adfb-4cd7-b514-f52f28ab614e/v /Users/mikesh/Programming/higblast/data/database-ba5de9d1-adfb-4cd7-b514-f52f28ab614e/d /Users/mikesh/Programming/higblast/data/database-ba5de9d1-adfb-4cd7-b514-f52f28ab614e/j\n" +
                "# Domain classification requested: imgt\n" +
                "\n" +
                "# V-(D)-J rearrangement summary for query sequence (Top V gene match, Top D gene match, Top J gene match, Chain type, stop codon, V-J frame, Productive, Strand).  Multiple equivalent top matches having the same score and percent identity, if present, are separated by a comma.\n" +
                "TRAV8-6*02,TRAV8-6*01\t.,.,.\tN/A\tVB\tNo\tN/A\tN/A\t+\n" +
                "\n" +
                "# V-(D)-J junction details based on top germline gene matches (V end, V-D junction, D region, D-J junction, J start).  Note that possible overlapping nucleotides at VDJ junction (i.e, nucleotides that could be assigned to either rearranging gene) are indicated in parentheses (i.e., (TACT)) but are not included under the V, D, or J gene itself\n" +
                "TGTGC\tGAGACTGATTAGGGACGA\tTTTTT\tN/A\tN/A\t\n" +
                "\n" +
                "# Alignment summary between query and top germline V gene hit (from, to, length, matches, mismatches, gaps, percent identity)\n" +
                "FR3-IMGT\t365\t385\t21\t18\t3\t0\t85.7\n" +
                "CDR3-IMGT (germline)\t386\t387\t2\t2\t0\t0\t100\n" +
                "Total\tN/A\tN/A\t23\t20\t3\t0\t87\n" +
                "\n" +
                "# Hit table (the first field indicates the chain type of the hit)\n" +
                "# Fields: subject id, q. start, query seq, s. start, subject seq\n" +
                "# 6 hits found\n" +
                "V\tTRAV8-6*02\t365\tGACACGGCTGTGTATTACTGTGC\t253\tGACACGGCTGAGTACTTCTGTGC\n" +
                "V\tTRAV8-6*01\t365\tGACACGGCTGTGTATTACTGTGC\t253\tGACACGGCTGAGTACTTCTGTGC\n" +
                "V\tTRAV16*01\t364\tAGACACGGCTGTGTATTACTGTGC\t240\tAGACTCAGCCATGTATTACTGTGC\n" +
                "D\t.\t406\tTTTTT\t25\tTTTTT\n" +
                "D\t.\t406\tTTTTT\t24\tTTTTT\n" +
                "D\t.\t406\tTTTTT\t23\tTTTTT"

        def mapping = parser.parse(chunk)

        assert mapping.vSegment.name == "TRAV8-6*01"
        assert mapping.dSegment.name == "."
        assert mapping.jSegment.name == "."
    }
}
