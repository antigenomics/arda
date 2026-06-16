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

package com.antigenomics.migmap.mutation

import com.antigenomics.migmap.blast.Alignment
import com.antigenomics.migmap.blast.BlastTestUtil
import com.antigenomics.migmap.mapping.RegionMarkup
import org.junit.Test

class MutationFormatterTest {
    @Test
    void ntMutationsTest() {
        def query = "CAGCTGGAGTTGGTACAGTCTGGGGCTGAGGAGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGATCCACATTCAGC" +
                "GGCCACTTTATGCACTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGGTGGATCAACTCTTACAGTGGTGCCACAAAGTAT" +
                "GCACAGAAGTTTCAGGGCAGGGTCACCATGACCAGGGACACGTCCATGACCACAATCTACATGGAGCTGAGCGGACTCACATCTGACGAC" +
                "ACGGCCGTGTATTTTTGTACCAGA"

        def segment = BlastTestUtil.annotatedSegmentDatabase.segments["IGHV1-2*02"]
        def regionMarkup = new RegionMarkup(75, 99, 150, 174, 287, 288)
        def alignment = new Alignment(0, query, 0, segment.sequence.substring(0, query.length()))
        def mutations = new MutationExtractor(segment, alignment, regionMarkup).mutations

        assert MutationFormatter.toStringNT(mutations).split("\t") ==
                ['S3:G>C,S6:C>G,S9:C>T,S14:G>A,S31:T>A',
                 'S79:A>C,S83:C>A,S88:C>G,S93:T>C,S97:A>T',
                 'S146:A>G',
                 'S156:C>T,S159:A>T,S169:G>C',
                 'S176:C>G,S227:C>G,S229:G>C,S234:G>A,S235:C>T,S252:A>G,S254:G>A,S257:G>C,S259:G>C,S283:A>T,S284:C>T',
                 'S288:G>A,S290:G>C']
    }

    @Test
    void aaMutationsTestShift1() {
        def query = "GCTGGAGTTGGTACAGTCTGGGGCTGAGGAGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGATCCACATTCAGC" +
                "GGCCACTTTATGCACTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGGTGGATCAACTCTTACAGTGGTGCCACAAAGTAT" +
                "GCACAGAAGTTTCAGGGCAGGGTCACCATGACCAGGGACACGTCCATGACCACAATCTACATGGAGCTGAGCGGACTCACATCTGACGAC" +
                "ACGGCCGTGTATTTTTGTACCAGA"

        def segment = BlastTestUtil.annotatedSegmentDatabase.segments["IGHV1-2*02"]
        def regionMarkup = new RegionMarkup(75 - 2, 99 - 2, 150 - 2, 174 - 2, 287 - 2, 288 - 2)
        def alignment = new Alignment(0, query, 2, segment.sequence.substring(2, 2 + query.length()))
        def mutations = new MutationExtractor(segment, alignment, regionMarkup).mutations

        assert MutationFormatter.toStringNT(mutations).split("\t") ==
                ['S3:G>C,S6:C>G,S9:C>T,S14:G>A,S31:T>A',
                 'S79:A>C,S83:C>A,S88:C>G,S93:T>C,S97:A>T',
                 'S146:A>G',
                 'S156:C>T,S159:A>T,S169:G>C',
                 'S176:C>G,S227:C>G,S229:G>C,S234:G>A,S235:C>T,S252:A>G,S254:G>A,S257:G>C,S259:G>C,S283:A>T,S284:C>T',
                 'S288:G>A,S290:G>C']
    }

    @Test
    void aaMutationsTestShift2() {
        def query = "CGCTGGAGTTGGTACAGTCTGGGGCTGAGGAGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGATCCACATTCAGC" +
                "GGCCACTTTATGCACTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGGTGGATCAACTCTTACAGTGGTGCCACAAAGTAT" +
                "GCACAGAAGTTTCAGGGCAGGGTCACCATGACCAGGGACACGTCCATGACCACAATCTACATGGAGCTGAGCGGACTCACATCTGACGAC" +
                "ACGGCCGTGTATTTTTGTACCAGA"

        def segment = BlastTestUtil.annotatedSegmentDatabase.segments["IGHV1-2*02"]
        def regionMarkup = new RegionMarkup(75 - 1, 99 - 1, 150 - 1, 174 - 1, 287 - 1, 288 - 1)
        def alignment = new Alignment(1, query.substring(1), 2, segment.sequence.substring(2, 2 + query.length() - 1))
        def mutations = new MutationExtractor(segment, alignment, regionMarkup).mutations

        assert MutationFormatter.toStringNT(mutations).split("\t") ==
                ['S3:G>C,S6:C>G,S9:C>T,S14:G>A,S31:T>A',
                 'S79:A>C,S83:C>A,S88:C>G,S93:T>C,S97:A>T',
                 'S146:A>G',
                 'S156:C>T,S159:A>T,S169:G>C',
                 'S176:C>G,S227:C>G,S229:G>C,S234:G>A,S235:C>T,S252:A>G,S254:G>A,S257:G>C,S259:G>C,S283:A>T,S284:C>T',
                 'S288:G>A,S290:G>C']
    }
}
