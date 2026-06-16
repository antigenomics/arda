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

import static com.antigenomics.migmap.mutation.SubRegion.*

class MutationExtractorTest {
    @Test
    void test() {
        def alignment = new Alignment(0,
                //0000000001  11111111122222222223333
                //1234567890  12345678901234567890123
                "AGATCGATCGA--CTGCTACGACTGCATGACTCAAT", 0,
                //0000000001111111111222222   2222333
                //1234567890123456789012345   6789012
                "AAATCGATCGAAACTGCTACGACTGC---ACTCAGT")

        def posInReadExp = [1, 11, 24, 32]
        def expectedCodes = ["S1:A>G", "D11:AA", "I26:ATG", "S31:G>A"]
        def mutations = MutationExtractor.extract(alignment).mutations

        println mutations

        assert expectedCodes.size() == mutations.size()

        mutations.eachWithIndex { it, i ->
            assert it.toString() == expectedCodes[i]
            assert it.posInRead == posInReadExp[i]
        }
    }

    @Test
    void markupTest() {
        def query = "CAGCTGGAGTTGGTACAGTCTGGGGCTGAGGAGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGATCCACATTCAGC" +
                "GGCCACTTTATGCACTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGGTGGATCAACTCTTACAGTGGTGCCACAAAGTAT" +
                "GCACAGAAGTTTCAGGGCAGGGTCACCATGACCAGGGACACGTCCATGACCACAATCTACATGGAGCTGAGCGGACTCACATCTGACGAC" +
                "ACGGCCGTGTATTTTTGTACCAGA"

        def segment = BlastTestUtil.annotatedSegmentDatabase.segments["IGHV1-2*02"]
        def regionMarkup = new RegionMarkup(75, 99, 150, 174, 287, 288)
        def alignment = new Alignment(0, query, 0, segment.sequence.substring(0, query.length()))
        def mutations = new MutationExtractor(segment, alignment, regionMarkup).mutations

        def expectedSubRegions = [FR1, FR1, FR1, FR1, FR1,
                                  CDR1, CDR1, CDR1, CDR1, CDR1,
                                  FR2,
                                  CDR2, CDR2, CDR2,
                                  FR3, FR3, FR3, FR3, FR3, FR3, FR3, FR3, FR3, FR3, FR3,
                                  CDR3, CDR3]

        mutations.eachWithIndex { it, i ->
            assert it.subRegion == expectedSubRegions[i]
        }
    }
}
