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

import com.antigenomics.migmap.mapping.Cdr3Markup
import com.antigenomics.migmap.mapping.Truncations
import org.junit.Test

class PSegmentTest {
    @Test
    void noPTest() {
        //                    V1     D0   D1           J0
        //                    |      |    |            |
        //             000000000011111111112222222222333333333344444
        //             012345678901234567890123456789012345678901234
        def cdr3seq = "TGTGCGAGGTGGCTTGGGGAAGACATTCGGACCTTTGACTCCTGG"
        def cdr3Markup = new Cdr3Markup(7, 14, 19, 32)
        def truncations = new Truncations(0, 0, 0, 0)

        def pSegments = PSegmentSearcher.search(cdr3Markup, truncations, cdr3seq)

        assert pSegments.pSegmentV == -1
        assert pSegments.pSegmentJ == -1
    }

    @Test
    void okPTest() {
        //                    V1     D0   D1           J0
        //                    |    . |    |            |
        //             000000000011111111112222222222333333333344444
        //             012345678901234567890123456789012345678901234
        def cdr3seq = "TGTGCGAGCTCGCTTGGGGAAGACATTCGAAGCTTTGACTCCTGG"
        def cdr3Markup = new Cdr3Markup(7, 14, 19, 32)
        def truncations = new Truncations(0, 0, 0, 0)

        def pSegments = PSegmentSearcher.search(cdr3Markup, truncations, cdr3seq)

        assert pSegments.pSegmentV == 12
        assert pSegments.pSegmentJ == 26
    }

    @Test
    void twoConsMMPTest() {
        //                    V1     D0   D1           J0
        //                    |      |    |            |
        //             000000000011111111112222222222333333333344444
        //             012345678901234567890123456789012345678901234
        def cdr3seq = "TGTGCGAGGTGGCTTGGGGAAGGGAGAAAAAGCTTTGACTCCTGG"
        //                                   GGAGTCAAAG
        def cdr3Markup = new Cdr3Markup(7, 14, 19, 32)
        def truncations = new Truncations(0, 0, 0, 0)

        def pSegments = PSegmentSearcher.search(cdr3Markup, truncations, cdr3seq)

        assert pSegments.pSegmentJ == 28
    }

    @Test
    void singleBadSepMMPTest() {
        //                    V1     D0   D1           J0
        //                    |      |    |            |
        //             000000000011111111112222222222333333333344444
        //             012345678901234567890123456789012345678901234
        def cdr3seq = "TGTGCGAGGTGGCTTGGGGAAGGGTGCCAAAGCTTTGACTCCTGG"
        //                                   GGAGTCAAAG
        def cdr3Markup = new Cdr3Markup(7, 14, 19, 32)
        def truncations = new Truncations(0, 0, 0, 0)

        def pSegments = PSegmentSearcher.search(cdr3Markup, truncations, cdr3seq)

        assert pSegments.pSegmentJ == 27
    }
}
