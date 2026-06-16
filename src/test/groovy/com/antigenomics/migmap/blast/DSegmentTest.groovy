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

import com.antigenomics.migmap.genomic.Segment
import org.junit.Ignore
import org.junit.Test

import static com.antigenomics.migmap.blast.BlastTestUtil.*

class DSegmentTest {
    @Test
    void igblastOnlyTest() {
        def seq = "GACACGGCCGTATATTTCTGTGCGAGGGATCAGGGGCCACGAGACCACCCCAGTTCATCGAATTGGGGCCAGGGAACCCTGGTCACCGTC"
        def read = toRead(seq)
        def chunk = toChunk(read)
        def mapping = toMapping(chunk)
        def readMapping = toReadMapping(mapping, read)

        assert readMapping.mapping.dSegment.name == "IGHD1-14*01"
        assert readMapping.mapping.cdr3Markup.dStart == 25
        assert readMapping.mapping.cdr3Markup.dEnd == 30
    }

    @Ignore
    @Test
    void dRecoveryTest() {
        //         ................A........--------------------------------------------
        //         ---------------------------------------------------------------------
        //         ---------------------------------------------........................
        //         ---------------------------------------------------------------------
        //         ---------------------------------------------------------------------
        //         000000000011111111112222222222333333333344444444445555555555666666666
        //         012345678901234567890123456789012345678901234567890123456789012345678
        def seq = "GACACGGCCGTATATTTCTGTGCGAGGGATCAGGGGCCACGAGCCGGCCAGGGAACCCTGGTCACCGTC"
        def read = toRead(seq)
        def chunk = toChunk(read)
        def mapping = toMapping(chunk)
        def readMapping = toReadMapping(mapping, read)

        assert mapping.dSegment == Segment.DUMMY_D
        assert readMapping.mapping.dSegment.name == "IGHD2-15*01"
        assert readMapping.mapping.cdr3Markup.dStart == 8
        assert readMapping.mapping.cdr3Markup.dEnd == 12
    }
}
