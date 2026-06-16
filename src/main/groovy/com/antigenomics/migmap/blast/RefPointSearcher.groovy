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
import groovy.transform.CompileStatic

import static com.antigenomics.migmap.pipeline.Util.GAP

@CompileStatic
class RefPointSearcher {
    final static int MAX_J_FR4_TRUNCATIONS = 5, MAX_V_FR4_TRUNCATIONS = 5

    RefPointSearcher() {
    }

    static int getCdr3End(Segment jSegment, Alignment alignment) {
        convertPositionJ(jSegment.referencePoint + 4,
                alignment.qstart, alignment.qseq, alignment.sstart, alignment.sseq,
                alignment.send)
    }

    static int getCdr3Start(Segment vSegment, Alignment alignment) {
        convertPositionV(vSegment.referencePoint - 3,
                alignment.qstart, alignment.qseq, alignment.sstart, alignment.sseq,
                alignment.send, alignment.qend)
    }

    static int convertPositionJ(int pos, int qstart, String qseq, int sstart, String sseq, int send) {
        if (pos > send) // incomplete, F/W out of scope
            return -1

        def truncations = sstart - pos
        if (truncations >= 0) // conserved residue truncated, alignment starts after F/W
            return truncations >= MAX_J_FR4_TRUNCATIONS ? -1 : qstart

        convertPosition(pos, qstart, qseq, sstart, sseq)
    }

    static int convertPositionV(int pos, int qstart, String qseq, int sstart, String sseq, int send, int qend) {
        if (sstart > pos)  // incomplete, C out of scope
            return -1

        def truncations = pos - send
        if (truncations >= 0) // conserved residue truncated, alignment ends before C
            return truncations >= MAX_V_FR4_TRUNCATIONS ? -1 : qend

        convertPosition(pos, qstart, qseq, sstart, sseq)
    }

    static int convertPosition(int pos, int qstart, String qseq, int sstart, String sseq) {
        int posRel = pos - sstart, posDelta = 0 // normal, count deletions
        for (int i = 0; i < posRel; i++) {
            if (qseq[i] == GAP)
                posDelta--
            else if (sseq[i] == GAP) {
                posDelta--
                posRel++
            }
        }

        posRel + posDelta + qstart
    }
}
