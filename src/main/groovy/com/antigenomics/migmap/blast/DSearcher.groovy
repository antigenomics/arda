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

@CompileStatic
class DSearcher {
    final int minHitSize, maxHitSize
    final List<DHitTracker> hitTrackers

    DSearcher(List<Segment> segments) {
        this(segments, 2, 4) // IgBlast has a minimum exact match size of 5
    }

    DSearcher(List<Segment> segments, int minHitSize, int maxHitSize) {
        this.hitTrackers = new LinkedList<>()
        for (Segment segment : segments) {
            hitTrackers.add(new DHitTracker(segment, minHitSize, maxHitSize))
        }
        this.minHitSize = minHitSize
        this.maxHitSize = maxHitSize
    }

    AuxiliaryDMapping scan(String cdr3seq, int vEndRel, int jStartRel) {
        int regionSize = jStartRel - vEndRel - 1
        if (regionSize < minHitSize)
            return null

        def subSeq = cdr3seq.substring(vEndRel + 1, jStartRel)
        def hitSegments = new ArrayList<>()

        // iterate from largest window to smallest one
        for (int i = maxHitSize; i >= minHitSize; i--) {
            // sliding window scan
            for (int j = 0; j < subSeq.length() - i; j++) {
                String kmer = subSeq.substring(j, j + i)

                boolean hasHit = false
                for (DHitTracker hitTracker : hitTrackers) {
                    if (hitTracker.hasHit(kmer)) {
                        hasHit = true
                        hitSegments.add(hitTracker.getSegment())
                    }
                }

                if (hasHit) {
                    int from = vEndRel + 1 + j;
                    return new AuxiliaryDMapping(from,
                            from + i,
                            // binomial P-value, uniform base distribution
                            1.0d - Math.pow(1.0d - Math.pow(0.25d, i), regionSize - i + 1),
                            hitSegments
                    )
                }
            }
        }

        null
    }
}
