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

import com.antigenomics.migmap.pipeline.Util
import com.antigenomics.migmap.genomic.Segment
import groovy.transform.CompileStatic

@CompileStatic
class DHitTracker {
    final Segment segment;
    final Set<String> kmers = new HashSet<>();

    DHitTracker(Segment segment, int minHitSize, int maxHitSize) {
        this.segment = segment;
        [segment.sequence, Util.revCompl(segment.sequence)].each {
            for (int i = minHitSize; i <= maxHitSize; i++) {
                for (int j = 0; j < it.length() - i; j++) {
                    kmers.add(it[j..<(j + i)])
                }
            }
        }
    }

    boolean hasHit(String queryKmer) {
        kmers.contains(queryKmer)
    }

    Segment getSegment() {
        segment
    }
}
