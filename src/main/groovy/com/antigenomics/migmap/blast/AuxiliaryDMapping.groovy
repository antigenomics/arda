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
import com.antigenomics.migmap.mapping.Cdr3Markup
import com.antigenomics.migmap.mapping.Mapping
import groovy.transform.CompileStatic

@CompileStatic
class AuxiliaryDMapping {
    final int dStart, dEnd
    final double score
    final ArrayList<Segment> dSegments

    AuxiliaryDMapping(int dStart, int dEnd, double score, ArrayList<Segment> dSegments) {
        this.dStart = dStart
        this.dEnd = dEnd
        this.score = score
        this.dSegments = dSegments
    }

    Mapping updateMapping(Mapping old) {
        Collections.sort(this.dSegments)
        new Mapping(old.vSegment, dSegments[0], old.jSegment,
                old.vStartInRef, old.vStartInQuery, old.regionMarkup,
                new Cdr3Markup(old.cdr3Markup.vEnd, dStart, dEnd, old.cdr3Markup.jStart),
                old.truncations, old.rc, old.complete, old.hasCdr3, old.inFrame, old.noStop, old.hasD, true, old.mutations)
    }
}
