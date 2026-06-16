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

package com.antigenomics.migmap.mapping

import com.antigenomics.migmap.genomic.Segment
import com.antigenomics.migmap.mutation.Mutation
import com.antigenomics.migmap.mutation.MutationFormatter
import groovy.transform.CompileStatic

@CompileStatic
class Mapping implements Serializable {
    final Segment vSegment, jSegment, dSegment
    final RegionMarkup regionMarkup
    final Cdr3Markup cdr3Markup
    final Truncations truncations
    final List<Mutation> mutations
    final int vStartInRef, vStartInQuery
    final boolean rc, complete, hasCdr3, inFrame, noStop, hasD, dFound

    Mapping(Segment vSegment, Segment dSegment, Segment jSegment, int vStartInRef, int vStartInQuery,
            RegionMarkup regionMarkup, Cdr3Markup cdr3Markup, Truncations truncations,
            boolean rc, boolean complete, boolean hasCdr3, boolean inFrame, boolean noStop,
            boolean hasD, boolean dFound,
            List<Mutation> mutations) {
        // needed to extract CDR3 sequence
        this.vStartInRef = vStartInRef
        this.vStartInQuery = vStartInQuery
        this.regionMarkup = regionMarkup
        this.rc = rc

        // main data
        this.vSegment = vSegment
        this.dSegment = dSegment
        this.jSegment = jSegment

        this.cdr3Markup = cdr3Markup
        this.truncations = truncations

        this.mutations = mutations

        // misc
        this.complete = complete
        this.hasCdr3 = hasCdr3
        this.inFrame = inFrame
        this.noStop = noStop
        this.hasD = hasD
        this.dFound = dFound
    }

    static final String OUTPUT_HEADER = "v.segment\td.segment\tj.segment\t" +
            "$RegionMarkup.OUTPUT_HEADER\t" +
            "$Cdr3Markup.OUTPUT_HEADER\t" +
            "$Truncations.OUTPUT_HEADER\t" +
            "$MutationFormatter.OUTPUT_HEADER_NT\t" +
            "rc\tcomplete\thas.cdr3\tin.frame\tno.stop"


    @Override
    public String toString() {
        [vSegment.toString(), dSegment.toString(), jSegment.toString(),
         regionMarkup.toString(),
         cdr3Markup.toString(),
         truncations.toString(),
         MutationFormatter.toStringNT(mutations),
         rc, complete, hasCdr3, inFrame, noStop].join("\t")
    }
}
