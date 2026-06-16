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

package com.antigenomics.migmap.clonotype

import com.antigenomics.migmap.pipeline.Util
import com.antigenomics.migmap.blast.PSegments
import com.antigenomics.migmap.genomic.Segment
import com.antigenomics.migmap.mapping.Cdr3Markup
import com.antigenomics.migmap.mapping.ReadMapping
import com.antigenomics.migmap.mapping.Truncations
import com.antigenomics.migmap.mutation.Mutation
import com.antigenomics.migmap.mutation.MutationConverter
import com.antigenomics.migmap.mutation.MutationFormatter
import groovy.transform.CompileStatic

@CompileStatic
class Clonotype implements Comparable<Clonotype>, Serializable {
    final String cdr3nt, cdr3aa
    final Segment vSegment, dSegment, jSegment
    final List<Mutation> mutations
    final long count
    final double freq
    final byte[] cdrInsertQual, mutationQual
    final boolean hasCdr3, complete, inFrame, noStop, canonical
    final Cdr3Markup cdr3Markup
    final Truncations truncations
    final PSegments pSegments
    final ReadMapping representativeMapping

    Clonotype(String cdr3nt, String cdr3aa,
              Segment vSegment, Segment dSegment, Segment jSegment,
              List<Mutation> mutations,
              long count, double freq,
              byte[] cdrInsertQual, byte[] mutationQual,
              boolean hasCdr3, boolean complete, boolean inFrame, boolean noStop, boolean canonical,
              Cdr3Markup cdr3Markup, Truncations truncations, PSegments pSegments) {
        this.cdr3nt = cdr3nt
        this.cdr3aa = cdr3aa
        this.vSegment = vSegment
        this.dSegment = dSegment
        this.jSegment = jSegment
        this.mutations = mutations
        this.count = count
        this.freq = freq
        this.cdrInsertQual = cdrInsertQual
        this.mutationQual = mutationQual
        this.hasCdr3 = hasCdr3
        this.complete = complete
        this.inFrame = inFrame
        this.noStop = noStop
        this.canonical = canonical
        this.cdr3Markup = cdr3Markup
        this.truncations = truncations
        this.pSegments = pSegments
        this.representativeMapping = null
    }

    Clonotype(ClonotypeKey key, ClonotypeData data, long total) {
        this.representativeMapping = key.representativeMapping
        this.cdr3nt = representativeMapping.cdr3nt
        this.cdr3aa = representativeMapping.cdr3aa
        this.vSegment = key.vSegment
        this.dSegment = key.dSegment
        this.jSegment = key.jSegment
        this.mutations = key.mutations
        this.count = data.count.longValue()
        this.freq = (double) count / total

        this.cdrInsertQual = new byte[data.cdrInsertQual.length()]
        def cdrQualCount = data.cdrQualCount.longValue()
        (0..<cdrInsertQual.length).each {
            def qual = (byte) ((double) data.cdrInsertQual.get(it) / cdrQualCount)
            cdrInsertQual[it] = qual
        }
        this.mutationQual = new byte[data.mutationQual.length()]
        (0..<mutationQual.length).each {
            def qual = (byte) ((double) data.mutationQual.get(it) / count)
            mutationQual[it] = qual
        }

        this.cdr3Markup = representativeMapping.mapping.cdr3Markup
        this.truncations = representativeMapping.mapping.truncations
        this.pSegments = representativeMapping.pSegments
        this.hasCdr3 = representativeMapping.mapping.hasCdr3
        this.complete = representativeMapping.mapping.complete

        this.inFrame = representativeMapping.inFrame
        this.noStop = representativeMapping.noStop
        this.canonical = representativeMapping.canonical

        MutationConverter.annotateMutationAa(this)
    }

    @Override
    int compareTo(Clonotype o) {
        -Long.compare(this.count, o.count)
    }

    static
    final String OUTPUT_HEADER = "freq\tcount\tv\td\tj\tcdr3nt\tcdr3aa\t" + MutationFormatter.OUTPUT_HEADER_NT +
            "\t" + MutationFormatter.OUTPUT_HEADER_AA +
            "\tcdr.insert.qual\tmutations.qual\t" + Cdr3Markup.OUTPUT_HEADER + "\t" + Truncations.OUTPUT_HEADER +
            "\t" + PSegments.OUTPUT_HEADER +
            "\thas.cdr3\tin.frame\tno.stop\tcomplete\tcanonical"

    @Override
    String toString() {
        [freq, count,
         vSegment.toString(), dSegment.toString(), jSegment.toString(),
         cdr3nt, cdr3aa,
         MutationFormatter.toStringNT(mutations),
         MutationFormatter.toStringAA(mutations),
         Util.qualToString(cdrInsertQual), Util.qualToString(mutationQual),
         cdr3Markup, truncations, pSegments,
         hasCdr3, inFrame, noStop, complete, canonical].join("\t")
    }
}
