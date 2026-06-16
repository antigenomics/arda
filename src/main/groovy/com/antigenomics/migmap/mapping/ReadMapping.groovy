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

import com.antigenomics.migmap.pipeline.Util
import com.antigenomics.migmap.blast.PSegments
import com.antigenomics.migmap.io.Read
import com.antigenomics.migmap.mutation.MutationFormatter
import groovy.transform.CompileStatic

@CompileStatic
class ReadMapping implements Serializable{
    final Read read
    final Mapping mapping
    final String cdr3nt, cdr3aa
    final byte[] mutationQual, cdrInsertQual
    final boolean canonical, inFrame, noStop
    final PSegments pSegments

    ReadMapping(Read read, Mapping mapping, String cdr3nt, String cdr3aa,
                byte[] mutationQual, byte[] cdrInsertQual,
                boolean canonical, boolean inFrame, boolean noStop,
                PSegments pSegments) {
        this.read = read
        this.mapping = mapping
        this.cdr3nt = cdr3nt
        this.cdr3aa = cdr3aa
        this.mutationQual = mutationQual
        this.cdrInsertQual = cdrInsertQual
        this.canonical = canonical
        this.inFrame = inFrame
        this.noStop = noStop
        this.pSegments = pSegments
    }

    byte getMinCdrInsertQual() {
        Util.minQual(cdrInsertQual)
    }

    byte getMinMutationQual() {
        Util.minQual(mutationQual)
    }

    boolean isMapped() {
        mapping != null
    }

    static
    final String OUTPUT_HEADER = "read.header\tcdr3nt\tcdr3aa\tcdr.insert.qual\tmutations.qual\t" +
            "$Mapping.OUTPUT_HEADER\t$MutationFormatter.OUTPUT_HEADER_AA\t$PSegments.OUTPUT_HEADER\tcanonical"

    @Override
    String toString() {
        [read.header, cdr3nt, cdr3aa, Util.qualToString(cdrInsertQual), Util.qualToString(mutationQual),
         mapping.toString(),
         MutationFormatter.toStringAA(mapping.mutations),
         pSegments.toString(), canonical].join("\t")
    }
}
