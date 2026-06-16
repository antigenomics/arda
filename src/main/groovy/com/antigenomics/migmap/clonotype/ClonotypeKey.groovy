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

import com.antigenomics.migmap.genomic.Segment
import com.antigenomics.migmap.mapping.ReadMapping
import com.antigenomics.migmap.mutation.Mutation
import groovy.transform.CompileStatic

@CompileStatic
class ClonotypeKey {
    final String cdr3nt
    protected final ReadMapping representativeMapping

    ClonotypeKey(ReadMapping readMapping) {
        this.cdr3nt = readMapping.cdr3nt
        this.representativeMapping = readMapping
    }

    Segment getvSegment() {
        representativeMapping.mapping.vSegment
    }

    Segment getdSegment() {
        representativeMapping.mapping.dSegment
    }

    Segment getjSegment() {
        representativeMapping.mapping.jSegment
    }

    List<Mutation> getMutations() {
        representativeMapping.mapping.mutations
    }

    @Override
    boolean equals(o) {
        if (this.is(o)) return true
        if (getClass() != o.class) return false

        ClonotypeKey that = (ClonotypeKey) o

        cdr3nt == that.cdr3nt &&
                vSegment == that.vSegment &&
                jSegment == that.jSegment &&
                mutations == that.mutations
    }

    @Override
    int hashCode() {
        int result
        result = cdr3nt.hashCode()
        result = 31 * result + vSegment.hashCode()
        result = 31 * result + jSegment.hashCode()
        31 * result + mutations.hashCode()
    }
}
