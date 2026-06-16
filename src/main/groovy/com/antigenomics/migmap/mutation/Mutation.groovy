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

package com.antigenomics.migmap.mutation

import com.antigenomics.migmap.genomic.Segment
import groovy.transform.CompileStatic

@CompileStatic
class Mutation implements Serializable {
    final int startInRead, endInRead
    final String ntFrom, ntTo
    String aaFrom, aaTo
    int aaPos

    int start, end
    final MutationType type
    Segment parent
    SubRegion subRegion

    static Mutation fromString(String signature, String aaSignature = null) {
        MutationType type = MutationType.byShortName(signature[0])

        def tmp = signature.substring(1).split(":")

        int pos = tmp[0].toInteger()
        String ntFrom = "", ntTo = ""
        switch (type) {
            case MutationType.Substitution:
                tmp = tmp[1].split(">")
                ntFrom = tmp[0]
                ntTo = tmp[1]
                break
            case MutationType.Deletion:
                ntFrom = tmp[1]
                break
            case MutationType.Insertion:
                ntTo = tmp[1]
                break
        }

        def mutation = new Mutation(type, pos, pos, -1, -1, ntFrom, ntTo)

        if (aaSignature) {
            def fromToAA = aaSignature.split(":")[1].split(">")
            switch (type) {
                case MutationType.Substitution:
                    mutation.aaFrom = fromToAA[0]
                    mutation.aaTo = fromToAA[1]
                    break
                case MutationType.Deletion:
                    mutation.aaFrom = fromToAA[0]
                    mutation.aaTo = ""
                    break
                case MutationType.Insertion:
                    mutation.aaFrom = ""
                    mutation.aaTo = fromToAA[0]
                    break
            }
            mutation.aaPos = (int) (pos / 3i)
        }

        mutation
    }

    Mutation(MutationType type,
             int start, int end, int startInRead, int endInRead,
             String ntFrom, String ntTo) {
        this.type = type
        this.end = end
        this.startInRead = startInRead
        this.endInRead = endInRead
        this.ntFrom = ntFrom
        this.ntTo = ntTo
        this.start = start
    }

    int getPos() {
        type == MutationType.Insertion ? end : start
    }

    int getPosInRead() {
        type == MutationType.Deletion ? endInRead : startInRead // deletion in ref == insertion in query
    }

    @Override
    public String toString() {
        // insertion occur before "position", so that seq[position, ...] is shifted forward
        type.shortName + pos + ":" +
                (type == MutationType.Insertion ? "" : ntFrom) +
                (type == MutationType.Substitution ? ">" : "") +
                (type == MutationType.Deletion ? "" : ntTo)
    }

    public String toStringAa() {
        // insertion occur before "position", so that seq[position, ...] is shifted forward
        type.shortName + aaPos + ":" +
                (type == MutationType.Insertion ? "" : aaFrom) +
                (type == MutationType.Substitution ? ">" : "") +
                (type == MutationType.Deletion ? "" : aaTo)
    }

    boolean equals(o) {
        if (this.is(o)) return true
        if (getClass() != o.class) return false

        Mutation mutation = (Mutation) o

        if (start != mutation.start) return false
        if (ntFrom != mutation.ntFrom) return false
        if (ntTo != mutation.ntTo) return false
        if (parent != mutation.parent) return false

        return true
    }

    int hashCode() {
        int result
        result = ntFrom.hashCode()
        result = 31 * result + ntTo.hashCode()
        result = 31 * result + start
        result = 31 * result + (parent != null ? parent.hashCode() : 0)
        return result
    }
}
