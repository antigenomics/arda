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

import com.antigenomics.migmap.mutation.Mutation
import com.antigenomics.migmap.mutation.MutationType
import com.antigenomics.migmap.mutation.SubRegion
import groovy.transform.CompileStatic

@CompileStatic
class RegionMarkup implements Serializable {
    final static RegionMarkup DUMMY = new RegionMarkup(-1, -1, -1, -1, -1, -1)

    final int cdr1Start, cdr1End, cdr2Start, cdr2End, cdr3Start, cdr3End // in read, used for extraction

    RegionMarkup(int cdr1Start, int cdr1End,
                 int cdr2Start, int cdr2End,
                 int cdr3Start, int cdr3End) {
        this.cdr1Start = cdr1Start
        this.cdr1End = cdr1End
        this.cdr2Start = cdr2Start
        this.cdr2End = cdr2End
        this.cdr3Start = cdr3Start
        this.cdr3End = cdr3End
    }

    SubRegion getSubRegion(Mutation mutation) {
        getSubRegion(mutation.type == MutationType.Deletion ? mutation.posInRead - 1 : mutation.posInRead)
    }

    SubRegion getSubRegion(int pos) {
        if (pos < cdr1Start) {
            return SubRegion.FR1
        } else if (pos < cdr1End) {
            return SubRegion.CDR1
        } else if (pos < cdr2Start) {
            return SubRegion.FR2
        } else if (pos < cdr2End) {
            return SubRegion.CDR2
        } else if (pos < cdr3Start) {
            return SubRegion.FR3
        } else {
            return SubRegion.CDR3
        }
    }

    boolean isComplete() {
        cdr1Start >= 0 && cdr1End >= 0 &&
                cdr2Start >= 0 && cdr2End >= 0 &&
                cdr3Start >= 0
    }

    static final String OUTPUT_HEADER = "cdr1.start.in.read\tcdr1.end.in.read\t" +
            "cdr2.start.in.read\tcdr2.end.in.read\t" +
            "cdr3.start.in.read\tcdr3.end.in.read"

    @Override
    public String toString() {
        [cdr1Start, cdr1End, cdr2Start, cdr2End, cdr3Start, cdr3End].join("\t")
    }
}
