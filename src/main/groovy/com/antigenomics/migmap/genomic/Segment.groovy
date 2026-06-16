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

package com.antigenomics.migmap.genomic

import com.antigenomics.migmap.pipeline.Util
import com.antigenomics.migmap.mapping.RegionMarkup
import groovy.transform.CompileStatic

@CompileStatic
class Segment implements Comparable<Segment>, Serializable {
    static
    final Segment DUMMY_J = new Segment(null, SegmentType.J, "DUMMY", Util.MY_NA, "AAAAAAAAAAAAAAAAAAAAAAAAA", -1),
                  DUMMY_D = new Segment(null, SegmentType.D, "DUMMY", Util.MY_NA, "AAAAAAAAAAAAAAAAAAAAAAAAA", -1)

    static boolean isDummy(Segment segment) {
        segment == null || segment.name == Util.MY_NA
    }

    final String gene, name, sequence, regexName
    final SegmentDatabase parent
    final int referencePoint
    final SegmentType type

    RegionMarkup regionMarkup = null

    Segment(SegmentDatabase parent, SegmentType type, String gene, String name, String sequence, int referencePoint) {
        this.parent = parent
        this.type = type
        this.gene = gene
        this.name = name
        this.regexName = name.replace(".", "\\.").replace("*", "\\*")
        this.sequence = sequence
        this.referencePoint = referencePoint
    }

    String toFastaString() {
        ">$name\n$sequence"
    }

    int getFrame() {
        referencePoint % 3
    }

    @Override
    boolean equals(o) {
        if (this.is(o)) return true
        if (getClass() != o.class) return false

        Segment segment = (Segment) o

        name == segment.name
    }

    @Override
    int hashCode() {
        name.hashCode()
    }

    @Override
    public String toString() {
        name
    }

    @Override
    int compareTo(Segment o) {
        this.name.compareTo(o.name)
    }
}
