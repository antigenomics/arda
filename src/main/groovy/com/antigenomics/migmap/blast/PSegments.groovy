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

import groovy.transform.CompileStatic

@CompileStatic
class PSegments implements Serializable {
    final int pSegmentV, pSegmentJ, pSegmentD5, pSegmentD3

    PSegments(int pSegmentV, int pSegmentD5, int pSegmentD3, int pSegmentJ) {
        this.pSegmentV = pSegmentV
        this.pSegmentD5 = pSegmentD5
        this.pSegmentD3 = pSegmentD3
        this.pSegmentJ = pSegmentJ
    }

    static final String OUTPUT_HEADER = "pol.v\tpol.d.5\tpol.d.3\tpol.j"

    @Override
    String toString() {
        [pSegmentV, pSegmentD5, pSegmentD3, pSegmentJ].join("\t")
    }
}
