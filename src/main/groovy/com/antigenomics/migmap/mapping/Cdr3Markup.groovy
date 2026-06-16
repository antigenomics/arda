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

import groovy.transform.CompileStatic

@CompileStatic
class Cdr3Markup implements Serializable {
    final int vEnd, dStart, dEnd, jStart   // within CDR3 coordinates, used for output

    Cdr3Markup(int vEnd, int dStart, int dEnd, int jStart) {
        this.vEnd = vEnd
        this.dStart = dStart
        this.dEnd = dEnd
        this.jStart = jStart
    }

    int getInsertSizeVJ() {
        dStart < 0 ? Math.max(0, jStart - vEnd - 1) : (insertSizeVD + insertSizeDJ)
    }

    int getInsertSizeVD() {
        dStart < 0 ? Math.max(0, dStart - vEnd - 1) : -1
    }

    int getInsertSizeDJ() {
        dStart < 0 ? Math.max(0, jStart - dEnd - 1) : -1
    }

    static final String OUTPUT_HEADER = "v.end.in.cdr3\td.start.in.cdr3\td.end.in.cdr3\tj.start.in.cdr3"

    @Override
    public String toString() {
        [vEnd, dStart, dEnd, jStart].join("\t")
    }
}
