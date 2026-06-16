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

import com.antigenomics.migmap.pipeline.Util
import groovy.transform.CompileStatic

@CompileStatic
class Alignment {
    final int qstart, sstart
    final int qLength, sLength
    final String sseq, qseq

    Alignment(int qstart, String qseq, int sstart, String sseq) {
        this.qstart = qstart
        this.sstart = sstart
        this.sseq = sseq
        this.qseq = qseq
        this.qLength = Util.removeGaps(qseq).length()
        this.sLength = Util.removeGaps(sseq).length()
    }

    int getSend() {
        sstart + sLength
    }

    int getQend() {
        qstart + qLength
    }
}
