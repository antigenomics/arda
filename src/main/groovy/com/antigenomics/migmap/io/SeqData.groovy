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

package com.antigenomics.migmap.io

import groovy.transform.CompileStatic

@CompileStatic
class SeqData {
    int nReads = 0, nEvents = 0
    final int[] qual
    final int seqId
    String computedQual = null

    SeqData(int seqId, String qual) {
        this.seqId = seqId
        this.qual = qual ? new int[qual.length()] : null
    }

    void append(String newQual, int count) {
        nEvents++
        nReads += count
        if (this.qual) {
            def qualCharArray = newQual.toCharArray()
            for (int i = 0; i < qualCharArray.length; i++)
                qual[i] += (qualCharArray[i] - 33) * count
        }
    }

    String computeQual() {
        if (qual) {
            if (computedQual)
                return computedQual
            else {
                for (int i = 0; i < qual.length; i++)
                    qual[i] /= nReads
                def qualCharArray = new char[qual.length]
                for (int i = 0; i < qual.length; i++)
                    qualCharArray[i] = (char) (qual[i] + 33)
                return computedQual = new String(qualCharArray)
            }
        }
        return null
    }
}
