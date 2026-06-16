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

import com.antigenomics.migmap.mapping.ReadMapping
import groovy.transform.CompileStatic

import java.util.concurrent.atomic.AtomicLong
import java.util.concurrent.atomic.AtomicLongArray

@CompileStatic
class ClonotypeData {
    final AtomicLong count, cdrQualCount
    final AtomicLongArray mutationQual, cdrInsertQual

    ClonotypeData(ReadMapping readMapping) {
        this.count = new AtomicLong(0)
        this.cdrQualCount = new AtomicLong(0)
        this.mutationQual = new AtomicLongArray(readMapping.mutationQual.length)
        this.cdrInsertQual = new AtomicLongArray(readMapping.cdrInsertQual.length)

        update(readMapping)
    }

    void update(ReadMapping readMapping) {
        count.incrementAndGet()
        readMapping.mutationQual.eachWithIndex { byte it, int i -> mutationQual.addAndGet(i, it) }
        if (readMapping.cdrInsertQual.length == cdrInsertQual.length()) {
            // protect against very rare cases of ambiguous D alignment
            // not a big deal as this quality is used just for display
            readMapping.cdrInsertQual.eachWithIndex { byte it, int i -> cdrInsertQual.addAndGet(i, it) }
            cdrQualCount.incrementAndGet()
        }
    }
}
