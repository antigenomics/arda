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

import com.antigenomics.migmap.genomic.Segment
import com.antigenomics.migmap.genomic.SegmentDatabase
import com.antigenomics.migmap.genomic.SegmentType
import com.antigenomics.migmap.mapping.Cdr3Markup
import groovy.transform.CompileStatic

@CompileStatic
class DSearcherBundle {
    private final Map<String, DSearcher> searchers = new HashMap<>()

    DSearcherBundle(SegmentDatabase segmentDatabase) {
        segmentDatabase.genes.each { gene ->
            searchers.put(gene,
                    new DSearcher(new ArrayList<Segment>(segmentDatabase.segments.values().findAll() { Segment it ->
                        it.type == SegmentType.D && it.gene == gene
                    })))
        }
    }

    AuxiliaryDMapping search(String gene, Cdr3Markup cdr3Markup, String cdr3nt) {
        if (cdr3Markup.jStart < 0)
            return null

        searchers[gene].scan(cdr3nt, cdr3Markup.vEnd, cdr3Markup.jStart)
    }
}
