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

import com.antigenomics.migmap.genomic.SegmentDatabase
import com.antigenomics.migmap.io.Read
import com.antigenomics.migmap.mapping.Mapping
import com.antigenomics.migmap.mapping.ReadMapping

class BlastTestUtil {
    static final SegmentDatabase annotatedSegmentDatabase
    static final BlastInstanceFactory blastInstanceFactory
    private static final BlastInstance _instance

    static {
        blastInstanceFactory = new BlastInstanceFactory("data/", "human", ["IGH", "IGL"])
        blastInstanceFactory.annotateV()
        annotatedSegmentDatabase = blastInstanceFactory.segmentDatabase
        _instance = blastInstanceFactory.create()
        _instance.put(null)
        _instance.close()
    }


    static Mapping toMapping(String chunk) {
        _instance.parser.parse(chunk)
    }

    static ReadMapping toReadMapping(Mapping mapping, Read read) {
        _instance.createReadMapping(mapping, read)
    }

    static Read toRead(String seq) {
        new Read("@test", seq, "I" * seq.length())
    }

    static String toChunk(Read read) {
        def instance = blastInstanceFactory.create()

        instance.put(read)
        instance.put(null)

        instance.nextChunk()
    }

    static ReadMapping fullMap(String seq) {
        def read = toRead(seq)
        toReadMapping(toMapping(toChunk(read)), read)
    }
}
