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

import com.antigenomics.migmap.pipeline.Util
import groovy.transform.CompileStatic

@CompileStatic
class FastqReader implements OutputPort<Read> {
    final BufferedReader reader

    FastqReader(String fileName, boolean resource = false) {
        reader = Util.getReader(fileName, resource)
    }

    @Override
    synchronized Read take() {
        def header = reader.readLine()
        if (!header) {
            return null
        }
        if (!header.startsWith("@")) {
            throw new RuntimeException("Bad FASTQ header: $header.")
        }

        def seq = reader.readLine()
        if (seq =~ /[^ATGCN]/) {
            throw new RuntimeException("Bad sequence: $seq.")
        }
        reader.readLine()

        def qual = reader.readLine()

        new Read(header, seq, qual)
    }
}
