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
class FastaReader implements OutputPort<Read> {
    final BufferedReader reader
    String header = ""

    FastaReader(String fileName, boolean resource = false) {
        reader = Util.getReader(fileName, resource)
    }

    @Override
    synchronized Read take() {
        if (header == "") {
            // first sequence
            header = reader.readLine()
        }

        if (!header) {
            // EOF
            return null
        }

        if (!header.startsWith(">")) {
            throw new RuntimeException("Bad FASTA header")
        }

        def seq = ""
        while (true) {
            def line = reader.readLine()
            if (line == null) {
                // EOF next read, return last sequence
                def _header = header
                header = null
                return new Read(_header, seq)
            } else if (line.startsWith(">")) {
                // reset header and return current read
                def _header = header
                header = line
                return new Read(_header, seq)
            } else {
                // sequence line
                seq += line
            }
        }
    }
}
