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

import org.junit.Test

class FastqReaderTest {
    @Test
    void test() {
        def reader = new FastqReader("sample.fastq.gz", true)

        def nReads = 0
        def read
        while ((read = reader.take()) != null) {
            assert read.header && read.header.length() > 0
            assert read.seq && read.seq.length() > 0
            assert read.qual && read.qual.length() > 0
            nReads++
        }

        assert nReads == 1000
    }
}
