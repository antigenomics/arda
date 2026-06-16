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

package com.antigenomics.migmap.pipeline

import com.antigenomics.migmap.clonotype.Clonotype
import com.antigenomics.migmap.blast.BlastTestUtil
import com.antigenomics.migmap.genomic.SegmentDatabase
import com.antigenomics.migmap.io.ClonotypeOutput
import com.antigenomics.migmap.io.FastqReader
import com.antigenomics.migmap.io.InputPortMerge
import com.antigenomics.migmap.io.PlainTextOutput
import com.antigenomics.migmap.io.ReadMappingOutput
import com.antigenomics.migmap.mapping.ReadMappingDetailsProvider
import com.antigenomics.migmap.mapping.ReadMappingFilter

class PipelineTestCache {
    static final PipelineTestCache INSTANCE = new PipelineTestCache()

    final Pipeline pipeline
    final ClonotypeOutput clonotypeOutput
    final SegmentDatabase segmentDatabase

    final File byReadOutputFile = new File("byread.tmp.txt"),
               clonotypeOutputFile = new File("clonotypes.tmp.txt")

    private PipelineTestCache() {
        def reader = new FastqReader("sample.fastq.gz", true)
        def filter = new ReadMappingFilter()
        clonotypeOutput = new ClonotypeOutput(new PlainTextOutput(new FileOutputStream(clonotypeOutputFile)))
        pipeline = new Pipeline(reader,
                BlastTestUtil.blastInstanceFactory,
                new InputPortMerge(
                        new ReadMappingOutput(new PlainTextOutput(new FileOutputStream(byReadOutputFile)),
                                new ReadMappingDetailsProvider()),
                        clonotypeOutput),
                filter)

        pipeline.run()

        segmentDatabase = pipeline.blastInstanceFactory.segmentDatabase

        clonotypeOutputFile.deleteOnExit()
        byReadOutputFile.deleteOnExit()
    }

    List<Clonotype> getClonotypes() {
        clonotypeOutput.clonotypeAccumulator.clonotypes
    }
}
