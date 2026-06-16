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

import com.antigenomics.migmap.blast.BlastTestUtil
import com.antigenomics.migmap.io.*
import com.antigenomics.migmap.mapping.ReadMapping
import com.antigenomics.migmap.mapping.ReadMappingDetailsProvider
import com.antigenomics.migmap.mapping.ReadMappingFilter
import org.junit.Test

import java.util.concurrent.atomic.AtomicInteger

class PipelineTest {
    private static class NullOutputStream extends OutputStream {
        static final NullOutputStream INSTANCE = new NullOutputStream()

        @Override
        void write(int b) throws IOException {
        }
    }

    @Test
    void sampleTest() {
        assert PipelineTestCache.INSTANCE.pipeline.inputCount == 1000
        assert PipelineTestCache.INSTANCE.pipeline.readMappingFilter.mappedRatio >= 0.95
        assert PipelineTestCache.INSTANCE.pipeline.readMappingFilter.noCdr3Ratio <= 0.1
        assert PipelineTestCache.INSTANCE.pipeline.readMappingFilter.incompleteRatio <= 0.1
        assert PipelineTestCache.INSTANCE.pipeline.readMappingFilter.nonCanonicalRatio <= 0.1
    }

    /*
    @Test
    void shmStatTest() {
        def bos = new ByteArrayOutputStream()
        def clonotypes = clonotypeOutput1.clonotypeAccumulator.clonotypes

        ClonotypeSerializer.save(clonotypes, bos)

        def bis = new ByteArrayInputStream(bos.toByteArray())

        def loadedClonotypes = ClonotypeSerializer.load(bis)

        def postAnalysis = new Analysis(loadedClonotypes)

        postAnalysis.generateHypermutationTable("temp.shm.txt")

        assert new File("temp.shm.txt").readLines().size() == 8865

        new File("temp.shm.txt").delete()
    }

    @Test
    void serializationTest() {
        def bos = new ByteArrayOutputStream()

        def clonotypes = clonotypeOutput1.clonotypeAccumulator.clonotypes

        ClonotypeSerializer.save(clonotypes, bos)

        def bis = new ByteArrayInputStream(bos.toByteArray())

        def loadedClonotypes = ClonotypeSerializer.load(bis)

        assert !loadedClonotypes.empty
        assert clonotypes.size() == loadedClonotypes.size()
    }
    */

    @Test
    void fastaTest() {
        def reader = new FastaReader("sample.fasta.gz", true)

        def pipeline = new Pipeline(reader,
                BlastTestUtil.blastInstanceFactory,
                DummyInputPort.INSTANCE,
                ReadMappingFilter.createDummy())

        pipeline.run()

        assert pipeline.readMappingFilter.total == 100
        assert pipeline.readMappingFilter.goodRatio == 1.0
    }

    @Test
    void dAssignmentTest() {
        def reader = new FastqReader("ambiguous_d.fastq.gz", true)

        def badDCount = new AtomicInteger(),
            wrongMappingCount = new AtomicInteger()

        def pipeline = new Pipeline(reader,
                BlastTestUtil.blastInstanceFactory,
                new InputPort<ReadMapping>() {
                    @Override
                    void put(ReadMapping obj) {
                        if (obj.mapped && obj.cdr3nt == "TGTGCGAGCGATCGGAACGGTATGGACGTCTGG") {
                            if (!obj.mapping.dSegment.name == "IGHD1-1*01")
                                badDCount.incrementAndGet()
                        } else {
                            wrongMappingCount.incrementAndGet()
                        }
                    }

                    @Override
                    void close() {

                    }
                },
                ReadMappingFilter.createDummy())

        pipeline.run()

        assert badDCount.get() == 0
        assert wrongMappingCount.get() == 0
    }

    @Test
    void badDataTest() {
        def reader = new FastqReader("bad_sample.fastq.gz", true)

        def filter = new ReadMappingFilter((byte) 20, true, true, true, true)

        def pipeline = new Pipeline(reader,
                BlastTestUtil.blastInstanceFactory,
                new ReadMappingOutput(new PlainTextOutput(NullOutputStream.INSTANCE),
                        new ReadMappingDetailsProvider()),
                filter)

        pipeline.run()

        assert filter.passed == filter.good
        assert pipeline.readMappingFilter.mappedRatio >= 0.90
        assert 1 - pipeline.readMappingFilter.noCdr3Ratio >= 0.03
    }
}
