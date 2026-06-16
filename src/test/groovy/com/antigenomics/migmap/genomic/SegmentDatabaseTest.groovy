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

package com.antigenomics.migmap.genomic

import com.antigenomics.migmap.blast.BlastInstanceFactory
import org.junit.Test

class SegmentDatabaseTest {
    @Test
    void massiveLoadTest() {
        ["human", "mouse", "rat", "rabbit", "rhesus_monkey"].each {
            def segmentDatabase = new SegmentDatabase("data/", it, ["TRA", "TRB", "TRG", "TRD",
                                                                    "IGH", "IGK", "IGL"])

            assert !segmentDatabase.segments.isEmpty()
            assert segmentDatabase.vSegments > 0
            assert segmentDatabase.dSegments > 0
            assert segmentDatabase.jSegments > 0
        }
    }

    @Test
    void loadTest() {
        def segmentDatabase = new SegmentDatabase("data/", "human", ["IGH"])

        assert !segmentDatabase.segments.isEmpty()

        def testSegments = ["IGHV4-59*09",
                            "IGHV3-23*01",
                            "IGHV3/OR16-8*01",
                            "IGHV3-73*02",
                            "IGHV4-39*06",
                            "IGHV3-30-3*03",
                            "IGHV1/OR15-2*01",
                            "IGHV4/OR15-8*02",
                            "IGHV4-28*05",
                            "IGHV1/OR15-3*02"]

        assert testSegments.each { segmentDatabase.segments.keySet().contains(it) }
    }

    @Test
    void makeDbTest() {
        def segmentDatabase = new SegmentDatabase("data/", "human", ["TRA"])
        segmentDatabase.makeBlastDb()

        assert new File(segmentDatabase.databaseTempPath).exists()
        assert new File("$segmentDatabase.databaseTempPath/v.nhr").exists()
        assert new File("$segmentDatabase.databaseTempPath/d.nhr").exists()
        assert new File("$segmentDatabase.databaseTempPath/j.nhr").exists()
    }

    @Test
    void annotationTestIG() {
        def factory = new BlastInstanceFactory("data/", "human", ["IGH", "IGK", "IGL"])
        factory.annotateV()

        assert factory.segmentDatabase.annotatedV == factory.segmentDatabase.vSegments
    }

    @Test
    void annotationTestTR() {
        def factory = new BlastInstanceFactory("data/", "human", ["TRA", "TRB", "TRG", "TRD"])
        factory.annotateV()

        assert factory.segmentDatabase.annotatedV + 2 == factory.segmentDatabase.vSegments // known 2 unannotateable cases
    }
}
