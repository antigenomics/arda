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

import com.antigenomics.migmap.mutation.SubRegion
import com.antigenomics.migmap.pipeline.PipelineTestCache
import com.antigenomics.migmap.clonotype.Clonotype
import com.antigenomics.migmap.mutation.Mutation
import org.junit.Test

class ClonotypeLoaderTest {
    @Test
    void decodeMutationTest() {
        def ntMut = "S5:G>A,S8:G>C,S41:T>G,S67:C>A,S70:C>T,S71:C>T"
        def aaMut = "S1:V>V,S2:Q>H,S13:P>P,S22:A>E,S23:A>V,S23:A>V"

        def mutations = ClonotypeLoader.decodeMutations(ntMut, aaMut, SubRegion.FR1, null, null, null, null)

        assert ntMut == mutations.collect { it.toString() }.join(",")
        assert aaMut == mutations.collect { it.toStringAa() }.join(",")
    }

    @Test
    void decodeMutationSingleTest() {
        def ntMut = "S5:G>A"
        def aaMut = "S1:V>V"

        def mutations = ClonotypeLoader.decodeMutations(ntMut, aaMut, SubRegion.FR1, null, null, null, null)

        assert ntMut == mutations.collect { it.toString() }.join(",")
        assert aaMut == mutations.collect { it.toStringAa() }.join(",")
    }

    @Test
    void decodeMutationEmptyTest() {
        def ntMut = ""
        def aaMut = ""

        def mutations = ClonotypeLoader.decodeMutations(ntMut, aaMut, SubRegion.FR1, null, null, null, null)

        assert ntMut == mutations.collect { it.toString() }.join(",")
        assert aaMut == mutations.collect { it.toStringAa() }.join(",")
    }

    @Test
    void loadSavedTest() {
        def originalClonotypes = PipelineTestCache.INSTANCE.clonotypes

        def loadedClonotypes = ClonotypeLoader.load(PipelineTestCache.INSTANCE.clonotypeOutputFile,
                PipelineTestCache.INSTANCE.segmentDatabase)

        assert originalClonotypes.size() == loadedClonotypes.size()

        loadedClonotypes.eachWithIndex { Clonotype loadedClonotype, int i ->
            def originalClonotype = originalClonotypes[i]
            assert originalClonotype.cdr3nt == loadedClonotype.cdr3nt
            assert originalClonotype.vSegment.name == loadedClonotype.vSegment.name
            assert originalClonotype.dSegment.name == loadedClonotype.dSegment.name
            assert originalClonotype.jSegment.name == loadedClonotype.jSegment.name

            assert originalClonotype.mutations.size() == loadedClonotype.mutations.size()

            loadedClonotype.mutations.eachWithIndex { Mutation loadedMutation, int j ->
                def originalMutation = originalClonotype.mutations[j]
                assert originalMutation.toString() == loadedMutation.toString()
                assert originalMutation.toStringAa() == loadedMutation.toStringAa()
                assert originalMutation.subRegion == loadedMutation.subRegion
            }

            assert originalClonotype.cdr3Markup.toString() == loadedClonotype.cdr3Markup.toString()
            assert originalClonotype.truncations.toString() == loadedClonotype.truncations.toString()
            assert originalClonotype.pSegments.toString() == loadedClonotype.pSegments.toString()
        }
    }
}
