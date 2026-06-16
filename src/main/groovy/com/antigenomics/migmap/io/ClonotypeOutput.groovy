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
import com.antigenomics.migmap.clonotype.Clonotype
import com.antigenomics.migmap.clonotype.ClonotypeAccumulator
import com.antigenomics.migmap.mapping.ReadMapping
import com.antigenomics.migmap.mapping.ReadMappingDetailsProvider
import groovy.transform.CompileStatic

@CompileStatic
class ClonotypeOutput implements InputPort<ReadMapping> {
    final ReadMappingDetailsProvider readMappingDetailsProvider
    final ClonotypeAccumulator clonotypeAccumulator
    final PlainTextOutput plainTextOutput

    ClonotypeOutput(PlainTextOutput plainTextOutput = StdOutput.INSTANCE,
                    ReadMappingDetailsProvider readMappingDetailsProvider = ReadMappingDetailsProvider.DUMMY) {
        this.plainTextOutput = plainTextOutput
        this.clonotypeAccumulator = new ClonotypeAccumulator()
        this.readMappingDetailsProvider = readMappingDetailsProvider
        if (plainTextOutput != StdOutput.INSTANCE)
            plainTextOutput.put(Clonotype.OUTPUT_HEADER + readMappingDetailsProvider.header)
    }

    @Override
    void put(ReadMapping readMapping) {
        clonotypeAccumulator.put(readMapping)
    }

    @Override
    void close() {
        Util.report("Sorting ${clonotypeAccumulator.clonotypeMap.size()} clonotype entries, " +
                "writing output.", 2)

        List<Clonotype> clonotypes = clonotypeAccumulator.clonotypeMap.collect {
            new Clonotype(it.key, it.value, clonotypeAccumulator.total)
        }.sort()

        clonotypes.each {
            plainTextOutput.put(it.toString() + readMappingDetailsProvider.getDetailsString(it.representativeMapping))
        }

        plainTextOutput.close()

        Util.report("Finished.", 2)
    }
}
