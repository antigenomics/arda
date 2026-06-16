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

package com.antigenomics.migmap.post

import com.antigenomics.migmap.mutation.Mutation
import com.milaboratory.core.sequence.NucleotideSequence
import groovy.transform.Canonical

@Canonical
class CorrectorClonotypeEntry {
    int count
    float freq
    List<String> data
    NucleotideSequence seq
    CorrectorClonotypeEntry parent
    Set<Mutation> mutations = new HashSet<>()

    void append(CorrectorClonotypeEntry other) {
        count += other.count
        freq += other.freq
    }
}
