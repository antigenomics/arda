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

package com.antigenomics.migmap.mutation

import groovy.transform.CompileStatic

@CompileStatic
class MutationFormatter {
    static
    final String OUTPUT_HEADER_NT = SubRegion.REGION_LIST.collect { SubRegion it -> "mutations.nt." + it }.join("\t"),
                 OUTPUT_HEADER_AA = SubRegion.REGION_LIST.collect { SubRegion it -> "mutations.aa." + it }.join("\t")

    static String toStringNT(List<Mutation> mutations) {
        // assume mutations are sorted

        def mutationStrings = [""] * SubRegion.REGION_LIST.length

        mutations.each {
            int order = it.subRegion.order

            if (mutationStrings[order].length() > 0)
                mutationStrings[order] += ","

            mutationStrings[order] += it.toString()
        }

        mutationStrings.join("\t")
    }

    static String toStringAA(List<Mutation> mutations) {
        if (mutations.any { it.aaFrom == null })
            return ""

        def mutationStrings = [""] * SubRegion.REGION_LIST.length

        mutations.each {
            int order = it.subRegion.order

            if (mutationStrings[order].length() > 0)
                mutationStrings[order] += ","

            mutationStrings[order] += it.toStringAa()
        }

        mutationStrings.join("\t")
    }
}
