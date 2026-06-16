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

package com.antigenomics.migmap.mapping

import com.antigenomics.migmap.io.Read
import groovy.transform.CompileStatic

import static com.antigenomics.migmap.pipeline.Util.translateLinear
import static java.lang.Math.max

class ReadMappingDetailsProvider {
    public static
    final List<String> ALLOWED_FIELDS = Collections.unmodifiableList(["fr1nt", "cdr1nt",
                                                                      "fr2nt", "cdr2nt",
                                                                      "fr3nt", "fr4nt",
                                                                      "contignt",
                                                                      "fr1aa", "cdr1aa",
                                                                      "fr2aa", "cdr2aa",
                                                                      "fr3aa", "fr4aa",
                                                                      "contigaa"])

    public static ReadMappingDetailsProvider DUMMY = new ReadMappingDetailsProvider([])

    private final List<String> fields

    ReadMappingDetailsProvider() {
        this(ALLOWED_FIELDS)
    }

    ReadMappingDetailsProvider(List<String> fields) {
        this.fields = fields.collect { it.toLowerCase() }
        def badFields = fields.findAll { !ALLOWED_FIELDS.contains(it) }
        if (!badFields.empty) {
            throw new RuntimeException("Bad fields supplied: ${badFields.join(",")}, " +
                    "allowed values: ${ALLOWED_FIELDS.join(",")}.")
        }
    }

    @CompileStatic
    String getHeader() {
        fields.empty ? "" : ("\t" + fields.join("\t"))
    }

    String getDetailsString(ReadMapping readMapping) {
        if (fields.empty) {
            return ""
        }

        def details = getDetails(readMapping)

        "\t" + fields.collect { details."$it" }.join("\t")
    }

    @CompileStatic
    static ReadMappingDetails getDetails(ReadMapping readMapping) {
        readMapping.mapped ?
                new ReadMappingDetails(readMapping.read, readMapping.mapping) :
                ReadMappingDetails.DUMMY
    }

    @CompileStatic
    static class ReadMappingDetails {
        final String seq
        final RegionMarkup readMarkup, referenceMarkup
        final int vStartInRef, vStartInQuery, cdr3Start, cdr3End

        static final ReadMappingDetails DUMMY = new ReadMappingDetails()

        ReadMappingDetails(Read read, Mapping mapping) {
            this.seq = read.seq
            this.vStartInRef = mapping.vStartInRef
            this.vStartInQuery = mapping.vStartInQuery
            this.referenceMarkup = mapping.vSegment.regionMarkup

            if (referenceMarkup == null) {
                throw new RuntimeException("No reference markup, " +
                        "it seems you're running on an unnatotated segment database, " +
                        "forgot factory.annotateV()? :)")
            }

            this.readMarkup = mapping.regionMarkup
            this.cdr3Start = max(readMarkup.cdr2End, readMarkup.cdr3Start)
            this.cdr3End = readMarkup.cdr3End
        }

        ReadMappingDetails() {
            this.seq = null
            this.vStartInRef = 480011
            this.vStartInQuery = vStartInRef
            this.cdr3Start = -1
            this.cdr3End = -1
            this.referenceMarkup = RegionMarkup.DUMMY
            this.readMarkup = RegionMarkup.DUMMY
        }

        int nCount(int pos) {
            max(vStartInRef - pos, 0)
        }

        int sCount(int pos) {
            max(vStartInQuery, pos)
        }

        String getFr1nt() {
            vStartInRef < referenceMarkup.cdr1Start && vStartInQuery < readMarkup.cdr1Start ?
                    'N' * nCount(0) + seq.substring(sCount(0), readMarkup.cdr1Start) :
                    'N' * max(0, referenceMarkup.cdr1Start)
        }

        String getCdr1nt() {
            vStartInRef < referenceMarkup.cdr1End && vStartInQuery < readMarkup.cdr1End ?
                    'N' * nCount(referenceMarkup.cdr1Start) + seq.substring(sCount(readMarkup.cdr1Start), readMarkup.cdr1End) :
                    'N' * (referenceMarkup.cdr1End - referenceMarkup.cdr1Start)
        }

        String getFr2nt() {
            vStartInRef < referenceMarkup.cdr2Start && vStartInQuery < readMarkup.cdr2Start ?
                    'N' * nCount(referenceMarkup.cdr1End) + seq.substring(sCount(readMarkup.cdr1End), readMarkup.cdr2Start) :
                    'N' * (referenceMarkup.cdr2Start - referenceMarkup.cdr1End)
        }

        String getCdr2nt() {
            vStartInRef < referenceMarkup.cdr2End && vStartInQuery < readMarkup.cdr2End ?
                    'N' * nCount(referenceMarkup.cdr2Start) + seq.substring(sCount(readMarkup.cdr2Start), readMarkup.cdr2End) :
                    'N' * (referenceMarkup.cdr2End - referenceMarkup.cdr2Start)
        }

        String getFr3nt() {
            vStartInRef < referenceMarkup.cdr3Start && vStartInQuery < readMarkup.cdr3Start ?
                    'N' * nCount(referenceMarkup.cdr2End) + seq.substring(sCount(readMarkup.cdr2End), cdr3Start) :
                    'N' * (referenceMarkup.cdr3Start - referenceMarkup.cdr2End)
        }

        String getFr4nt() {
            cdr3End >= 0 ? seq.substring(cdr3End) : ""
        }

        String getContignt() {
            fr1nt + cdr1nt + fr2nt + cdr2nt + fr3nt + downstream
        }

        String getFr1aa() {
            translateLinear(fr1nt)
        }

        String getCdr1aa() {
            translateLinear(cdr1nt)
        }

        String getFr2aa() {
            translateLinear(fr2nt)
        }

        String getCdr2aa() {
            translateLinear(cdr2nt)
        }

        String getFr3aa() {
            translateLinear(fr3nt)
        }

        String getFr4aa() {
            translateLinear(fr4nt)
        }

        String getContigaa() {
            translateLinear(contignt)
        }

        String getDownstream() {
            cdr3Start >= 0 ? seq.substring(cdr3Start) : ""
        }
    }
}