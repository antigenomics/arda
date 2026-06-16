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

package com.antigenomics.migmap.blast

import com.antigenomics.migmap.genomic.Segment
import com.antigenomics.migmap.genomic.SegmentDatabase
import com.antigenomics.migmap.mapping.Cdr3Markup
import com.antigenomics.migmap.mapping.Mapping
import com.antigenomics.migmap.mapping.RegionMarkup
import com.antigenomics.migmap.mapping.Truncations
import com.antigenomics.migmap.mutation.MutationExtractor
import groovy.transform.CompileStatic

import java.util.concurrent.atomic.AtomicInteger
import java.util.regex.Pattern

import static com.antigenomics.migmap.pipeline.Util.BLAST_NA
import static com.antigenomics.migmap.pipeline.Util.groomMatch

@CompileStatic
class BlastParser {
    final static int MIN_CDR3_LEN = 6
    final AtomicInteger total = new AtomicInteger(),
                        noMatch = new AtomicInteger(),
                        vNotFound = new AtomicInteger(),
                        noCdr3 = new AtomicInteger()
    final Pattern dSummaryPattern =
            //                 V     D     J  chain stop frame (prod) strand
            Pattern.compile(/# V-.+\n(\S+)\t(\S+)\t(\S+)\tV\S\t(.+)\t(\S+)\t\S+\t(\S+)\n/),
                  noDSummaryPattern =
                          //                  V     J    chain   stop frame (prod)  strand
                          Pattern.compile(/# V-.+\n(\S+)\t(\S+)\tV\S\t(\S+)\t(\S+)\t\S+\t(\S+)\n/),
                  cdr1Pattern = Pattern.compile(/# Alignment summary(?:.+\n)+CDR1-IMGT\t([0-9]+)\t([0-9]+)/),
                  cdr2Pattern = Pattern.compile(/# Alignment summary(?:.+\n)+CDR2-IMGT\t([0-9]+)\t([0-9]+)/),
                  cdr3Pattern = Pattern.compile(/# Alignment summary(?:.+\n)+CDR3-IMGT \(germline\)\t([0-9]+)\t([0-9]+)/)


    final SegmentDatabase segmentDatabase

    public BlastParser(SegmentDatabase segmentDatabase) {
        this.segmentDatabase = segmentDatabase
    }

    static Alignment createAlignment(List<String> hitChunk) {
        new Alignment(hitChunk[0].toInteger() - 1, hitChunk[1], hitChunk[2].toInteger() - 1, hitChunk[3])
    }

    Mapping createMapping(List<String> summary,
                          List<String> cdr1Bounds, List<String> cdr2Bounds, List<String> cdr3Bounds,
                          Segment vSegment, Segment dSegment, Segment jSegment,
                          Alignment vAlignment, Alignment jAlignment, Alignment dAlignment,
                          boolean hasD) {
        def rc = summary[-1] != "+", inFrame = summary[-2] != "Out-of-frame", noStop = summary[-3] != "Yes"

        int cdr1Start = -1, cdr1End = -1,
            cdr2Start = -1, cdr2End = -1,
            cdr3Start = -1, cdr3End = -1

        if (cdr1Bounds) {
            cdr1Start = cdr1Bounds[0].toInteger() - 1
            cdr1End = cdr1Bounds[1].toInteger()
        }

        if (cdr2Bounds) {
            cdr2Start = cdr2Bounds[0].toInteger() - 1
            cdr2End = cdr2Bounds[1].toInteger()
        }

        // - Find CDR3 end using J reference point manually
        if (cdr3Bounds) {
            cdr3Start = cdr3Bounds[0].toInteger() - 4
        } else {
            // try rescue CDR3
            cdr3Start = RefPointSearcher.getCdr3Start(vSegment, vAlignment)
        }

        if (jSegment != Segment.DUMMY_J && cdr3Start >= 0) {
            cdr3End = RefPointSearcher.getCdr3End(jSegment, jAlignment)
        }

        def regionMarkup = new RegionMarkup(cdr1Start, cdr1End, cdr2Start, cdr2End, cdr3Start, cdr3End)
        def cdr3Len = cdr3End - cdr3Start
        def hasCdr3 = cdr3Start >= 0 && (cdr3End < 0 || cdr3Len >= MIN_CDR3_LEN),
            complete = cdr3End >= 0 && hasCdr3

        // - Markup of V/D/J within CDR3
        int vCdr3End = -1, dCdr3Start = -1, dCdr3End = -1, jCdr3Start = -1,
            vDel = -1, dDel5 = -1, dDel3 = -1, jDel = -1

        if (hasCdr3) {
            vCdr3End = vAlignment.qend - cdr3Start
            vDel = vSegment.sequence.length() - vAlignment.send

            if (jSegment != Segment.DUMMY_J) {
                jCdr3Start = jAlignment.qstart - cdr3Start
                jDel = jAlignment.sstart
            }

            if (dSegment != Segment.DUMMY_D) {
                dCdr3Start = dAlignment.qstart - cdr3Start
                dCdr3End = dCdr3Start + dAlignment.qLength
                dDel5 = dAlignment.sstart
                dDel3 = dSegment.sequence.length() - dAlignment.send
            }

            if (dCdr3Start < 0 || dCdr3End >= cdr3Len) {
                // Handle rare D out of CDR3 cases
                dSegment = Segment.DUMMY_D
                dCdr3Start = -1
                dCdr3End = -1
                dDel5 = -1
                dDel3 = -1
            }
        } else {
            noCdr3.incrementAndGet()
        }

        Cdr3Markup cdr3Markup = new Cdr3Markup(vCdr3End, dCdr3Start, dCdr3End, jCdr3Start)
        def truncations = new Truncations(vDel, dDel5, dDel3, jDel)

        // Finally, deal with hypermutations
        // offset for converting coordinate in read to coordinate in germline V
        def mutationExtractor = new MutationExtractor(vSegment, vAlignment, regionMarkup)

        if (hasCdr3) {
            if (dSegment != Segment.DUMMY_D) {
                mutationExtractor.extractD(dSegment, dAlignment)
            }
            if (jSegment != Segment.DUMMY_J) {
                mutationExtractor.extractJ(jSegment, jAlignment)
            }
        }

        return new Mapping(vSegment, dSegment, jSegment, vAlignment.sstart, vAlignment.qstart,
                regionMarkup, cdr3Markup, truncations,
                rc, complete, hasCdr3, inFrame, noStop, hasD, dSegment != Segment.DUMMY_D,
                mutationExtractor.mutations)
    }

    Mapping parse(String chunk) {
        total.incrementAndGet()

        // Rearrangement summary
        boolean hasD = true // tells if the chain has D segment, while dFound tells whether it was found

        def summary = groomMatch(chunk =~ dSummaryPattern)

        if (summary == null) {
            hasD = false
            summary = groomMatch(chunk =~ noDSummaryPattern)
        }

        if (summary == null) {
            noMatch.incrementAndGet()
            return null
        }

        // Information on segments mapped
        // - Segment names, can be multiple of them
        def vSegmentNames = summary[0],
            jSegmentNames = hasD ? summary[2] : summary[1],
            dSegmentNames = hasD ? summary[1] : BLAST_NA

        boolean dFound = dSegmentNames != BLAST_NA,
                vFound = vSegmentNames != BLAST_NA,
                jFound = jSegmentNames != BLAST_NA

        if (!vFound) {
            vNotFound.incrementAndGet()
            return null
        }

        Segment vSegment = segmentDatabase.segments[vSegmentNames.split(",").sort()[0]],
                dSegment = dFound ? segmentDatabase.segments[dSegmentNames.split(",").sort()[0]] : Segment.DUMMY_D,
                jSegment = jFound ? segmentDatabase.segments[jSegmentNames.split(",").sort()[0]] : Segment.DUMMY_J

        // - Alignments for V, D and J segments, remember here and further BLAST coordinates are 1-based

        Alignment vAlignment = createAlignment(groomMatch(chunk =~
                //                                            qstart     qseq        sstart     sseq
                /# Hit table(?:.+\n)+V\t$vSegment.regexName\t([0-9]+)\t([ATGCN-]+)\t([0-9]+)\t([ATGCN-]+)/)),
                  dAlignment = dFound ? createAlignment(groomMatch(chunk =~
                          /# Hit table(?:.+\n)+D\t$dSegment.regexName\t([0-9]+)\t([ATGCN-]+)\t([0-9]+)\t([ATGCN-]+)/)) : null,
                  jAlignment = jFound ? createAlignment(groomMatch(chunk =~
                          /# Hit table(?:.+\n)+J\t$jSegment.regexName\t([0-9]+)\t([ATGCN-]+)\t([0-9]+)\t([ATGCN-]+)/)) : null

        // Deduce CDR/FW markup
        // - CDR1,2 and CDR3(start only) coords
        def cdr1Bounds = groomMatch(chunk =~ cdr1Pattern),
            cdr2Bounds = groomMatch(chunk =~ cdr2Pattern),
            cdr3Bounds = groomMatch(chunk =~ cdr3Pattern)

        createMapping(summary,
                cdr1Bounds, cdr2Bounds, cdr3Bounds,
                vSegment, dSegment, jSegment,
                vAlignment, jAlignment, dAlignment,
                hasD)
    }
}