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

import com.antigenomics.migmap.pipeline.Util
import com.antigenomics.migmap.genomic.SegmentDatabase
import com.antigenomics.migmap.io.InputPort
import com.antigenomics.migmap.io.OutputPort
import com.antigenomics.migmap.io.Read
import com.antigenomics.migmap.mapping.Mapping
import com.antigenomics.migmap.mapping.ReadMapping
import com.antigenomics.migmap.mutation.Mutation
import com.antigenomics.migmap.mutation.MutationConverter
import com.antigenomics.migmap.mutation.MutationType
import groovy.transform.CompileStatic

import java.util.regex.Pattern

class BlastInstance implements OutputPort<ReadMapping>, InputPort<Read> {
    final Queue<Read> reads = new LinkedList<>() // todo: optimize?
    final Process proc
    final BufferedReader reader
    final PrintWriter writer
    final BlastParser parser
    private final DSearcherBundle auxDSearcher
    final boolean byRead

    protected boolean last = false

    protected BlastInstance(Process proc, BlastParser parser, SegmentDatabase segmentDatabase,
                            boolean byRead) {
        this.proc = proc
        this.reader = proc.in.newReader()
        this.writer = proc.out.newPrintWriter()
        this.parser = parser
        this.auxDSearcher = new DSearcherBundle(segmentDatabase)
        this.byRead = byRead
    }

    static String getHeader(String chunk) {
        (chunk =~ /# Query: (.+)/)[0][1]
    }

    @CompileStatic
    String nextChunk() {
        if (last) {
            // finished
            return null
        }

        String chunk = "", line

        while (true) {
            line = reader.readLine()

            if (line == null) {
                if (chunk.length() == 0) {
                    return null // empty input
                }
                last = true // warn of EOF
                break // went to EOF
            } else if (line.startsWith("# IGBLASTN") && chunk.length() > 0) {
                break // went to the beginning of next chunk
            } else {
                chunk += line + "\n"
            }
        }

        chunk
    }

    ReadMapping createReadMapping(Mapping mapping, Read read) {
        String cdr3nt, cdr3aa
        byte[] mutationQual, cdrInsertQual
        boolean canonical, inFrame, noStop
        PSegments pSegments

        if (mapping) {
            if (mapping.rc) {
                // don't forget to reverse complement
                read = read.rc
            }

            def seq = read.seq
            def regionMarkup = mapping.regionMarkup

            mutationQual = new byte[mapping.mutations.size()]

            mapping.mutations.eachWithIndex { Mutation it, Integer i ->
                mutationQual[i] = it.type == MutationType.Substitution ? read.qualAt(it.posInRead) : Util.MAX_QUAL
            }

            if (mapping.hasCdr3) {
                if (mapping.complete) {
                    cdr3nt = seq.substring(regionMarkup.cdr3Start, regionMarkup.cdr3End)
                    cdr3aa = Util.translateCdr(cdr3nt)

                    // Try refining D segment assignment
                    /*if (mapping.hasD && !mapping.dFound) {
                        // search only small inserts, IgBlast handles the rest OK
                        def searchResult = auxDSearcher.search(mapping.vSegment.gene,
                                mapping.cdr3Markup, cdr3nt)

                        if (searchResult != null &&
                                searchResult.score <= 0.2) {
                            // require match of considerable size 2 of 4, 3 of 14, ...
                            mapping = searchResult.updateMapping(mapping)
                        }
                    }*/

                    def cdrMarkup = mapping.cdr3Markup

                    // Check for P-segments
                    pSegments = PSegmentSearcher.search(mapping.cdr3Markup, mapping.truncations, cdr3nt)

                    // Quality of N-nucleotides
                    int i = 0
                    if (mapping.dFound) {
                        cdrInsertQual = new byte[Math.max(0, cdrMarkup.dStart - cdrMarkup.vEnd) +
                                Math.max(0, cdrMarkup.jStart - cdrMarkup.dEnd)]

                        if (cdrMarkup.vEnd <= cdrMarkup.dStart) {
                            ((regionMarkup.cdr3Start + cdrMarkup.vEnd)..<(regionMarkup.cdr3Start + cdrMarkup.dStart)).each {
                                cdrInsertQual[i++] = read.qualAt(it)
                            }
                        }

                        if (cdrMarkup.dEnd <= cdrMarkup.jStart) {
                            ((regionMarkup.cdr3Start + cdrMarkup.dEnd)..<(regionMarkup.cdr3Start + cdrMarkup.jStart)).each {
                                cdrInsertQual[i++] = read.qualAt(it)
                            }
                        }
                    } else {
                        cdrInsertQual = new byte[Math.max(0, cdrMarkup.jStart - cdrMarkup.vEnd)]
                        if (cdrMarkup.vEnd <= cdrMarkup.jStart) {
                            ((regionMarkup.cdr3Start + cdrMarkup.vEnd)..<(regionMarkup.cdr3Start + cdrMarkup.jStart)).each {
                                cdrInsertQual[i++] = read.qualAt(it)
                            }
                        }
                    }
                    canonical = Util.isCanonical(cdr3nt)
                } else {
                    cdr3nt = seq.substring(regionMarkup.cdr3Start)
                    cdr3aa = Util.translateLinear(cdr3nt)
                    cdrInsertQual = new byte[cdr3nt.length()]
                    (0..<cdr3nt.length()).each {
                        cdrInsertQual[it] = read.qualAt(regionMarkup.cdr3Start + it)
                    }
                    canonical = false
                    pSegments = new PSegments(-1, -1, -1, -1)
                }

                inFrame = mapping.inFrame && !cdr3aa.contains("?")
                noStop = mapping.noStop && !cdr3aa.contains("*")
            } else {
                cdr3nt = Util.MY_NA
                cdr3aa = Util.MY_NA
                mutationQual = new byte[0]
                cdrInsertQual = new byte[0]
                canonical = false
                inFrame = mapping.inFrame
                noStop = mapping.noStop
                pSegments = new PSegments(-1, -1, -1, -1)
            }
        } else {
            cdr3nt = null
            cdr3aa = null
            mutationQual = null
            cdrInsertQual = null
            canonical = false
            inFrame = false
            noStop = false
            pSegments = null
        }

        ReadMapping readMapping = new ReadMapping(read, mapping,
                cdr3nt, cdr3aa,
                mutationQual, cdrInsertQual,
                canonical, inFrame, noStop,
                pSegments)

        if (byRead) {
            MutationConverter.annotateMutationAa(readMapping)
        }

        readMapping
    }

    @Override
    @CompileStatic
    ReadMapping take() {
        def chunk = nextChunk(), read

        if (chunk && (read = reads.poll())) { // second condition protects from empty input
            return createReadMapping(parser.parse(chunk), read)
        } else {
            close()
            return null
        }
    }

    static final Pattern ALLOWED_BASES = ~/^[ATGCUatgcuRrYySsWwKkMmBbDdHhVvNn]+$/,
                         AMBIGUOUS_BASE = ~/[^ATGCUatgcu]/;

    @CompileStatic
    static boolean isGoodBlastSeq(String seq) {
        if (seq.length() == 0 || !seq =~ ALLOWED_BASES) {
            return false
        }

        int goodBases = seq.length() - seq.replaceAll(AMBIGUOUS_BASE, "").length()

        if ((goodBases / (double) seq.length()) >= 0.7) {
            return false
        }

        true
    }

    @Override
    @CompileStatic
    void put(Read input) {
        if (input) { // should contain at least one non-ambiguous base, otherwise BLAST crashes
            if (isGoodBlastSeq(input.seq)) {
                reads.add(input)
                writer.println(">" + input.header + "\n" + input.seq)
            }
        } else {
            // no more reads
            writer.close()
        }
    }

    @Override
    void close() {
        proc.waitFor()

        if (proc.exitValue() > 0) {
            Util.error("IgBlast failed with CODE${proc.exitValue()}:\n${proc.getErrorStream()}", 2)
        }
    }
}
