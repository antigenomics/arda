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

import com.antigenomics.migmap.pipeline.ExecutionUtil
import com.antigenomics.migmap.genomic.SegmentDatabase

class BlastInstanceFactory {
    final List cmdLine
    final File dir
    final List<String> env
    final BlastParser parser
    final String species
    final List<String> genes
    final SegmentDatabase segmentDatabase
    final boolean allAlleles, useKabat, byRead

    BlastInstanceFactory(String dataBundlePath,
                         String species, List<String> genes,
                         boolean allAlleles = true, boolean useKabat = false,
                         String customDatabaseFileName = null,
                         String databaseTempPath = null,
                         boolean byRead = true) {
        this.segmentDatabase = new SegmentDatabase(dataBundlePath, species, genes, allAlleles,
                customDatabaseFileName, databaseTempPath)
        this.parser = new BlastParser(segmentDatabase)
        this.species = species
        this.genes = genes
        this.allAlleles = allAlleles
        this.useKabat = useKabat
        this.byRead = byRead

        def seqtype = genes[0].startsWith("TR") ? "TCR" : "Ig"

        def OPTS = ["SEGM_OPT"  :
                            ["v", "d", "j"].collect { segment ->
                                "-germline_db_${segment.toUpperCase()} ${new File(segmentDatabase.databaseTempPath).absolutePath}/$segment"
                            },

                    "AUX_OPT"   :
                            "-auxiliary_data $dataBundlePath/optional_file/${species}_gl.aux",

                    "ORG_OPT"   :
                            "-organism $species -ig_seqtype $seqtype",

                    "DOMAIN_OPT":
                            "-domain_system ${useKabat ? "kabat" : "imgt"}"
        ]

        def REPORT_OPT = ["-num_alignments_V 3", "-num_alignments_D 3", "-num_alignments_J 3"],
            OUTFMT_OPT = ["-outfmt", "7 sseqid qstart qseq sstart sseq"],
            OUTPUT_OPT = "-out -"

        this.cmdLine = [[ExecutionUtil.igBlast, OPTS.values(), REPORT_OPT, OUTPUT_OPT].
                                flatten().join(" ").split(" "),
                        OUTFMT_OPT
        ].flatten()

        this.env = ["IGDATA=\"$dataBundlePath\""]
        this.dir = new File(dataBundlePath) // <- this is most crucial here. Blast will not find internal_data otherwise

        segmentDatabase.makeBlastDb()
    }

    void annotateV() {
        segmentDatabase.annotateV(this)
    }

    BlastInstance create() {
        new BlastInstance(cmdLine.execute(env, dir), parser, segmentDatabase, byRead)
    }
}
