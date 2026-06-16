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

package com.antigenomics.migmap

import com.antigenomics.migmap.genomic.SegmentDatabase
import com.antigenomics.migmap.io.ClonotypeLoader
import com.antigenomics.migmap.pipeline.Util
import com.antigenomics.migmap.post.analysis.Analysis

def ALLOWED_SPECIES = [SegmentDatabase.SPECIES_ALIAS.keySet(), SegmentDatabase.SPECIES_ALIAS.values()].flatten(),
    ALLOWED_CHAINS = ["TRA", "TRB", "TRG", "TRD", "IGH", "IGL", "IGK"]

def cli = new CliBuilder(usage: "Analyze [options] input.txt output_prefix")

cli.R(args: 1, argName: "chain1,...",
        "Receptor gene and chain. Several chains can be specified, separated with commas. " +
                "Allowed values: $ALLOWED_CHAINS. [required]")
cli.S(args: 1, argName: "name",
        "Species. Allowed values: $ALLOWED_SPECIES. [required]")

def opt = cli.parse(args)

if (opt.h || opt == null || opt.arguments().size() != 2) {
    cli.usage()
    System.exit(3)
}

def species = (String) opt.S, genes = ((String) opt.R).split(",") as List<String>

if (!ALLOWED_SPECIES.contains(species)) {
    Util.error("Unknown species $species.", 3)
}

boolean tr = false, ig = false
genes.each { gene ->
    if (!ALLOWED_CHAINS.contains(gene)) {
        Util.error("Unknown gene $gene parameter.", 3)
    }
    if (gene.toUpperCase().startsWith("TR")) {
        tr = true
    }
    if (gene.toUpperCase().startsWith("IG")) {
        ig = true
    }
}

if (tr && ig) {
    Util.error("Mixing TCR and IG genes is not allowed.", 3)
}

//

def segmentDatabase = new SegmentDatabase(".", species, genes)

def inputFile = opt.arguments()[0], outputPrefix = opt.arguments()[1]

def clonotypes = ClonotypeLoader.load(new File(inputFile), segmentDatabase)

def postAnalysis = new Analysis(clonotypes)

postAnalysis.generateHypermutationTable(outputPrefix + ".shm.txt")
postAnalysis.generateClonotypeTree(outputPrefix)