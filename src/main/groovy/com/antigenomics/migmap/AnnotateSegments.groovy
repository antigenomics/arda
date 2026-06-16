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

import com.antigenomics.migmap.blast.BlastInstanceFactory
import com.antigenomics.migmap.io.Read
import com.antigenomics.migmap.pipeline.ExecutionUtil
import com.antigenomics.migmap.pipeline.Util

import java.util.regex.Pattern

def ALLOWED_CHAINS = ["TRA", "TRB", "TRG", "TRG", "IGH", "IGL", "IGK"],
    ALLOWED_SPECIES = ["human", "mouse", "rat", "rabbit", "rhesus_monkey"],
    HOME = new File(this.class.protectionDomain.codeSource.location.path).parent.replaceAll("%20", " ")

def cli = new CliBuilder(usage: "AnnotateSegments input.txt output.txt")

// Runtime
cli._(longOpt: "blast-dir", args: 1, argName: "path",
        "Path to folder that contains 'igblastn' and 'makeblastdb' binaries. " +
                "[default = assume they are added to \$PATH and execute them directly]")
cli._(longOpt: "data-dir", args: 1, argName: "path",
        "Path to folder that contains data bundle (internal_data/ and optional_file/ directories). " +
                "[default = \$install_dir/data/]")
cli.p(args: 1, argName: "int",
        "Number of cores to use. [default = all available processors]")

// Mapping
cli.R(args: 1, argName: "chain1,...",
        "Receptor gene and chain. Several chains can be specified, separated with commas. " +
                "Allowed values: $ALLOWED_CHAINS. [required]")
cli.S(args: 1, argName: "name",
        "Species. Allowed values: $ALLOWED_SPECIES. [required]")
cli._(longOpt: "use-kabat",
        "Will use KABAT nomenclature for CDR/FW partitioning. " +
                "[default = use IMGT nomenclature]")

// Misc
cli.h("Display this help message")

// PARSE ARGUMENTS

def opt = cli.parse(args)

if (opt.h || opt == null || opt.arguments().size() != 2 || !opt.R || !opt.S) {
    cli.usage()
    System.exit(3)
}

// I/O

def inputFileName = opt.arguments()[0], outputFileName = opt.arguments()[1]

if (!new File(inputFileName).exists()) {
    Util.error("Input file $inputFileName does not exist.", 3)
}

// Ensure all subdirs are created for output
def parentFolder = new File(outputFileName).absoluteFile.parentFile
if (parentFolder) {
    parentFolder.mkdirs()
}

// Runtime

if (opt.'blast-dir') {
    ExecutionUtil.BLAST_HOME = new File((String) opt.'blast-dir').absolutePath
}

ExecutionUtil.checkBlastBinaries()

def dataDir = new File((String) (opt.'data-dir' ?: "$HOME/data/")).absolutePath

["internal_data/", "optional_file/"].each {
    if (!new File("$dataDir/$it").exists()) {
        Util.error("Internal data folder $it missing in $dataDir.", 3)
    }
}

// Mapping

def species = (String) opt.S, genes = ((String) opt.R).split(",") as List<String>,
    useKabat = (boolean) opt.'use-kabat'

if (!ALLOWED_SPECIES.contains(species)) {
    Util.error("Unknown species $species.", 3)
}

boolean tr = false, ig = false
genes.each { gene ->
    if (!ALLOWED_CHAINS.contains(gene)) {
        Util.error("Unknown gene $gene.", 3)
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

// Annotate input table

// Deal with V segments first

Util.report("Annotating V segments...", 2)

def vReferencePoints = new HashMap<Integer, Integer>()

def blastInstanceFactory = new BlastInstanceFactory(dataDir, species, genes, true, useKabat, null, null, false)

def instance = blastInstanceFactory.create()

def readThread = new Thread(new Runnable() {
    @Override
    void run() {
        def readMapping
        while ((readMapping = instance.take()) != null) {
            if (readMapping.mapped && readMapping.mapping.regionMarkup.cdr3Start >= 0) {
                vReferencePoints[readMapping.read.header.toInteger()] = readMapping.mapping.regionMarkup.cdr3Start + 3
            }
        }
    }
})

readThread.start()

def i = 0
new File(inputFileName).splitEachLine("[\t ]+") { splitLine ->
    if (!splitLine[0].startsWith("#")) {
        def seq = splitLine[5], segmentTypeStr = splitLine[2]

        if (segmentTypeStr.startsWith("V")) {
            instance.put(new Read(i.toString(), seq))
        }
    }
    i++
}

instance.put(null)

readThread.join()

// Read one more time & assign reference points

Util.report("Annotating J segments and writing output.", 2)

def J_PATTERN = Pattern.compile("T(?:GG|T[TC])GG.{4}GG.")

def getJRef = { String seq ->
    def matcher = J_PATTERN.matcher(seq)
    matcher.find() ? matcher.start() - 1 : -2
}

i = 0
new File(outputFileName).withPrintWriter { pw ->
    new File(inputFileName).splitEachLine("[\t ]+") { splitLine ->
        def refPoint = -1

        if (!splitLine[0].startsWith("#")) {
            def seq = splitLine[5], segmentTypeStr = splitLine[2]

            if (segmentTypeStr.startsWith("V")) {
                refPoint = vReferencePoints[i] ?: -2
            } else if (segmentTypeStr.startsWith("J")) {
                refPoint = getJRef(seq)
            }

            if (refPoint > -2) {
                splitLine[4] = refPoint.toString()
            } else {
                splitLine[0] = "#" + splitLine[0]
            }
        }

        pw.println(splitLine.join("\t"))
        i++
    }
}

Util.report("Finished.", 2)