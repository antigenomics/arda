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

import com.antigenomics.migmap.pipeline.Util
import com.antigenomics.migmap.mutation.Mutation
import com.antigenomics.migmap.post.CorrectorClonotypeEntry
import com.milaboratory.core.sequence.NucleotideSequence
import com.milaboratory.core.tree.SequenceTreeMap
import com.milaboratory.core.tree.TreeSearchParameters
import com.milaboratory.util.Factory

def DEFAULT_ERROR_RATE = "0.01", DEFAULT_DEPTH = "1"

def cli = new CliBuilder(usage: "Correct [options] clonotype_table.txt output.txt")

cli.r(args: 1, argName: "float", "Maximum child-to-parent ratio for a single error. [default=$DEFAULT_ERROR_RATE]")
cli.d(args: 1, argName: "int", "Scan depth, maximum number of allowed mismatches. [default=$DEFAULT_DEPTH]")
cli._(longOpt: "contig", "Use entire contigs for error correction. " +
        "May not work properly for non-full-length sequencing and does not rely on parsimony principle.")

// Misc
cli.h("Display this help message")

// PARSE ARGUMENTS

def opt = cli.parse(args)

if (opt.h || opt == null || opt.arguments().size() != 2) {
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

def errorRate = Math.log((opt.r ?: DEFAULT_ERROR_RATE).toFloat()), depth = (opt.d ?: DEFAULT_DEPTH).toInteger(),
    contig = (boolean) opt.'contig'

// Parsing required columns, here we operate on full contigs

def countColIndex = -1, freqColIndex = -1, contigNtIndex = -1,
    cdr3NtIndex = -1,
    mutationsFr1Index = -1, mutationsCdr1Index = -1,
    mutationsFr2Index = -1, mutationsCdr2Index = -1,
    mutationsFr3Index = -1, mutationsFr4Index = -1

def header
def parseColumns = {
    freqColIndex = header.findIndexOf { it.toLowerCase() == "freq" }
    countColIndex = header.findIndexOf { it.toLowerCase() == "count" }
    contigNtIndex = header.findIndexOf { it.toLowerCase() == "contignt" }
    cdr3NtIndex = header.findIndexOf { it.toLowerCase() == "cdr3nt" }
    mutationsFr1Index = header.findIndexOf { it.toLowerCase() == "mutations.nt.fr1" }
    mutationsCdr1Index = header.findIndexOf { it.toLowerCase() == "mutations.nt.cdr1" }
    mutationsFr2Index = header.findIndexOf { it.toLowerCase() == "mutations.nt.fr2" }
    mutationsCdr2Index = header.findIndexOf { it.toLowerCase() == "mutations.nt.cdr2" }
    mutationsFr3Index = header.findIndexOf { it.toLowerCase() == "mutations.nt.fr3" }
    mutationsFr4Index = header.findIndexOf { it.toLowerCase() == "mutations.nt.fr4" }

    if ([freqColIndex, countColIndex, contigNtIndex, cdr3NtIndex,
         mutationsFr1Index, mutationsCdr1Index, mutationsFr2Index, mutationsCdr2Index,
         mutationsFr3Index, mutationsFr4Index].any { it < 0 }) {
        Util.error("One or more the critical columns " +
                "('freq', 'count', 'contignt', 'cdr3nt', 'mutations.nt.fr1/cdr1/fr2/cdr2/fr3/fr4') " +
                "are missing.", 3)
    }
}


def stm = new SequenceTreeMap<NucleotideSequence, List<CorrectorClonotypeEntry>>(NucleotideSequence.ALPHABET)
def searchParams = new TreeSearchParameters(depth, depth, depth, depth)

new File(inputFileName).splitEachLine("\t") { List<String> splitLine ->

    if (!header) {
        header = splitLine
        parseColumns()
    } else {
        def seq = contig ? new NucleotideSequence(splitLine[contigNtIndex]) :
                new NucleotideSequence(splitLine[cdr3NtIndex])
        def entry = new CorrectorClonotypeEntry(splitLine[countColIndex].toInteger(),
                splitLine[freqColIndex].toFloat(),
                splitLine, seq)
        if (!contig) {
            splitLine[[mutationsFr1Index, mutationsCdr1Index,
                       mutationsFr2Index, mutationsCdr2Index,
                       mutationsFr3Index, mutationsFr4Index]].collect { it.split(",") }.flatten().each { String it ->
                if (it.length() > 0) {
                    entry.mutations.add(Mutation.fromString(it))
                }
            }
        }

        stm.createIfAbsent(seq, new Factory<List<CorrectorClonotypeEntry>>() {
            @Override
            List<CorrectorClonotypeEntry> create() {
                new ArrayList<CorrectorClonotypeEntry>()
            }
        }).add(entry)
    }
}

def countNonCdr3Mutations = { CorrectorClonotypeEntry parent, CorrectorClonotypeEntry child ->
    if (!parent.mutations.findAll {
        !child.mutations.contains(it) &&
                child.data[contigNtIndex][it.pos] != "N"
    }.empty)
        return -1 // parsimony principle

    child.mutations.findAll {
        !parent.mutations.contains(it) &&
                parent.data[contigNtIndex][it.pos] != "N"
    }.size()
}

stm.values().each { children ->
    children.each { child ->
        def iter = stm.getNeighborhoodIterator(child.seq, searchParams)
        float bestParentScore = 1.0

        def bestParent = (CorrectorClonotypeEntry) null, parents

        while ((parents = iter.next())) {
            parents.each { parent ->
                int mutationsCount = iter.mutationsCount

                if (!contig) {
                    int additionalMutations = countNonCdr3Mutations(parent, child)
                    if (additionalMutations < 0)
                        return
                    mutationsCount += additionalMutations
                }

                if (parent.count > child.count) {
                    def score = Math.log(child.count / (float) parent.count) / mutationsCount

                    if (score < bestParentScore) {
                        bestParent = parent
                        bestParentScore = score
                    }
                }
            }
        }

        if (bestParentScore <= errorRate) {
            child.parent = bestParent
        }
    }
}

def sortedEntries = stm.values().flatten().sort { CorrectorClonotypeEntry it -> -it.count }

sortedEntries.reverseEach { CorrectorClonotypeEntry it ->
    if (it.parent != null) {
        it.parent.append(it)
    }
}

new File(outputFileName).withPrintWriter { pw ->
    pw.println(header.join("\t"))

    sortedEntries.each { CorrectorClonotypeEntry it ->
        if (it.parent == null) {
            it.data[countColIndex] = it.count.toString()
            it.data[freqColIndex] = it.freq.toString()
            pw.println(it.data.join("\t"))
        }
    }
}


