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
import com.antigenomics.migmap.post.ContigMergerClonotypeEntry

def cli = new CliBuilder(usage: "MergeContigs [options] clonotype_table.txt output.txt")

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

// Parsing required columns, here we operate on full contigs

def countColIndex = -1, freqColIndex = -1, contigNtIndex = -1

def header
def parseColumns = {
    freqColIndex = header.findIndexOf { it.toLowerCase() == "freq" }
    countColIndex = header.findIndexOf { it.toLowerCase() == "count" }
    contigNtIndex = header.findIndexOf { it.toLowerCase() == "contignt" }

    if ([freqColIndex, countColIndex, contigNtIndex].any { it < 0 }) {
        Util.error("One or more the critical columns " +
                "('freq', 'count', 'contignt') " +
                "are missing.", 3)
    }
}


def entries = new ArrayList<ContigMergerClonotypeEntry>()

new File(inputFileName).splitEachLine("\t") { List<String> splitLine ->

    if (!header) {
        header = splitLine
        parseColumns()
    } else {
        def entry = new ContigMergerClonotypeEntry(splitLine[countColIndex].toInteger(),
                splitLine[freqColIndex].toFloat(), splitLine[contigNtIndex].replaceAll("N", ""),
                splitLine)

        entries.add(entry)
    }
}

// Assume already sorted

for (int i = 0; i < entries.size(); i++) {
    def parent = entries[i]
    if (!parent.appended) {
        for (int j = i + 1; j < entries.size(); j++) {
            def child = entries[j]

            if (!child.appended &&
                    child.contingWithoutNs.contains(parent.contingWithoutNs) ||
                    parent.contingWithoutNs.contains(child.contingWithoutNs)) {
                child.appended = true
                parent.append(child)
            }
        }
    }
}

new File(outputFileName).withPrintWriter { pw ->
    pw.println(header.join("\t"))
    entries.sort { -it.count }.each {
        if (!it.appended) {
            it.data[countColIndex] = it.count.toString()
            it.data[freqColIndex] = it.freq.toString()
            pw.println(it.data.join("\t"))
        }
    }
}