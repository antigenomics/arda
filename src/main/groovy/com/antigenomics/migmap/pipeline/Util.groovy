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

package com.antigenomics.migmap.pipeline

import com.antigenomics.migmap.genomic.SegmentDatabase
import groovy.transform.CompileStatic

import java.util.concurrent.TimeUnit
import java.util.regex.Matcher
import java.util.zip.GZIPInputStream
import java.util.zip.GZIPOutputStream


class Util {
    static final int N_THREADS = Runtime.runtime.availableProcessors()
    static final String MY_NA = ".", BLAST_NA = "N/A"
    static final byte MAX_QUAL = 40, MIN_QUAL = 2
    static final char GAP = '-'
    static final String MAX_QUAL_SYMBOL = (char) (MAX_QUAL + 33)
    static final int QUAL_OFFSET = 33
    static int VERBOSITY_LEVEL = 2

    private static Date start = null

    @CompileStatic
    private static String timePassed(long millis) {
        String.format("%02dm%02ds",
                TimeUnit.MILLISECONDS.toMinutes(millis),
                TimeUnit.MILLISECONDS.toSeconds(millis) -
                        TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(millis))
        )
    }

    @CompileStatic
    static void report(String message, int verbosity = 1) {
        if (verbosity <= VERBOSITY_LEVEL) {
            if (start == null)
                start = new Date()

            System.err.println("[${new Date()} MIGMAP +" + timePassed(new Date().time - start.time) + "] $message")
        }
    }

    @CompileStatic
    static double toPercent(double ratio) {
        ((int) (ratio * 10000) / 100)
    }

    @CompileStatic
    static void error(String message, int code) {
        // code
        // 1 - runtime error in wrapper
        // 2 - blast error
        // 3 - error in main script
        System.err.println("[${new Date()} MIGMAP ERROR] $message")
        SegmentDatabase.clearTemporaryFiles()
        System.exit(code)
    }

    @CompileStatic
    static byte minQual(byte[] qual) {
        byte min = MAX_QUAL
        for (byte q : qual) {
            min = (byte) Math.min(q, min)
        }
        min
    }

    @CompileStatic
    static String qualToString(byte[] qual) {
        qual.collect { (char) ((int) it + 33) }.join("")
    }

    @CompileStatic
    static InputStream getStream(String fname, boolean resource) {
        resource ? resourceGetHelper(fname) : new FileInputStream(fname)
    }

    @CompileStatic
    private static InputStream resourceGetHelper(String fname) {
        def stream = Util.class.classLoader.getResourceAsStream(fname)
        if (!stream) {
            fname = ["./build/resources/main/" + fname,
                     "./build/resources/test/" + fname].find { String it -> new File(it).exists() }
            return new FileInputStream(fname)
        } else {
            return stream
        }
    }

    @CompileStatic
    static BufferedReader getReader(String fname, boolean resource) {
        def inputStream = getStream(fname, resource)
        new BufferedReader(new InputStreamReader(fname.endsWith(".gz") ? new GZIPInputStream(inputStream) : inputStream))
    }

    @CompileStatic
    static BufferedWriter getWriter(String outfile) {
        def compressed = outfile.endsWith(".gz")
        new BufferedWriter(new OutputStreamWriter(compressed ?
                new GZIPOutputStream(new FileOutputStream(outfile)) : new FileOutputStream(outfile)))
    }


    static List<String> groomMatch(Matcher matcher) {
        matcher.size() > 0 ? matcher[0][1..-1] : null//[]
    }


    static final char ntA = 'A', ntT = 'T', ntG = 'G', ntC = 'C', ntN = 'N',
                      nta = 'a', ntt = 't', ntg = 'g', ntc = 'c', ntn = 'n'

    @CompileStatic
    static char compl(char nt) {
        switch (nt) {
            case ntA:
                return ntT
            case ntT:
                return ntA
            case ntG:
                return ntC
            case ntC:
                return ntG
            case ntN:
                return ntN
            case nta:
                return ntt
            case ntt:
                return nta
            case ntg:
                return ntc
            case ntc:
                return ntg
            case ntn:
                return ntn
            default:
                return ntN
        }
    }

    @CompileStatic
    static String revCompl(String seq) {
        char[] chars = seq.reverse().toCharArray()
        for (int i = 0; i < chars.length; i++) {
            chars[i] = compl(chars[i])
        }
        new String(chars)
    }

    @CompileStatic
    static String codon2aa(String codon) {
        String codonUpper = codon.toUpperCase()
        switch (codonUpper) {
            case 'TTT': return 'F'
            case 'TTC': return 'F'
            case 'TTA': return 'L'
            case 'TTG': return 'L'
            case 'TCT': return 'S'
            case 'TCC': return 'S'
            case 'TCA': return 'S'
            case 'TCG': return 'S'
            case 'TAT': return 'Y'
            case 'TAC': return 'Y'
            case 'TAA': return '*'
            case 'TAG': return '*'
            case 'TGT': return 'C'
            case 'TGC': return 'C'
            case 'TGA': return '*'
            case 'TGG': return 'W'
            case 'CTT': return 'L'
            case 'CTC': return 'L'
            case 'CTA': return 'L'
            case 'CTG': return 'L'
            case 'CCT': return 'P'
            case 'CCC': return 'P'
            case 'CCA': return 'P'
            case 'CCG': return 'P'
            case 'CAT': return 'H'
            case 'CAC': return 'H'
            case 'CAA': return 'Q'
            case 'CAG': return 'Q'
            case 'CGT': return 'R'
            case 'CGC': return 'R'
            case 'CGA': return 'R'
            case 'CGG': return 'R'
            case 'ATT': return 'I'
            case 'ATC': return 'I'
            case 'ATA': return 'I'
            case 'ATG': return 'M'
            case 'ACT': return 'T'
            case 'ACC': return 'T'
            case 'ACA': return 'T'
            case 'ACG': return 'T'
            case 'AAT': return 'N'
            case 'AAC': return 'N'
            case 'AAA': return 'K'
            case 'AAG': return 'K'
            case 'AGT': return 'S'
            case 'AGC': return 'S'
            case 'AGA': return 'R'
            case 'AGG': return 'R'
            case 'GTT': return 'V'
            case 'GTC': return 'V'
            case 'GTA': return 'V'
            case 'GTG': return 'V'
            case 'GCT': return 'A'
            case 'GCC': return 'A'
            case 'GCA': return 'A'
            case 'GCG': return 'A'
            case 'GAT': return 'D'
            case 'GAC': return 'D'
            case 'GAA': return 'E'
            case 'GAG': return 'E'
            case 'GGT': return 'G'
            case 'GGC': return 'G'
            case 'GGA': return 'G'
            case 'GGG': return 'G'
            default:
                if (codonUpper.contains("N") && codonUpper.length() == 3)
                    return "X" // undefined
                else
                    return '?' // incomplete/missing
        }
    }

    @CompileStatic
    static String translateCdr(String seq) {
        def aaSeq = ""
        def oof = seq.size() % 3
        if (oof > 0) {
            def mid = (int) (seq.size() / 2)
            seq = seq.substring(0, mid) + ("?" * (3 - oof)) + seq.substring(mid, seq.length())
        }

        def leftEnd = -1, rightEnd = -1
        for (int i = 0; i <= seq.size() - 3; i += 3) {
            def codon = seq.substring(i, i + 3)
            if (codon.contains("?")) {
                leftEnd = i
                break
            }
            aaSeq += codon2aa(codon)
        }

        if (oof == 0)
            return aaSeq

        def aaRight = ""
        for (int i = seq.size(); i >= 3; i -= 3) {
            def codon = seq.substring(i - 3, i)
            if (codon.contains("?")) {
                rightEnd = i
                break
            }
            aaRight += codon2aa(codon)
        }

        aaSeq + seq.substring(leftEnd, rightEnd).toLowerCase() + aaRight.reverse()
    }

    @CompileStatic
    static String removeGaps(String seq) {
        seq.replaceAll("-", "")
    }

    @CompileStatic
    static String translateLinear(String seq) {
        def aaSeq = ""

        for (int i = 0; i <= seq.size() - 3; i += 3) {
            def codon = seq.substring(i, i + 3)
            aaSeq += codon2aa(codon)
        }

        aaSeq
    }

    @CompileStatic
    static boolean isCanonical(String cdr3) {
        cdr3 =~ /^TG[TC].+(?:TGG|TT[TC])$/
    }

}
