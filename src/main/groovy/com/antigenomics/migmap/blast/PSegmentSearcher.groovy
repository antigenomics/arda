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
import com.antigenomics.migmap.mapping.Cdr3Markup
import com.antigenomics.migmap.mapping.Truncations
import groovy.transform.CompileStatic

@CompileStatic
public class PSegmentSearcher {
    final static int P_SEGM_MIN_LEN = 3

    static PSegments search(Cdr3Markup cdr3Markup, Truncations truncations, String cdr3seq) {
        boolean dFound = cdr3Markup.dStart >= 0
        int vp3 = -1, dp5 = -1, dp3 = -1, jp5 = -1

        if (truncations.vDel == 0) {
            vp3 = searchRight(cdr3Markup.vEnd + 1, dFound ? cdr3Markup.dStart : cdr3Markup.jStart, 0, cdr3seq)
        }

        if (truncations.dDel5 == 0) {
            dp5 = searchLeft(cdr3Markup.dStart - 1, cdr3Markup.vEnd, cdr3Markup.dEnd, cdr3seq)
        }

        if (truncations.dDel3 == 0) {
            dp3 = searchRight(cdr3Markup.dEnd + 1, cdr3Markup.jStart, cdr3Markup.dStart, cdr3seq)
        }

        if (truncations.jDel == 0) {
            jp5 = searchLeft(cdr3Markup.jStart - 1, dFound ? cdr3Markup.vEnd : cdr3Markup.dEnd, cdr3seq.length(), cdr3seq)
        }

        new PSegments(
                vp3 - cdr3Markup.vEnd >= P_SEGM_MIN_LEN ? vp3 : -1,
                cdr3Markup.dStart - dp5 >= P_SEGM_MIN_LEN ? dp5 : -1,
                dp3 - cdr3Markup.dEnd >= P_SEGM_MIN_LEN ? dp3 : -1,
                cdr3Markup.jStart - jp5 >= P_SEGM_MIN_LEN ? jp5 : -1)
    }

    static int searchRight(int from, int to, int from2, String cdr3seq) {
        int lastMatch = -1, lastMm = -2
        for (int i = from; i < to; i++) {
            int j = 2 * from - i - 1
            if (j < from2)
                break
            if (cdr3seq.charAt(i) == Util.compl(cdr3seq.charAt(j))) {
                lastMatch = i
            } else {
                if (i - lastMm <= 2) { // no more than 2 consequent mismatches or 2 mismatches separated by 1 base
                    lastMatch = lastMm - 1
                    break
                }
                lastMm = i
            }
        }
        lastMatch
    }

    static int searchLeft(int from, int to, int from2, String cdr3seq) {
        int lastMatch = -1, lastMm = cdr3seq.length() + 1
        for (int i = from; i > to; i--) {
            int j = 2 * from - i + 1
            if (j >= from2)
                break
            if (cdr3seq.charAt(i) == Util.compl(cdr3seq.charAt(j))) {
                lastMatch = i
            } else {
                if (lastMm - i <= 2) { // no more than 2 consequent mismatches or 2 mismatches separated by 1 base
                    lastMatch = lastMm + 1
                    break
                }
                lastMm = i
            }
        }
        lastMatch
    }

}
