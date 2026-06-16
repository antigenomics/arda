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

import groovy.transform.CompileStatic

@CompileStatic
class Truncations implements Serializable {
    final int vDel, dDel5, dDel3, jDel

    Truncations(int vDel, int dDel5, int dDel3, int jDel) {
        this.vDel = vDel
        this.dDel5 = dDel5
        this.dDel3 = dDel3
        this.jDel = jDel
    }

    static final String OUTPUT_HEADER = "v.del\td.del.5\td.del.3\tj.del"

    @Override
    public String toString() {
        [vDel, dDel5, dDel3, jDel].join("\t")
    }
}
