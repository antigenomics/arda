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

package com.antigenomics.migmap.mutation

import groovy.transform.CompileStatic

@CompileStatic
enum SubRegion {
    FR1(0), CDR1(1), FR2(2), CDR2(3), FR3(4), CDR3(5), FR4(6)

    final int order

    SubRegion(int order) {
        this.order = order
    }
    
    static final SubRegion[] REGION_LIST = [FR1, CDR1, FR2, CDR2, FR3, CDR3, FR4] as SubRegion[]
}