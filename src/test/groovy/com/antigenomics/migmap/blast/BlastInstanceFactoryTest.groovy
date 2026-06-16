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
import com.antigenomics.migmap.pipeline.Util
import org.junit.Test

class BlastInstanceFactoryTest {
    @Test
    void executeTest() {
        // Should fail with error = 2 if software is not present
        println ExecutionUtil.makeDb.execute().text
        println ExecutionUtil.igBlast.execute().text
    }

    @Test
    void createTest() {
        def factory = new BlastInstanceFactory("data/", "human", ["IGH"])

        def instances = (1..Util.N_THREADS).collect { factory.create() }
        instances.each { BlastInstance it -> it.put(null); it.close() }
    }
}
