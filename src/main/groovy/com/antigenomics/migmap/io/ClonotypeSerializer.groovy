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

package com.antigenomics.migmap.io

import com.antigenomics.migmap.clonotype.Clonotype

class ClonotypeSerializer {
    static void save(List<Clonotype> clonotypes, OutputStream outputStream) {
        def out = new ObjectOutputStream(outputStream)
        out.writeObject(clonotypes)
        out.close()
    }

    static List<Clonotype> load(InputStream inputStream) {
        new ObjectInputStream(inputStream).readObject() as List<Clonotype>
    }
}
