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

package com.antigenomics.migmap.post.analysis

import com.antigenomics.migmap.pipeline.Util
import com.antigenomics.migmap.clonotype.Clonotype
import com.antigenomics.migmap.mutation.Mutation
import com.milaboratory.core.alignment.Alignment
import com.milaboratory.core.sequence.NucleotideSequence
import com.milaboratory.core.tree.SequenceTreeMap
import com.milaboratory.core.tree.TreeSearchParameters
import com.milaboratory.util.Factory
import groovyx.gpars.GParsPool

import java.util.concurrent.ConcurrentHashMap

class ClonotypeTree {
    static final TreeSearchParameters germOnlyDiffSearchParams = new TreeSearchParameters(5, 0, 0, 5)

    // Grouping by CDR3
    static List<List<Clonotype>> groupClonotypesByCdr3(Iterable<Clonotype> sample) {
        def cTree = new SequenceTreeMap<NucleotideSequence,
                CdrNode>(NucleotideSequence.ALPHABET)

        Util.report("- Pre-group by CDR3")
        sample.each { Clonotype clonotype ->
            def cdr3Nt = new NucleotideSequence(clonotype.cdr3nt)
            def cdrNode = cTree.createIfAbsent(cdr3Nt, new Factory<CdrNode>() {
                @Override
                CdrNode create() {
                    new CdrNode(clonotype)
                }
            })
            cdrNode.clonotypes.add(clonotype)
        }

        Util.report("- Finding CDR3 differences explained by SHMs in germline region")

        def cdrNodeMappingsMap = new ConcurrentHashMap<CdrNode, Set<CdrNode>>()

        GParsPool.withPool(Util.N_THREADS) {
            cTree.values().eachParallel { CdrNode cdrNode ->
                def niter = cTree.getNeighborhoodIterator(cdrNode.cdr3,
                        germOnlyDiffSearchParams)

                def prevCdr3Hash = new HashSet<NucleotideSequence>()
                prevCdr3Hash.add(cdrNode.cdr3)
                def mappings = new HashSet<CdrNode>()

                def nextCdrNode
                while ((nextCdrNode = niter.next()) != null) {
                    if (!prevCdr3Hash.contains(nextCdrNode.cdr3)) {
                        prevCdr3Hash.add(nextCdrNode.cdr3)

                        def alignment = niter.currentAlignment

                        if (cdr3MatchWithoutGermlineMutations(cdrNode.representative,
                                nextCdrNode.representative,
                                nextCdrNode.cdr3,
                                alignment)) {
                            mappings.add(nextCdrNode)
                        }
                    }
                }
                cdrNodeMappingsMap.put(cdrNode, mappings)
            }
        }

        Util.report("- Aggregating")
        def usedNodes = new HashSet<CdrNode>()

        def addMappings = { Set<CdrNode> nodes ->
            def origSize = nodes.size()

            nodes.collect().each {
                nodes.addAll(cdrNodeMappingsMap[it])
            }

            origSize < nodes.size()
        }

        def groupedClonotypes = new ArrayList<List<Clonotype>>()

        cdrNodeMappingsMap.keySet().each {
            if (!usedNodes.contains(it)) {
                def mappings = new HashSet<CdrNode>()
                mappings.add(it)

                while (addMappings(mappings));

                usedNodes.addAll(mappings)

                def clonotypeGroup = new ArrayList<Clonotype>()

                mappings.each {
                    clonotypeGroup.addAll(it.clonotypes)
                }

                groupedClonotypes.add(clonotypeGroup)
            }
        }

        groupedClonotypes
    }

    static boolean inGermline(int position, Clonotype clonotype) {
        clonotype.cdr3Markup.vEnd >= position || clonotype.cdr3Markup.jStart <= position ||
                (clonotype.cdr3Markup.dEnd >= position && clonotype.cdr3Markup.dStart <= position)
    }

    static boolean cdr3MatchWithoutGermlineMutations(Clonotype clonotype1, Clonotype clonotype2,
                                                     NucleotideSequence cdr3nt2,
                                                     Alignment<NucleotideSequence> aln) {
        def mut12 = aln.absoluteMutations, mut21 = aln.invert(cdr3nt2).absoluteMutations

        // in germline in at least one of sequences

        for (int i = 0; i < mut12.size(); i++) {
            if (!inGermline(mut12.getPositionByIndex(i), clonotype1) &&
                    !inGermline(mut21.getPositionByIndex(i), clonotype2))
                return false
        }
        true
    }

    // Graph

    static List<Edge> getEdges(List<Clonotype> clonotypes) {
        // Try connect clonotypes
        def edgeList = new ArrayList<Edge>()
        def parentNodeMap = new HashMap<Clonotype, Set<Clonotype>>()
        for (int i = 0; i < clonotypes.size(); i++) {
            for (int j = i + 1; j < clonotypes.size(); j++) {
                def edge = connect(clonotypes[i], clonotypes[j])
                if (edge) {
                    edgeList.add(edge)
                    def parentList
                    parentNodeMap.put(edge.to,
                            parentList = (parentNodeMap[edge.to] ?: new HashSet<>()))
                    parentList.add(edge.from)
                }
            }
        }

        // Remove redundancy
        edgeList.findAll { edge ->
            //
            //      2 <-- 3
            //      ^     ^
            //       \   /
            //         1
            //
            // remove 1->2
            // any parent(2) that is not 1 has 1 as parent
            !parentNodeMap[edge.to].any { it != edge.from && hasParent(parentNodeMap[it], edge.from) }
        }
    }

    private static boolean hasParent(Set<Clonotype> clonotypes, Clonotype clonotype) {
        clonotypes && clonotypes.contains(clonotype)
    }

    static Set<String> getIntersection(List<Mutation> mutations1, List<Mutation> mutations2) {
        def set = new HashSet<String>(mutations1*.toString())
        new HashSet<>(mutations2*.toString().findAll { set.contains(it.toString()) })
    }

    static Edge connect(Clonotype clonotype1, Clonotype clonotype2) {
        def mutations1 = clonotype1.mutations,
            mutations2 = clonotype2.mutations

        def intersection = mutations1.size() > mutations2.size() ?
                getIntersection(mutations1, mutations2) : getIntersection(mutations2, mutations1)

        intersection.size() == mutations1.size() ?
                new Edge(clonotype1, clonotype2, mutations2.findAll {
                    !intersection.contains(it)
                }) :
                (intersection.size() == mutations2.size() ?
                        new Edge(clonotype2, clonotype1, mutations1.findAll {
                            !intersection.contains(it)
                        }) :
                        null)
    }
}
