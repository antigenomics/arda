package com.antigenomics.migmap.post.analysis

import com.antigenomics.migmap.pipeline.Util
import com.antigenomics.migmap.clonotype.Clonotype
import com.antigenomics.migmap.genomic.Segment
import com.antigenomics.migmap.mutation.Mutation

class Analysis {
    final Map<Clonotype, Integer> sample

    Analysis(List<Clonotype> clonotypes) {
        int index = 0
        this.sample = clonotypes.collectEntries { [(it): index++] }
    }

    void generateClonotypeTree(String fileNamePrefix) {
        Util.report("Grouping clonotypes by CDR3")
        def clusters = ClonotypeTree.groupClonotypesByCdr3(sample.keySet())

        // This is for coloring : similar CDR3 gets similar indices
        int index = 0
        def sortByCdr3Indices = new HashSet<String>(sample.keySet()*.cdr3aa)
                .sort()
                .collectEntries { [(it): index++] }

        Util.report("Building edges and writing network")
        new File(fileNamePrefix + ".net.txt").withPrintWriter { pwNet ->
            pwNet.println("from\tto\tinteraction")
            new File(fileNamePrefix + ".edge.txt").withPrintWriter { pwEdge ->
                pwEdge.println("edge\tshm.count\tshm.count.neg\treplacement.ratio\tdifference")
                new File(fileNamePrefix + ".node.txt").withPrintWriter { pwNode ->
                    pwNode.println("node\tcdr3aa\tv\tj\tcdr3.code\tcluster.index\tfreq")
                    clusters.eachWithIndex { cluster, clusterIndex ->
                        def edges = ClonotypeTree.getEdges(cluster)

                        def nodes = new HashSet<Clonotype>()
                        edges.each { edge ->
                            pwNet.println(sample[edge.from] + "\t" + sample[edge.to] + "\tshm")
                            pwEdge.println(sample[edge.from] + " (shm) " + sample[edge.to] + "\t" +
                                    edge.difference.size() + "\t" + (-edge.difference.size()) + "\t" +
                                    calcReplacementFraction(edge.difference) + "\t" +
                                    edge.difference*.toStringAa().join(","))
                            nodes.add(edge.from)
                            nodes.add(edge.to)
                        }

                        nodes.each {
                            pwNode.println(sample[it] + "\t" +
                                    it.cdr3aa + "\t" + it.vSegment.name + "\t" + it.jSegment.name + "\t" +
                                    sortByCdr3Indices[it.cdr3aa] + "\t" + clusterIndex + "\t" + it.freq)
                        }
                    }
                }
            }
        }

        Util.report("Done")
    }


    static double calcReplacementFraction(Collection<Mutation> mutations) {
        mutations.findAll { it.aaFrom != it.aaTo }.size() / (double) mutations.size()
    }

    void generateHypermutationTable(String fileName) {
        def segmentMap = new HashMap<Segment, Counter>()

        Util.report("Summarizing hypermutations")
        sample.keySet().each { Clonotype clonotype ->
            Counter counter
            [clonotype.vSegment, clonotype.dSegment, clonotype.jSegment].findAll { !Segment.isDummy(it) }.each {
                segmentMap.put(it, counter = (segmentMap[it] ?: new Counter()))
                counter.update(clonotype)
            }
        }

        Util.report("Writing statistics")
        new File(fileName).withPrintWriter { pw ->
            pw.println("clonotype.id\tclonotype.v\tclonotype.j\tsegment\tsegment.name\tregion\tmutation.type\t" +
                    "pos.nt\tfrom.nt\tto.nt\tpos.aa\tfrom.aa\tto.aa\t" +
                    "count.clonotypes\tcount.reads\tcount.freq\t" +
                    "total.clonotypes\ttotal.count\ttotal.freq")
            sample.each {
                def clonotype = it.key, ind = it.value
                clonotype.mutations.each { mutation ->
                    def segmentCounter = segmentMap[mutation.parent]
                    pw.println([ind,
                                clonotype.vSegment.name,
                                clonotype.jSegment.name,
                                mutation.parent.type,
                                mutation.parent.name,
                                mutation.subRegion,
                                mutation.type,
                                mutation.pos,
                                mutation.ntFrom,
                                mutation.ntTo,
                                mutation.aaPos,
                                mutation.aaFrom,
                                mutation.aaTo,
                                1,
                                clonotype.count,
                                clonotype.freq,
                                segmentCounter.numberOfClonotypes,
                                segmentCounter.count,
                                segmentCounter.frequency
                    ].join("\t"))
                }
            }
        }
        Util.report("Done")
    }
}
