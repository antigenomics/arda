package com.antigenomics.migmap.mutation

import com.antigenomics.migmap.pipeline.Util
import com.antigenomics.migmap.clonotype.Clonotype
import com.antigenomics.migmap.mapping.ReadMapping

class MutationConverter {
    static void annotateMutationAa(Clonotype clonotype) {
        if (!clonotype.mutations.empty && clonotype.mutations.first().aaFrom != null)
            return // already annotated
        annotateMutationAa(clonotype.representativeMapping)
    }

    static void annotateMutationAa(ReadMapping readMapping) {
        if (readMapping.mapping) {
            annotateMutationAa(
                    readMapping.mapping.mutations,
                    readMapping.mapping.rc ? readMapping.read.rc.seq : readMapping.read.seq,
                    readMapping.mapping.vStartInRef, readMapping.mapping.vStartInQuery)
        }
    }

    static void annotateMutationAa(List<Mutation> mutations, String rawQuery, int vStartInRef, int vStartInQuery) {
        String ref = Util.translateLinear('N' * vStartInRef + mutateBack(rawQuery, mutations).substring(vStartInQuery)),
               query = Util.translateLinear('N' * vStartInRef + rawQuery.substring(vStartInQuery))

        mutations.each {
            annotateMutationAa(it, vStartInQuery, vStartInRef, query, ref)
        }
    }

    static void annotateMutationAa(Mutation mutation, int vStartInQuery, int vStartInRef, String query, String ref) {
        int startInQuery = (int) ((mutation.startInRead - vStartInQuery + vStartInRef) / 3),
            endInQuery = (int) ((mutation.endInRead - vStartInQuery + vStartInRef - 1) / 3),
            startInRef = (int) (mutation.start / 3),
            endInRef = (int) ((mutation.end - 1) / 3)

        mutation.aaPos = (int) (mutation.pos / 3)
        mutation.aaFrom = mutation.type == MutationType.Insertion ? "" :
                (endInRef >= ref.length() || startInRef < 0 ? "X" :
                ref[startInRef..endInRef])
        mutation.aaTo = mutation.type == MutationType.Deletion ? "" :
                (endInQuery >= ref.length() || startInQuery < 0 ? "X" :
                query[startInQuery..endInQuery])
    }

    static String mutateBack(String readSeq, List<Mutation> mutations) {
        List<String> mutatedSeq = readSeq.toCharArray().collect { it.toString() }

        mutations.each { mut ->
            switch (mut.type) {
                case MutationType.Substitution:
                    mutatedSeq[mut.startInRead] = mut.ntFrom
                    break
                case MutationType.Deletion:
                    mutatedSeq[mut.endInRead] = mut.ntFrom + mutatedSeq[mut.endInRead]
                    break
                case MutationType.Insertion:
                    (mut.startInRead..<mut.endInRead).each { int it -> mutatedSeq[it] = "" }
                    break
            }
        }

        mutatedSeq.join("")
    }
}
