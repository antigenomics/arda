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

import com.antigenomics.migmap.blast.BlastInstance
import com.antigenomics.migmap.blast.BlastInstanceFactory
import com.antigenomics.migmap.io.InputPort
import com.antigenomics.migmap.io.OutputPort
import com.antigenomics.migmap.io.Read
import com.antigenomics.migmap.mapping.ReadMapping
import com.antigenomics.migmap.mapping.ReadMappingFilter
import groovy.transform.CompileStatic

import java.util.concurrent.atomic.AtomicLong

@CompileStatic
class Pipeline {
    private final AtomicLong inputCounter = new AtomicLong()
    final long limit
    final int nThreads
    final OutputPort<Read> input
    final BlastInstanceFactory blastInstanceFactory
    final InputPort<ReadMapping> output
    final ReadMappingFilter readMappingFilter

    Pipeline(OutputPort<Read> input, BlastInstanceFactory blastInstanceFactory,
             InputPort<ReadMapping> output, ReadMappingFilter readMappingFilter,
             long limit = -1, int nThreads = Util.N_THREADS) {
        this.input = input
        this.blastInstanceFactory = blastInstanceFactory
        this.output = output
        this.readMappingFilter = readMappingFilter
        this.limit = limit < 0 ? Long.MAX_VALUE : limit
        this.nThreads = nThreads
    }

    void run() {
        def threads = new Thread[2 * nThreads]
        def instances = new BlastInstance[nThreads]

        Util.report("Started analysis", 2)

        // read >> blast instance threads
        (0..<nThreads).each { int it ->
            BlastInstance instance = blastInstanceFactory.create()
            instances[it] = instance
            threads[it] = new Thread(new Runnable() {
                @Override
                void run() {
                    def read
                    while (((read = input.take()) != null) &&
                            (inputCounter.incrementAndGet() <= limit)
                    ) {
                        instance.put(read)
                    }
                    instance.put(null) // finished
                }
            })
        }

        // NOTE: here blast instance acts like a buffer

        // blast instance >> output threads
        (0..<nThreads).each {
            BlastInstance instance = instances[it]
            threads[nThreads + it] = new Thread(new Runnable() {
                @Override
                void run() {
                    ReadMapping result
                    while ((result = instance.take()) != null) {
                        if (readMappingFilter.pass(result)) {
                            output.put(result)
                        }
                    }
                }
            })
        }

        def reporter = new Thread(new Runnable() {
            @Override
            void run() {
                try {
                    while (!Thread.currentThread().isInterrupted()) {
                        Util.report("Loaded $inputCount reads. " +
                                (readMappingFilter.total > 0L ? readMappingFilter.toProgressString() : ""), 2)
                        Thread.sleep(10000)
                    }
                } catch (InterruptedException e) {

                }
            }
        })

        threads.each { Thread it -> it.start() }

        reporter.daemon = true
        reporter.start()

        threads.each { Thread it -> it.join() }

        reporter.interrupt()

        Util.report("Finished analysis. ${readMappingFilter.toProgressString()}", 2)

        // Close all ports

        readMappingFilter.unmappedInputPort.close()
        output.close()
    }

    final Thread.UncaughtExceptionHandler h = new Thread.UncaughtExceptionHandler() {
        void uncaughtException(Thread t, Throwable e) {
            throw new RuntimeException("Error in pipeline", e)
        }
    }

    long getInputCount() {
        inputCounter.get()
    }
}
