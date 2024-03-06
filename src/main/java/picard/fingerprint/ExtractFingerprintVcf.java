/*
 * The MIT License
 *
 * Copyright (c) 2020 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package picard.fingerprint;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.Hidden;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Map;
 
 /**
  * Program to create a fingerprint for the <b>contaminating</b> sample when the level of contamination is both known and
  * uniform in the genome.
  *
  * @author Yossi Farjoun
  */
 @CommandLineProgramProperties(
         summary = "Computes/Extracts the fingerprint genotype likelihoods from the supplied file." +
                 "It is given as a list of PLs at the fingerprinting sites.",
         oneLineSummary = "Computes a fingerprint from the input file.",
         programGroup = DiagnosticsAndQCProgramGroup.class)
 
 public class ExtractFingerprintVcf extends CommandLineProgram {
 
     @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input VCF file.")
     public String INPUT;

     @Argument(doc = "Input VCF index file.")
     public String INPUT_INDEX_PATH = null;
 
     @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output directory for fingerprint files (VCF).")
     public File OUTPUT;
 
     @Argument(shortName = "H", doc = "A file of haplotype information. The file lists a set of SNPs, optionally arranged in high-LD blocks, to be used for fingerprinting. See " +
             "https://software.broadinstitute.org/gatk/documentation/article?id=9526 for details.")
     public File HAPLOTYPE_MAP;
 
     @Argument(doc = "The maximum number of reads to use as evidence for any given locus. This is provided as a way to limit the " +
             "effect that any given locus may have.")
     public int LOCUS_MAX_READS = 50;
 
     @Hidden
     @Argument(doc = "When true code will check for readability on input files (this can be slow on cloud access)")
     public boolean TEST_INPUT_READABILITY = true;
 
     private static final Log log = Log.getInstance(ExtractFingerprint.class);

    // Do we actually need a reference?
    @Override
    protected boolean requiresReference() {
        return true;
    }

    @Override
    protected int doWork() {

    final Path inputPath;
    try {
        inputPath = IOUtil.getPath(INPUT);
    } catch (IOException e) {
        throw new RuntimeException(e);
    }
    if (TEST_INPUT_READABILITY) {
        IOUtil.assertFileIsReadable(inputPath);
    }
    final Path indexPath;
    try {
        indexPath = IOUtil.getPath(INPUT_INDEX_PATH);
    } catch (IOException e) {
        throw new RuntimeException(e);
    }
    if (TEST_INPUT_READABILITY) {
        IOUtil.assertFileIsReadable(indexPath);
    }

    IOUtil.assertFileIsReadable(referenceSequence.getReferenceFile());
    IOUtil.assertDirectoryIsWritable(OUTPUT);
    IOUtil.assertFileIsReadable(HAPLOTYPE_MAP);

    final FingerprintChecker checker = new FingerprintChecker(HAPLOTYPE_MAP);

    checker.setLocusMaxReads(LOCUS_MAX_READS);
    checker.setValidationStringency(VALIDATION_STRINGENCY);

    final Map<FingerprintIdDetails, Fingerprint> oneFileFingerprints;

    oneFileFingerprints = checker.fingerprintVcf(inputPath, indexPath, true);

    for (Map.Entry<FingerprintIdDetails, Fingerprint> entry : oneFileFingerprints.entrySet()) {
        String sample = entry.getKey().getSample();
        System.out.println("Sample: " + sample);

        Fingerprint sampleFingerprint = entry.getValue();
        try {
            File outputFile = new File(OUTPUT, sample + ".fingerprint.vcf");
            FingerprintUtils.writeFingerPrint(sampleFingerprint, outputFile, referenceSequence.getReferenceFile(), sample, "PLs derived from " + INPUT + " using an assumed contamination of 0");
        } catch (Exception e) {
            log.error(e);
        }
    }

    System.out.println(inputPath);
    return 0;
    }
 }
 