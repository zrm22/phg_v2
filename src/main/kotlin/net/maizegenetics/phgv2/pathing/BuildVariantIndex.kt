package net.maizegenetics.phgv2.pathing

import biokotlin.util.GetVCFVariants
import biokotlin.util.SimpleVariant
import biokotlin.util.bufferedReader
import biokotlin.util.getAllVCFFiles
import biokotlin.util.validateVCFs
import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.github.ajalt.clikt.parameters.types.int
import com.google.common.collect.Range
import com.google.common.collect.RangeMap
import com.google.common.collect.TreeRangeMap
import htsjdk.variant.variantcontext.Allele
import htsjdk.variant.variantcontext.Genotype
import htsjdk.variant.variantcontext.GenotypeBuilder
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.VariantContextBuilder
import htsjdk.variant.variantcontext.writer.Options
import htsjdk.variant.variantcontext.writer.VariantContextWriter
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder
import htsjdk.variant.vcf.VCFFileReader
import htsjdk.variant.vcf.VCFHeaderLine
import kotlinx.coroutines.CoroutineScope
import kotlinx.coroutines.Deferred
import kotlinx.coroutines.Dispatchers
import kotlinx.coroutines.Job
import kotlinx.coroutines.async
import kotlinx.coroutines.channels.Channel
import kotlinx.coroutines.channels.SendChannel
import kotlinx.coroutines.channels.consumeEach
import kotlinx.coroutines.joinAll
import kotlinx.coroutines.launch
import kotlinx.coroutines.runBlocking
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.utils.AltHeaderMetaData
import net.maizegenetics.phgv2.utils.Position
import net.maizegenetics.phgv2.utils.VariantContextUtils
import net.maizegenetics.phgv2.utils.createGenericHeader
import org.apache.logging.log4j.LogManager
import java.io.File

class BuildVariantIndex : CliktCommand(help="Build PHG Index from the Variants") {
    private val myLogger = LogManager.getLogger(BuildVariantIndex::class.java)

    //Need either a directory of the gvcf files.
    //The gvcfs are stored in dbPath
    // val uri = dbPath + "/gvcf_dataset"
    val dbPath by option(help = "Tile DB URI")
        .required() //Needs to be required now due to the agc archive

    val indexFile by option(help = "The full path of the variant index file.")
        .required()

    val gvcfPath by option(help = "The full path of the variant gvcf file.")
        .required()

    val numThreads by option(help = "Number of threads to use for processing gvcf files.")
        .int()
        .default(5)

    val outputFile by option(help = "Output file")
        .required()


    override fun run() {
        myLogger.info("Building Variant Index")
        buildVariantIndex(dbPath, indexFile, gvcfPath, numThreads, outputFile)
    }

    fun buildVariantIndex(dbPath: String, indexFile: String, gvcfPath: String, numThreads: Int = 5, outputFile:String) {
        //Implement function to build the variant index from the tileDB
        //Build the phg mapping for SampleName + Position Map -> HaplotypeID
        //TODO uncomment this after we test speed.
//        val samplePositionToHapIdMap = buildSamplePositionToHapIdMap(dbPath)

        //Now we need to walk through the gvcf files and build a SNP list
//        val variants = buildVariantSet(gvcfPath)
        val variants = processAllGVCFsMultithread(
            File(gvcfPath).walkTopDown().filter{
                it.isFile && (it.name.endsWith(".gvcf") || it.name.endsWith(".gvcf.gz") || it.name.endsWith("g.vcf") || it.name.endsWith(".g.vcf.gz"))
            }.map { it.absolutePath }.toList(),
            numThreads
        )

        println("Total Variants Found: ${variants.size}")
        //Then we can do a 2nd pass through the gvcf files and build the index using those
        buildMergedVariantContexts(variants, gvcfPath, outputFile)
    }

    /**
     * Function to build a lookup map that will allow us to quickly find the HaplotypeID for a given SampleName and Position
     */
    fun buildSamplePositionToHapIdMap(dbPath: String): Map<String, RangeMap<Position, String>> {
        //Build a graph from the tileDB at dbPath
        val graph = HaplotypeGraph("$dbPath/hvcf_files")

        val sampleGametes = graph.sampleGametesInGraph()

        val refRanges = graph.ranges()
        val samplePositionToHapIdMap = mutableMapOf<String, RangeMap<Position, String>>()

        for(sample in sampleGametes) {
            samplePositionToHapIdMap[sample.name] = TreeRangeMap.create<Position, String>()
        }

        for(refRange in refRanges) {
            for(sample in sampleGametes) {
                val hapId = graph.sampleToHapId(refRange, sample)
                if(hapId != null) {
                    samplePositionToHapIdMap[sample.name]!!.put(refRange.toClosedGuavaRange(), hapId)
                }
            }
        }

        return samplePositionToHapIdMap
    }

    fun buildVariantSet(gvcfPath: String): Set<Position> {
        //walk through all files in the gvcfPath and build a set of all variants
        val variantSet = File(gvcfPath).walkTopDown().filter{
            it.isFile && (it.name.endsWith(".gvcf") || it.name.endsWith(".gvcf.gz") || it.name.endsWith("g.vcf") || it.name.endsWith(".g.vcf.gz"))
        }.flatMap { file ->
            myLogger.info("Processing file: ${file.absolutePath}")
            val startTime = System.nanoTime()
//            val variants = processSingleGVCFForVariants(file.absolutePath)
            val variants = processSingleGVCFForVariantsSimple(file.absolutePath)
            val endTime = System.nanoTime()
            myLogger.info("Finished processing file: ${file.absolutePath} in ${(endTime - startTime) / 1E9} seconds")
            variants
        }.toSet()


        return variantSet
    }

    //This one tests out VariantContext reading.  This will be slower likely
    fun processSingleGVCFForVariants(gvcfFile: String): Set<Position> {
        //Implement function to read a gvcf file and extract the variant positions
        val variantPositions = mutableSetOf<Position>()
        //Read through the gvcf file and extract the variant positions
        VCFFileReader(File(gvcfFile), false).use { vcfReader ->
            val iterator = vcfReader.iterator()
            while (iterator.hasNext()) {
                val variant = iterator.next()
                //Skipping indels here
                if(variant.isVariant && !variant.isIndel) {
                    val pos = Position(variant.contig, variant.start)
                    variantPositions.add(pos)
                }
            }
        }
        //This is a placeholder implementation
        return variantPositions
    }

    fun processSingleGVCFForVariantsSimple(gvcfFile: String): Set<Position> {
        val variantPositions = mutableSetOf<Position>()
        bufferedReader(gvcfFile).use { reader ->
            var line = reader.readLine()
            while(line != null) {
                //Metadata lines start with ##
                if(line.startsWith("#")) {
                    line = reader.readLine()
                    continue
                }
                //otherwise it's a variant line.  We need to parse it down, check to see if its a SNP and then add it to the set
                val stringParsed = line.split("\t")
                val contig = stringParsed[0]
                val pos = stringParsed[1].toInt()
                val ref = stringParsed[3]
                val alts = stringParsed[4].split(",").filter { it != "<NON_REF>" } //Alt can be multiple alleles

                if(alts.isEmpty()) {
                    line = reader.readLine()
                    continue
                }
                //Check to see if it's a SNP
                val isSNP = ref.length == 1 && alts.all { it.length == 1  }
                if(isSNP) {
                    variantPositions.add(Position(contig, pos))
                }
                line = reader.readLine()
            }
        }
        return variantPositions
    }

    fun processAllGVCFsMultithread(gvcfFiles : List<String>, numThreads: Int): Set<Position> {
        //Implement function to process all gvcf files in parallel
        val inputFileChannel = Channel<String>(100)
        val outputChannel = Channel<Set<Position>>(100)
        val variantSet = mutableSetOf<Position>()
        runBlocking {
            for(inputFile in gvcfFiles) {
                launch(Dispatchers.IO) {
                    //Add each gvcf file to the processing queue
                    inputFileChannel.send(inputFile)
                }
            }
            inputFileChannel.close()

            val processGVCFsJobList: MutableList<Job> = mutableListOf()

            repeat(numThreads) {
                processGVCFsJobList.add(launch(Dispatchers.IO) {
                    processSingleGVCFMultithread(inputFileChannel, outputChannel)
                })
            }
            //Launch a job to process the output channel

            val processPositionsJob = launch(Dispatchers.Default) {
                processPositions(outputChannel, variantSet)
            }

            //Wait for all gvcf processing jobs to finish
            processGVCFsJobList.joinAll()
            outputChannel.close()
            //Wait for the position processing job to finish
            processPositionsJob.join()
        }
        return variantSet
    }

    suspend fun processSingleGVCFMultithread(inputFileChannel: Channel<String>, outputChannel: SendChannel<Set<Position>>) {
        for(inputFile in inputFileChannel) {

            val variantPositions = mutableSetOf<Position>()
            myLogger.info("Processing file: $inputFile")
            var counter = 0
            bufferedReader(inputFile).use { reader ->

                val startTime = System.nanoTime()
                var line = reader.readLine()
                while(line != null) {
//                    counter++
//                    if(counter%10000 == 0) {
//                        myLogger.info("Processed $counter lines in file: $inputFile")
//                    }
                    //Metadata lines start with ##
                    if(line.startsWith("#")) {
                        line = reader.readLine()
                        continue
                    }
                    //otherwise it's a variant line.  We need to parse it down, check to see if its a SNP and then add it to the set
                    val stringParsed = line.split("\t")
                    val contig = stringParsed[0]
                    val pos = stringParsed[1].toInt()
                    val ref = stringParsed[3]
                    val alts = stringParsed[4].split(",").filter { it != "<NON_REF>" } //Alt can be multiple alleles

                    if(alts.isEmpty()) {
                        line = reader.readLine()
                        continue
                    }
                    //Check to see if it's a SNP
                    val isSNP = ref.length == 1 && alts.all { it.length == 1  }
                    if(isSNP) {
//                        outputChannel.send(Position(contig, pos))
                        variantPositions.add(Position(contig, pos))
                    }
                    line = reader.readLine()
                }
                outputChannel.send(variantPositions)
                val endTime = System.nanoTime()
                myLogger.info("Finished processing file: $inputFile in ${(endTime - startTime) / 1E9} seconds")
            }
        }
    }

//    suspend fun processPositions(positionChannel: Channel<Position>, variantSet: MutableSet<Position>) {
//        for(position in positionChannel) {
//            variantSet.add(position)
//        }
//    //        positionChannel.consumeEach { position ->
////            variantSet.add(position)
////        }
//    }

    suspend fun processPositions(positionChannel: Channel<Set<Position>>, variantSet: MutableSet<Position>) {
        for(positions in positionChannel) {
            variantSet.addAll(positions)
        }
        //        positionChannel.consumeEach { position ->
//            variantSet.add(position)
//        }
    }



    //Now need a function to go through each position, get the snp value from the gvcf and build a merged Variant context to export
    fun buildMergedVariantContexts(variantPosSet: Set<Position>, gvcfPath: String, outputFile : String) {
        //loop through each gvcf file and build a VCFFileReader for each
        val fileReaders = File(gvcfPath).walkTopDown().filter{
            it.isFile && (it.name.endsWith(".gvcf") || it.name.endsWith(".gvcf.gz") || it.name.endsWith("g.vcf") || it.name.endsWith(".g.vcf.gz"))
        }.map { Pair(it.nameWithoutExtension, VCFFileReader(it,true)) }

        val sampleNames = fileReaders.map { it.first }.toList() //Assuming all files have the same sample name for now

        val sortedPositions = variantPosSet.sorted()
        VariantContextWriterBuilder()
            .unsetOption(Options.INDEX_ON_THE_FLY)
            .setOutputFile(File(outputFile))
            .setOutputFileType(VariantContextWriterBuilder.OutputType.VCF)
            .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
            .build().use { writer ->

                val header = createGenericHeader(sampleNames,emptySet<VCFHeaderLine>())

                writer.writeHeader(header)

                val reportingStepSize = sortedPositions.size/100

                for((index,position) in sortedPositions.withIndex()) {
                    if(index % reportingStepSize == 0) {
                        myLogger.info("Processing position ${index} / ${sortedPositions.size}")
                    }
                    val variants = extractOutVariantsForPosition(fileReaders, position)

                    //loop through the variants and merge them together
                    val mergedVariantContext = mergeGVCFContexts(variants)

                    writer.add(mergedVariantContext)
                }
            }


        //TODO implement
    }

    private fun extractOutVariantsForPosition(
        fileReaders: Sequence<Pair<String, VCFFileReader>>,
        position: Position
    ): List<VariantContext> {
        val variants = fileReaders.mapNotNull { (name, fileReader) ->
            val variantContextIterator = fileReader.query(position.contig, position.position, position.position)
            val vcList = mutableListOf<VariantContext>()
            while (variantContextIterator.hasNext()) {
                val vc = variantContextIterator.next()
                vcList.add(vc)
            }
            if (vcList.isEmpty()) {
                myLogger.info("Found ${vcList.size} variants in ${vcList.size} file at position $position")
                null
            } else if (vcList.size >= 1) {
                myLogger.info("Found ${vcList.size} variants in ${vcList.size} file at position $position")
                vcList.first()
            } else {
                vcList.first()
            }
        }

        return variants.toList()
    }

    fun mergeGVCFContexts(variants: List<VariantContext>): VariantContext {
        //need to determine the ref and alt alleles across all variants
        //Only the SNPs should have the correct reference allele as the REFBLOCKs will only have one allele

        var refAllele = ""
        val altAlleles = mutableSetOf<String>()

        val vcb = VariantContextBuilder()

        val refList = mutableListOf<String>()
        val genotypeList = mutableListOf<Genotype>()

        for(variant in variants) {
            val variantRef = variant.reference.baseString
            val variantAlts = variant.alternateAlleles.map { it.baseString }.filter { it != "<NON_REF>" }
            //Check if it's a SNP
            val vcType = getVariantType(variant)
            val genotypeBuilder = GenotypeBuilder()

            when(vcType) {
                "SNP" -> {
                    refAllele = variantRef
                    altAlleles.addAll(variantAlts)
                    genotypeBuilder.alleles(variant.genotypes.first().alleles)
                    genotypeBuilder.name(variant.genotypes.first().sampleName)
                    genotypeList.add(genotypeBuilder.make())
                }
                "REFBLOCK" -> {
                    //Do nothing for now
                    refList.add(variant.genotypes.first().sampleName) //Set get a list of names that we need to set to the ref allele later
                }
                "INSERTION" -> {
                    //set the alt allele to be <INS>
                    altAlleles.add(Allele.SV_SIMPLE_INS.baseString)
                    genotypeBuilder.alleles(listOf(Allele.SV_SIMPLE_INS))
                    genotypeBuilder.name(variant.genotypes.first().sampleName)
                    genotypeList.add(genotypeBuilder.make())
                }
                "DELETION" -> {
                    //set the alt allele to be <DEL>
                    altAlleles.add(Allele.SV_SIMPLE_DEL.baseString)
                    genotypeBuilder.alleles(listOf(Allele.SV_SIMPLE_DEL))
                    genotypeBuilder.name(variant.genotypes.first().sampleName)
                    genotypeList.add(genotypeBuilder.make())
                }
                else -> {
                    //TODO figure out complex handling
                }
            }
        }

        //Need to build ref genotypes
        for(refName in refList) {
            val genotypeBuilder = GenotypeBuilder()
            val refAlleleObj = Allele.create(refAllele, true)
            genotypeBuilder.alleles(listOf(refAlleleObj, refAlleleObj))
            genotypeBuilder.name(refName)
            genotypeList.add(genotypeBuilder.make())
        }

        return vcb.chr(variants.first().contig)
            .start(variants.first().start.toLong())
            .stop(variants.first().end.toLong())
            .alleles(listOf(Allele.create(refAllele, true)) + altAlleles.map { Allele.create(it, false) })
            .genotypes(genotypeList).make()

    }

    fun getVariantType(variant: VariantContext): String {
        val refBases = variant.reference.baseString
        val altBases = variant.alternateAlleles.map { it.baseString }.filter { it != "<NON_REF>" }

        return when {
            altBases.isEmpty() -> "REFBLOCK"
            refBases.length == 1 && altBases.all { it.length == 1 } -> "SNP"
            refBases.length < altBases.maxOf { it.length } -> "INSERTION"
            refBases.length > altBases.minOf { it.length } -> "DELETION"
            else -> "COMPLEX"
        }
    }


}