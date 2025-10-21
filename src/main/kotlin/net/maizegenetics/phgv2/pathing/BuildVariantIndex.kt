package net.maizegenetics.phgv2.pathing

import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.options.default
import com.github.ajalt.clikt.parameters.options.option
import com.github.ajalt.clikt.parameters.options.required
import com.google.common.collect.Range
import com.google.common.collect.RangeMap
import com.google.common.collect.TreeRangeMap
import net.maizegenetics.phgv2.api.HaplotypeGraph
import net.maizegenetics.phgv2.utils.Position
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

    override fun run() {
        myLogger.info("Building Variant Index")
        buildVariantIndex(dbPath, indexFile)
    }

    fun buildVariantIndex(dbPath: String, indexFile: String) {
        //Implement function to build the variant index from the tileDB
        //Build the phg mapping for SampleName + Position Map -> HaplotypeID
        val samplePositionToHapIdMap = buildSamplePositionToHapIdMap(dbPath)

        //Now we need to walk through the gvcf files and build a SNP list
        val variants = buildVariantSet(dbPath)
        //Then we can do a 2nd pass through the gvcf files and build the index using those
    }

    /**
     * Function to build a lookup map that will allow us to quickly find the HaplotypeID for a given SampleName and Position
     */
    fun buildSamplePositionToHapIdMap(dbPath: String): Map<String, RangeMap<Position, String>> {
        //Build a graph from the tileDB at dbPath
        val graph = HaplotypeGraph("$dbPath/hvcf_dataset")

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

    fun buildVariantSet(dbPath: String): Set<Position> {
        //Build a graph from the tileDB at dbPath
        val gvcfPath = "$dbPath/gvcf_dataset"

        //walk through all files in the gvcfPath and build a set of all variants
        File(gvcfPath).walkTopDown().filter{
            it.isFile && (it.name.endsWith(".gvcf") || it.name.endsWith(".gvcf.gz") || it.name.endsWith("g.vcf") || it.name.endsWith(".g.vcf.gz"))
        }.forEach { file ->
            myLogger.info("Processing file: ${file.absolutePath}")
        }

        val graph = HaplotypeGraph("$dbPath/gvcf_dataset")

        val refRanges = graph.ranges()
        val variantSet = mutableSetOf<Position>()

        for(refRange in refRanges) {
            val variantsInRange = graph.variantsInRange(refRange)
            variantSet.addAll(variantsInRange)
        }

        return variantSet
    }
}