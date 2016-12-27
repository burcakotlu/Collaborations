/**
 * 
 */
package hacettepe.lgmd;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import jaxbxjctool.NCBIEutils;

import org.apache.log4j.Logger;

import remap.Remap;
import rsat.GenerationofSequencesandMatricesforSNPs;
import rsat.PositionFrequencyAndLogoMatrices;
import rsat.SNPInformation;
import rsat.UserDefinedObservedAlleles;
import ui.GlanetRunner;

import common.Commons;

import enumtypes.CommandLineArguments;

/**
 * @author Burçak Otlu
 * @date Dec 23, 2016
 * @project Glanet 
 *
 */
public class GenerationofSequencesForSNPsandMatricesforAllTFs {
	
	final static Logger logger = Logger.getLogger(GenerationofSequencesForSNPsandMatricesforAllTFs.class);


	// TF starts
	public static void readGivenSNPsAndWriteSequencesandMatricesForAllTFs(
			Map<String, String> chrName2RefSeqIdforLatestAssemblyReturnedByNCBIEutilsMap, 
			String forRSAFolder,
			Map<String,String> snpKey_chrName_1BasedStart_1BasedEnd_2_ObservedAllelesMap, 
			Map<String, String> tfName2PfmMatrices,
			Map<String, String> tfName2LogoMatrices, 
			String enrichmentType) {
		
		String chrNameWithPreceedingChr = null;
		String chrNameWithoutPreceedingChr = null;

		int snpOneBasedStart;
		int snpOneBasedEnd;

		String tfName;
		String snpKey_chrName_1BasedStart_1BasedEnd = null;
		String slashSeparatedUserDefinedObservedAlleles = null;
		List<String> observedAllelesList = null;
				
		int indexofFirstUnderscore = -1;
		int indexofSecondUnderscore = -1;
		
		String directoryBase = Commons.TF_PFM_AND_LOGO_Matrices + System.getProperty( "file.separator");
		
		SNPInformation snpInformation = null;
		String snpDirectory;		
		String fastaFile;
		String snpReferenceSequence = null;
		int alteredSequenceCount;
			
		// 10 March 2014
		// Pay attention, there can be more than two observed alleles such as
		// A\tG\tT\t-\tACG
		// Pay attention, for the same chrName and ChrPosition there can be more than one rsIDs
		// Therefore each rsInformation can have observedAlleles String. It is rare but it is possible.
		
		/****************************************************************************************/
		/********************* Write TF PFM and Logo Matrices starts*****************************/
		/****************************************************************************************/			
		for( Map.Entry<String, String> pfmEntry : tfName2PfmMatrices.entrySet()){
				tfName = pfmEntry.getKey();
				
				PositionFrequencyAndLogoMatrices.writeMatrixFile( 
						forRSAFolder, 
						directoryBase, 
						tfName,
						Commons.PFM_MATRICES + Commons.UNDERSCORE + tfName, 
						pfmEntry.getValue());
		}// End of for PFM

		for( Map.Entry<String, String> logoEntry : tfName2LogoMatrices.entrySet()){
				tfName = logoEntry.getKey();
				
				PositionFrequencyAndLogoMatrices.writeMatrixFile(
						forRSAFolder, 
						directoryBase, 
						tfName,
						Commons.LOGO_MATRICES + Commons.UNDERSCORE + tfName,
						logoEntry.getValue());
		
		}// End of for LOGO
		/****************************************************************************************/
		/********************* Write TF PFM and Logo Matrices ends*******************************/
		/****************************************************************************************/			
	
		for(Map.Entry<String, String> entry: snpKey_chrName_1BasedStart_1BasedEnd_2_ObservedAllelesMap.entrySet()){
			
			snpKey_chrName_1BasedStart_1BasedEnd = entry.getKey();
			slashSeparatedUserDefinedObservedAlleles = entry.getValue();
			
			indexofFirstUnderscore = snpKey_chrName_1BasedStart_1BasedEnd.indexOf('_');
			indexofSecondUnderscore = snpKey_chrName_1BasedStart_1BasedEnd.indexOf('_',indexofFirstUnderscore+1);
			
			chrNameWithPreceedingChr = snpKey_chrName_1BasedStart_1BasedEnd.substring(0, indexofFirstUnderscore);
			chrNameWithoutPreceedingChr = chrNameWithPreceedingChr.substring( 3);

			snpOneBasedStart = Integer.parseInt(snpKey_chrName_1BasedStart_1BasedEnd.substring(indexofFirstUnderscore + 1, indexofSecondUnderscore));
			snpOneBasedEnd = Integer.parseInt(snpKey_chrName_1BasedStart_1BasedEnd.substring(indexofSecondUnderscore+1));
			
			snpDirectory = forRSAFolder + Commons.SNPs + System.getProperty("file.separator") + snpKey_chrName_1BasedStart_1BasedEnd + System.getProperty( "file.separator");;
			
			

			snpInformation = new SNPInformation();
			snpInformation.setChrNameWithoutPreceedingChr(chrNameWithoutPreceedingChr);
			snpInformation.setOneBasedStart(snpOneBasedStart);
			snpInformation.setOneBasedEnd(snpOneBasedEnd);
			
			
			
			// Get Fasta File for each SNP
			// Get SNP Reference DNA Sequence from fasta file for each SNP
			fastaFile = GenerationofSequencesandMatricesforSNPs.getDNASequence(
					chrNameWithoutPreceedingChr,
					snpOneBasedStart - Commons.NUMBER_OF_BASES_BEFORE_SNP_POSITION,
					snpOneBasedEnd + Commons.NUMBER_OF_BASES_AFTER_SNP_POSITION,
					chrName2RefSeqIdforLatestAssemblyReturnedByNCBIEutilsMap);
			
			snpReferenceSequence = GenerationofSequencesandMatricesforSNPs.getDNASequenceFromFastaFile(fastaFile);
			snpInformation.setSnpReferenceSequence(snpReferenceSequence);

			
			if(snpReferenceSequence.isEmpty()){
				//Alert me and Skip this snp
				System.out.println("Snp Reference Sequence is empty explore why: " + "chr" + chrNameWithoutPreceedingChr + "_" + snpOneBasedStart + "_" + snpOneBasedEnd);
				continue;
			}
			
			/*****************************************************************/
			/******** Write SNP Reference DNA Sequence starts ****************/
			/*****************************************************************/
			GenerationofSequencesandMatricesforSNPs.writeSequenceFile( 
					snpDirectory, 
					Commons.SNP_REFERENCE_SEQUENCE + "_" + entry.getKey(),
					Commons.SNP_REFERENCE_SEQUENCE + "_" + entry.getKey(),
					fastaFile);
			/*****************************************************************/
			/******** Write SNP Reference DNA Sequence ends ******************/
			/*****************************************************************/
			
			//Write snp altered sequences
			//Write Observed Alleles
			GenerationofSequencesandMatricesforSNPs.writeObservedAllelesFile(
					snpDirectory,
					Commons.OBSERVED_ALLELES,
					slashSeparatedUserDefinedObservedAlleles);
			
			//Create SNP Altered Sequences			
			observedAllelesList = GenerationofSequencesandMatricesforSNPs.convertSlashSeparatedObservedAllelesIntoAStringList(slashSeparatedUserDefinedObservedAlleles);
			
			GenerationofSequencesandMatricesforSNPs.createSNPAlteredSequencesUsingUserDefinedObservedAlleles(
					snpInformation,
					observedAllelesList,
					chrName2RefSeqIdforLatestAssemblyReturnedByNCBIEutilsMap);
			//Write SNP Altered Sequences starts
			/*******************************************************************/
			/************* Write SNP Altered Sequences starts ******************/
			/*******************************************************************/					
			alteredSequenceCount = 1;
			
			for(Entry<String, String> name2Sequence: snpInformation.getAlteredSequenceName2SequenceMap().entrySet()){
				
				String alteredSequenceName = name2Sequence.getKey();
				String alteredSequence = name2Sequence.getValue();
				
				GenerationofSequencesandMatricesforSNPs.writeSequenceFile(
						snpDirectory,
						Commons.SNP_ALTERED_SEQUENCE + alteredSequenceCount,
						Commons.SNP_ALTERED_SEQUENCE + Commons.UNDERSCORE + alteredSequenceName + Commons.UNDERSCORE + entry.getKey(),
						alteredSequence);
				
				alteredSequenceCount++;
				
				alteredSequenceName = null;
				alteredSequence = null;
				
			}//End of for each alteredSequence entry
			/*******************************************************************/
			/************* Write SNP Altered Sequences ends ********************/
			/*******************************************************************/
			
			//Write TF extended sequences
			/*******************************************************************/
			/************* Write TF Extended Sequence starts *******************/
			/*******************************************************************/					
			//Write only one TF extended sequence file +/- 200 bp around SNP position
			GenerationofSequencesandMatricesforSNPs.writeTFExtendedPeakSequenceFile( 
					snpDirectory, 
					snpOneBasedStart,
					snpOneBasedEnd,
					chrNameWithoutPreceedingChr, 
					chrName2RefSeqIdforLatestAssemblyReturnedByNCBIEutilsMap);
						
			/*******************************************************************/
			/************* Write TF Extended Sequence ends *********************/
			/*******************************************************************/	
						
		}//End of FOR each snp
						
	}
	
	
	public static void callNCBIREMAPAndGenerateTargetCoordinatesInLatestAssembly( 
			String dataFolder,
			String outputFolder, 
			Map<String,String> snpKey_chrName_1BasedStart_1BasedEnd_2_TargetMap,
			Map<String,String> snpKey_chrName_1BasedStart_1BasedEnd_2_ObservedAllelesMap,
			String remapInputFile_0Based_Start_EndExclusive_GRCh37_P13_BED_FILE,
			String latestAssembyNameReturnedByNCBIEutils,
			Map<String, String> assemblyName2RefSeqAssemblyIDMap) {

		
		String forRSA_REMAP_Folder = outputFolder + Commons.FOR_RSA + System.getProperty( "file.separator") + Commons.NCBI_REMAP + System.getProperty( "file.separator");

		//GLANET internal data is in Commons.GRCH37_P13 coordinates
		String sourceReferenceAssemblyID = assemblyName2RefSeqAssemblyIDMap.get(Commons.GRCH37_P13);
		String targetReferenceAssemblyID = assemblyName2RefSeqAssemblyIDMap.get(latestAssembyNameReturnedByNCBIEutils);;
			
		String merge = Commons.NCBI_REMAP_API_MERGE_FRAGMENTS_DEFAULT_ON;
		String allowMultipleLocation = Commons.NCBI_REMAP_API_ALLOW_MULTIPLE_LOCATIONS_TO_BE_RETURNED_DEFAULT_ON;
		double minimumRatioOfBasesThatMustBeRemapped = Commons.NCBI_REMAP_API_MINIMUM_RATIO_OF_BASES_THAT_MUST_BE_REMAPPED_DEFAULT_0_POINT_5_;
		double maximumRatioForDifferenceBetweenSourceLengtheAndTargetLength = Commons.NCBI_REMAP_API_MAXIMUM_RATIO_FOR_DIFFERENCE_BETWEEN_SOURCE_LENGTH_AND_TARGET_LENGTH_DEFAULT_2;

		String inputFormat = Commons.BED;

		if( GlanetRunner.shouldLog())logger.info( "******************************************************************************");

		Remap.remap(
				dataFolder, 
				sourceReferenceAssemblyID, 
				targetReferenceAssemblyID,
				forRSA_REMAP_Folder + remapInputFile_0Based_Start_EndExclusive_GRCh37_P13_BED_FILE,
				forRSA_REMAP_Folder + Commons.REMAP_DUMMY_OUTPUT_FILE,
				forRSA_REMAP_Folder + Commons.REMAP_REPORT_CHRNAME_1Based_START_END_XLS_FILE,
				forRSA_REMAP_Folder + Commons.REMAP_DUMMY_GENOME_WORKBENCH_PROJECT_FILE, 
				merge, 
				allowMultipleLocation,
				minimumRatioOfBasesThatMustBeRemapped, 
				maximumRatioForDifferenceBetweenSourceLengtheAndTargetLength,
				inputFormat, 
				Commons.REMAP_FROM_GRCh37p13_TO_LATEST_ASSEMBLY_RETURNED_BY_NCBIEUTILS_FOR_RSA);

		Remap.fillConversionMap(
				forRSA_REMAP_Folder, 
				Commons.REMAP_REPORT_CHRNAME_1Based_START_END_XLS_FILE,
				snpKey_chrName_1BasedStart_1BasedEnd_2_TargetMap,
				Commons.UNDERSCORE);
		
		if( GlanetRunner.shouldLog())logger.info( "******************************************************************************");

	}



	
	public static void main(String[] args) {
				
		String glanetFolder = args[CommandLineArguments.GlanetFolder.value()];

		// jobName starts
		String jobName = args[CommandLineArguments.JobName.value()].trim();
		if( jobName.isEmpty()){
			jobName = Commons.NO_NAME;
		}
		// jobName ends

		String dataFolder = glanetFolder + Commons.DATA + System.getProperty( "file.separator");
		String outputFolder = args[CommandLineArguments.OutputFolder.value()];

		//String givenInputDataFolder = outputFolder + Commons.GIVENINPUTDATA + System.getProperty( "file.separator");
		String forRSAFolder = outputFolder + Commons.FOR_RSA + System.getProperty( "file.separator");
		
		
		//Left here
		//Which input file should I use?
		String givenInputFileName = args[CommandLineArguments.InputFileNameWithFolder.value()];

		//Uplift starts
		//Check the assembly of the given input file hg19
		//If it is not equal to the latest assembly hg38
		//Then uplift hg19 into hg38.
		Map<String,String> snpKey_chrName_1BasedStart_1BasedEnd_2_ObservedAllelesMap = new HashMap<String,String>();

		UserDefinedObservedAlleles.fillSNPKey2ObservedAllelesMap(
				givenInputFileName,
				snpKey_chrName_1BasedStart_1BasedEnd_2_ObservedAllelesMap);
		
		Map<String,String> snpKey_chrName_1BasedStart_1BasedEnd_2_TargetMap = new HashMap<String,String>();
		
		Map<String,String> target_snpKey_chrName_1BasedStart_1BasedEnd_2_ObservedAllelesMap = new HashMap<String,String>();

		//Generate remap input file
		//Fill keys of source_1BasedStart_1BasedEnd2TargetMap
		Remap.generateREMAPInputFile(
				forRSAFolder,
				snpKey_chrName_1BasedStart_1BasedEnd_2_ObservedAllelesMap, 
				snpKey_chrName_1BasedStart_1BasedEnd_2_TargetMap,
				Commons.REMAP_INPUT_FILE_0BASED_START_ENDEXCLUSIVE_GRCH37_P13_COORDINATES_BED_FILE);
		
		String latestAssemblyNameReturnedByNCBIEutils = null;
		Map<String, String> assemblyName2RefSeqAssemblyIDMap  = null;
		
	
		//Means that there are no source genomic loci to be remapped to target genomic loci.
		if (!snpKey_chrName_1BasedStart_1BasedEnd_2_TargetMap.isEmpty()){
			
			/***************************************************************************************/
			/***************************************************************************************/
			/***************************************************************************************/
			latestAssemblyNameReturnedByNCBIEutils = NCBIEutils.getLatestAssemblyNameReturnedByNCBIEutils();
			/***************************************************************************************/
			/***************************************************************************************/
			/***************************************************************************************/

			/***************************************************************************************/
			/***************************************************************************************/
			/***************************************************************************************/
			assemblyName2RefSeqAssemblyIDMap = new HashMap<String, String>();
			
			Remap.remap_show_batches(dataFolder, Commons.NCBI_REMAP_API_SUPPORTED_ASSEMBLIES_FILE);
			
			Remap.fillAssemblyName2RefSeqAssemblyIDMap(
					dataFolder, 
					Commons.NCBI_REMAP_API_SUPPORTED_ASSEMBLIES_FILE,
					assemblyName2RefSeqAssemblyIDMap);
			/***************************************************************************************/
			/***************************************************************************************/
			/***************************************************************************************/

			callNCBIREMAPAndGenerateTargetCoordinatesInLatestAssembly(
					dataFolder, 
					outputFolder,
					snpKey_chrName_1BasedStart_1BasedEnd_2_TargetMap, 
					snpKey_chrName_1BasedStart_1BasedEnd_2_ObservedAllelesMap,
					Commons.REMAP_INPUT_FILE_0BASED_START_ENDEXCLUSIVE_GRCH37_P13_COORDINATES_BED_FILE,
					latestAssemblyNameReturnedByNCBIEutils,
					assemblyName2RefSeqAssemblyIDMap);
			
		}//End of IF
		//Uplift ends 
		
	
		//Fill target_snpKey_chrName_1BasedStart_1BasedEnd_2_ObservedAllelesMap
		UserDefinedObservedAlleles.fillTargetSNPKey2ObservedAllelesMap(
				snpKey_chrName_1BasedStart_1BasedEnd_2_ObservedAllelesMap,
				snpKey_chrName_1BasedStart_1BasedEnd_2_TargetMap,
				target_snpKey_chrName_1BasedStart_1BasedEnd_2_ObservedAllelesMap);
		

		// pfm matrices
		//TODO Can you also use ENCODE motifs?
		//String encodeMotifsInputFileName = Commons.ENCODE_MOTIFS;
		String jasparCoreInputFileName = Commons.JASPAR_CORE;


		// Example Data
		// 7 NC_000007.13 GRCh37
		// Chromosome 7 CM000669.2 = NC_000007.14 0 GRCh37
		Map<String, String> chrName2RefSeqIdforLatestAssemblyReturnedByNCBIEutilsMap = new HashMap<String, String>();			
		
//		/***************************************************************************************/
//		/***********************************Part1 starts****************************************/
//		/***************************************************************************************/
//		String latestAssemblyNameReturnedByNCBIEutils = null;
//		latestAssemblyNameReturnedByNCBIEutils = NCBIEutils.getLatestAssemblyNameReturnedByNCBIEutils();
//		/***************************************************************************************/
//		/***********************************Part1 ends******************************************/
//		/***************************************************************************************/
		
		
//		/***************************************************************************************/
//		/***********************************Part2 starts****************************************/
//		/***************************************************************************************/
//		Map<String, String> assemblyName2RefSeqAssemblyIDMap = new HashMap<String, String>();
//		
//		//No need to call this again and again.
//		//Since it is already called in GenerationofAllTFAnnotationsFileInGRCh37p13AndInLatestAssembly class
//		//Remap.remap_show_batches(dataFolder, Commons.NCBI_REMAP_API_SUPPORTED_ASSEMBLIES_FILE);
//		
//		Remap.fillAssemblyName2RefSeqAssemblyIDMap(
//				dataFolder, 
//				Commons.NCBI_REMAP_API_SUPPORTED_ASSEMBLIES_FILE,
//				assemblyName2RefSeqAssemblyIDMap);
//		/***************************************************************************************/
//		/***********************************Part2 ends******************************************/
//		/***************************************************************************************/
	
		/***************************************************************************************/
		/***********************************Part3 starts****************************************/
		/***************************************************************************************/
		String refSeqAssemblyID = NCBIEutils.getRefSeqAssemblyID(latestAssemblyNameReturnedByNCBIEutils, assemblyName2RefSeqAssemblyIDMap);
		/***************************************************************************************/
		/***********************************Part3 ends******************************************/
		/***************************************************************************************/
						
		/***************************************************************************************/
		/***********************************Part4 starts****************************************/
		/***************************************************************************************/
		// Download from  ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/RefSeqAssemblyID.assembly.txt
		String assemblyReportFileName = Commons.ASSEMBLY_REPORTS +  refSeqAssemblyID + Commons.ASSEMBLY_REPORTS_FILE_EXTENSION ;
		NCBIEutils.getAssemblyReport(refSeqAssemblyID, dataFolder, assemblyReportFileName);
		/***************************************************************************************/
		/***********************************Part4 ends******************************************/
		/***************************************************************************************/
						
		/***************************************************************************************/
		/***********************************Part5 starts****************************************/
		/***************************************************************************************/
		NCBIEutils.fillChrName2RefSeqIDMap(dataFolder, assemblyReportFileName, chrName2RefSeqIdforLatestAssemblyReturnedByNCBIEutilsMap);
		/***************************************************************************************/
		/***********************************Part5 ends******************************************/
		/***************************************************************************************/

		// Construct pfm matrices from encode-motif.txt file
		// A tf can have more than one pfm matrices
		// Take the transpose of given matrices in encode-motif.txt
		// Write the matrices in tab format for RSAT tool
		Map<String, String> tfName2PfmMatrices = new HashMap<String, String>();

		Map<String, String> tfName2LogoMatrices = new HashMap<String, String>();

		// Construct position frequency matrices from Encode Motifs
		//PositionFrequencyAndLogoMatrices.constructPfmMatricesfromEncodeMotifs( dataFolder, encodeMotifsInputFileName, tfName2PfmMatrices);

		// Construct logo matrices from Encode Motifs
		//PositionFrequencyAndLogoMatrices.constructLogoMatricesfromEncodeMotifs( dataFolder, encodeMotifsInputFileName, tfName2LogoMatrices);

		
		// Construct position frequency matrices from Jaspar Core
		// Construct logo matrices from Jaspar Core
		PositionFrequencyAndLogoMatrices.constructPfmMatricesandLogoMatricesfromJasparCoreWithSpecificTFName(
				dataFolder, 
				jasparCoreInputFileName, 
				tfName2PfmMatrices,
				tfName2LogoMatrices);
		
		
			
		// Generate DNA Sequences for all given SNPs and TF Matrices for all TF Elements
		//update this method
		//left here
		readGivenSNPsAndWriteSequencesandMatricesForAllTFs(
				chrName2RefSeqIdforLatestAssemblyReturnedByNCBIEutilsMap, 
				forRSAFolder,
				target_snpKey_chrName_1BasedStart_1BasedEnd_2_ObservedAllelesMap, 
				tfName2PfmMatrices, 
				tfName2LogoMatrices,
				Commons.TF);


		

	}

}
