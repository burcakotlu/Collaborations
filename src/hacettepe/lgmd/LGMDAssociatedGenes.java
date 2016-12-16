/**
 * 
 */
package hacettepe.lgmd;

import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntObjectHashMap;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import jaxbxjctool.AugmentationofGivenIntervalwithRsIds;
import jaxbxjctool.AugmentationofGivenRsIdwithInformation;
import jaxbxjctool.RsInformation;
import augmentation.humangenes.HumanGenesAugmentation;
import auxiliary.FileOperations;
import auxiliary.HG38RefSeqGeneInformation;

/**
 * @author Burçak Otlu
 * @date Nov 30, 2016
 * @project Collaborations 
 *
 *
 * Limb Girdle Muscular Dystrophy Associated Genes
 * These genes are from https://www.omim.org/phenotypicSeries/PS253600
 * 
 */
public class LGMDAssociatedGenes {
	
	
	public static void addGeneSymbols(List<String> LGMDAssociatedGeneSymbolsList,String geneSymbol){
		if (!LGMDAssociatedGeneSymbolsList.contains(geneSymbol)){
			LGMDAssociatedGeneSymbolsList.add(geneSymbol);
		}
	}
	
	public static void readLGMDAssociatedGenes(
			String OMIMLGMDDirectory,
			String LGMDAssociatedGenesFileName,
			List<String> LGMDAssociatedGeneSymbolsList){
		
		FileReader fileReader = null;
		BufferedReader bufferedReader = null;
		
		String strLine = null;
		String genes = null;
		String gene = null;
		
		int indexofFirstTab = -1;
		int indexofSecondTab = -1;
		int indexofThirdTab = -1;
		int indexofFourthTab = -1;
		int indexofFifhTab = -1;
		
		int indexofFormerComma = -1;
		int indexofLatterComma = -1;
		
		
		try {
			fileReader = FileOperations.createFileReader(OMIMLGMDDirectory + LGMDAssociatedGenesFileName);
			bufferedReader = new BufferedReader(fileReader);
						
			while((strLine= bufferedReader.readLine())!=null){
				
				//Example line
				//1p34.1	"Muscular dystrophy-dystroglycanopathy (limb-girdle), type C, 3"	3	613157	"POMGNT1, MEB, MDDGA3, MDDGB3, MDDGC3, RP76"	606822

				indexofFirstTab= strLine.indexOf("\t");
				indexofSecondTab= strLine.indexOf("\t",indexofFirstTab+1);
				indexofThirdTab= strLine.indexOf("\t",indexofSecondTab+1);
				indexofFourthTab= strLine.indexOf("\t",indexofThirdTab+1);
				indexofFifhTab= strLine.indexOf("\t",indexofFourthTab+1);
				
				//example genes 
				//"POMGNT1, MEB, MDDGA3, MDDGB3, MDDGC3, RP76"
				genes = strLine.substring(indexofFourthTab+1, indexofFifhTab);
				
				//Get rid of double quotas at the beginning and at the end
				genes = genes.substring(1, genes.length()-1);
				
				indexofFormerComma = genes.indexOf(",");
				
				//First gene
				if (indexofFormerComma>0){
					gene = genes.substring(0,indexofFormerComma).trim().toUpperCase();
					indexofLatterComma = genes.indexOf(",",indexofFormerComma+1);	
					
					if (!LGMDAssociatedGeneSymbolsList.contains(gene)){
						LGMDAssociatedGeneSymbolsList.add(gene);
					}
					
				}else{
					//There is only one gene
					gene = genes.trim().toUpperCase();
					if (!LGMDAssociatedGeneSymbolsList.contains(gene)){
						LGMDAssociatedGeneSymbolsList.add(gene);
					}
				}
				
				//Middle gene or genes
				while(indexofLatterComma>0){
					gene = genes.substring(indexofFormerComma+1,indexofLatterComma).trim().toUpperCase();
					if (!LGMDAssociatedGeneSymbolsList.contains(gene)){
						LGMDAssociatedGeneSymbolsList.add(gene);
					}

					indexofFormerComma = indexofLatterComma;
					indexofLatterComma = genes.indexOf(",",indexofFormerComma+1);										
				}//End of while
				
				//Last gene
				if (indexofFormerComma>0){
					gene = genes.substring(indexofFormerComma+1).trim().toUpperCase();
					if (!LGMDAssociatedGeneSymbolsList.contains(gene)){
						LGMDAssociatedGeneSymbolsList.add(gene);
					}

				}

			}//End of while
			
			//gene symbol are between 4th tab and 5th tab
			
			//Close
			bufferedReader.close();
			
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static void getLGMDAssociatedGenesEntrezIDs(
			String dataFolder,
			List<String> LGMDAssociatedGenes,
			List<Integer> LGMDAssociatedGeneIdsList){
		
		String geneSymbol = null;
		List<Integer> geneIDList = null;
		int geneID= -1;
		
		TIntObjectMap<String> geneEntrezId2GeneOfficialSymbolMap = new TIntObjectHashMap<String>();
		Map<String,List<Integer>> geneSymbol2ListofGeneIdMap = new HashMap<String,List<Integer>>();
		
		HumanGenesAugmentation.fillGeneId2GeneHugoSymbolMap(
				dataFolder, 
				geneEntrezId2GeneOfficialSymbolMap);
		
		HumanGenesAugmentation.fillGeneSymbol2ListofGeneIDMap(
				dataFolder,
				geneSymbol2ListofGeneIdMap);
		
		
		for(Iterator<String> itr =LGMDAssociatedGenes.iterator();itr.hasNext();){
			
			geneSymbol = itr.next();
			geneIDList = geneSymbol2ListofGeneIdMap.get(geneSymbol);
			
			if (geneIDList!=null){
				
				//System.out.println(geneSymbol + "\t" + geneIDList.get(0) + "\t" + "contains.");
				
				for(Iterator<Integer> itr2 =geneIDList.iterator();itr2.hasNext();){
					geneID = itr2.next();
					
					if (!LGMDAssociatedGeneIdsList.contains(geneID)){
						LGMDAssociatedGeneIdsList.add(geneID);
					}
					
				}
			
			}
			
		}//End of FOR		
		
	}
	
	public static void checkAnyCommonGenesBetweenLGMDAssociatedGenesAndGLANETAnalysis(
			String folder,
			String fileName,
			String commonsFile,
			List<String> LGMDAssociatedGeneSymbols,
			List<Integer> LGMDAssociatedGeneIds){
				
		FileReader fileReader = null;
		BufferedReader bufferedReader = null;
		String strLine = null;
		
		String geneSymbol = null;
		
		FileWriter fileWriter = null;
		BufferedWriter bufferedWriter = null;
		
		
		//Read File
		//Check whether there is any common geneID and geneSymbol
		
		try {
			fileReader = FileOperations.createFileReader(folder + fileName);
			bufferedReader = new BufferedReader(fileReader);
			
			fileWriter =FileOperations.createFileWriter(folder + commonsFile);
			bufferedWriter = new BufferedWriter(fileWriter);
			
			//Don't Skip header line since every file may not have header line
			
			//Example lines
			//GivenIntervalNumber	GivenInteval	#ofExonOverlaps	ExonOverlapsInformation	#ofIntronOverlaps	IntronOverlapsInformation	#of5p1Overlaps	5p1OverlapsInformation	#of5p2Overlaps	5p2OverlapsInformation	#of5dOverlaps	5dOverlapsInformation	#of3p1Overlaps	3p1OverlapsInformation	#of3p2Overlaps	3p2OverlapsInformation	#of3dOverlaps	3dOverlapsInformation	#ofGeneOverlapsPerGivenInterval
			//1	chr1_14541_14541	1	[653635_WASH7P_EXON_1] 	0		0		0		1	[79501_OR4F5_5D] 	1	[100287102_DDX11L1_3P1] 	2	[102465909_MIR6859-2_3P2] [102466751_MIR6859-1_3P2] 	2	[641702_FAM138F_3D] [645520_FAM138A_3D] 	7
						
			while((strLine=bufferedReader.readLine())!=null){
				
				for(Iterator<String> itr=LGMDAssociatedGeneSymbols.iterator();itr.hasNext();){
					geneSymbol = itr.next();					
					if(strLine.contains("_" +geneSymbol + "_")){
						bufferedWriter.write(geneSymbol + " --> "+ strLine + System.getProperty("line.separator"));
					}
				}//End of for each gene symbol
				
//				for(Iterator<Integer> itr=LGMDAssociatedGeneIds.iterator();itr.hasNext();){
//					geneID = itr.next();					
//					if(strLine.contains("[" +geneID + "_")){
//						System.out.println(geneID + " --> "+ strLine);
//					}
//				}//End of for each gene id
								
			}//End of reading file
			
			//Close
			bufferedReader.close();
			bufferedWriter.close();
		
		} catch (IOException e) {
			e.printStackTrace();
		}				
	}
	
	
	public static void readIntervalsAndWriteSNPs(
			String LGMDRelatedGenesFolder,
			String LGMDRelatedGenesIntervalsFile,
			String SNPsForLGMDRelatedGenesIntervalsFile){
		
		FileReader fileReader = null;
		BufferedReader bufferedReader = null;
		
		FileWriter fileWriter = null;
		BufferedWriter bufferedWriter = null;
		
		String strLine = null;
		int indexofFirstTab = -1;
		int indexofSecondTab = -1;
		
		String chrName= null;
		int start= -1;
		int end = -1;
		
		AugmentationofGivenIntervalwithRsIds augmentationofGivenIntervalwithRsIds  = null;
		AugmentationofGivenRsIdwithInformation augmentationofGivenRsIdwithInformation = null;
		
	
		List<Integer> rsIDList = null;
		
		Integer rsID = -1;
		
		RsInformation rsInformation = null;
				
		try {
			
			augmentationofGivenIntervalwithRsIds = new AugmentationofGivenIntervalwithRsIds();
			augmentationofGivenRsIdwithInformation = new AugmentationofGivenRsIdwithInformation();
		
			fileReader = FileOperations.createFileReader(LGMDRelatedGenesFolder + LGMDRelatedGenesIntervalsFile);
			bufferedReader = new BufferedReader(fileReader);
			
			fileWriter = FileOperations.createFileWriter(LGMDRelatedGenesFolder + SNPsForLGMDRelatedGenesIntervalsFile);
			bufferedWriter  = new BufferedWriter(fileWriter);
	
			while((strLine = bufferedReader.readLine())!=null){
				
				indexofFirstTab = strLine.indexOf('\t');
				indexofSecondTab = strLine.indexOf('\t', indexofFirstTab+1);
				
				chrName = strLine.substring(0, indexofFirstTab);
				start = Integer.parseInt(strLine.substring(indexofFirstTab+1, indexofSecondTab));
				end = Integer.parseInt(strLine.substring(indexofSecondTab+1));
				
				rsIDList = augmentationofGivenIntervalwithRsIds.getRsIdsInAGivenInterval(chrName.substring(3),start,end);
				
				if (rsIDList!=null && !rsIDList.isEmpty()){
					
					//rsIDStringList = new ArrayList<String>();
					
					for(Iterator<Integer> itr= rsIDList.iterator();itr.hasNext();){
						rsID = itr.next();
						 
						rsInformation = augmentationofGivenRsIdwithInformation.getInformationforGivenRsId(rsID.toString());

						if (rsInformation!=null){
							bufferedWriter.write("rs" + rsInformation.getRsId() +  System.getProperty("line.separator"));
							
						}
						
					}//End of FOR
					
//					rsInformationList = augmentationofGivenRsIdwithInformation.getInformationforGivenRsIdList(rsIDStringList,ncbiEutilStatistics);
//					
//					if (rsInformationList==null  || rsInformationList.isEmpty()){
//						System.out.println(rsIDList.size());						
//					}
//					
//					for(Iterator<RsInformation> itr=rsInformationList.iterator();itr.hasNext();){
//						rsInformation = itr.next();
//						bufferedWriter.write("rs" +  rsInformation.getRsId() +  System.getProperty("line.separator"));
//								
//					}//End of FOR
											
				}//End of IF
								
			}//End of WHILE
			
			//Close
			bufferedReader.close();
			bufferedWriter.close();
		
		} catch (Exception e1) {
			e1.printStackTrace();
		}
				
		
	}

	
	
	public static void generateIntervalsForGivenGeneSymbols(
			String LGMDRelatedGenesFolder,
			String LGMDRelatedGenesIntervalsFile,
			String LGMDRelatedGenesSymbolsThatAreFoundFile,
			String LGMDRelatedGenesSymbolsThatAreNotFoundFile,
			List<String> LGMDAssociatedGeneSymbolsList,
			TObjectIntMap<String> geneSymbol2GeneInternalNumberMap,		
			TIntObjectMap<HG38RefSeqGeneInformation> geneInternalNumber2GeneInformationMap,
			int numberofBasesBeforeTSS,
			int numberofBasedAfterTSS){
		
		String geneSymbol = null;
		int geneNumberInternal = -1;
		HG38RefSeqGeneInformation geneInformation = null;
		
		char geneStrand;
		
		FileWriter fileWriter = null;
		BufferedWriter bufferedWriter = null;
		
		FileWriter fileWriterGeneSymbolsThatAreFound = null;
		BufferedWriter bufferedWriterGeneSymbolsThatAreFound = null;
		
		FileWriter fileWriterGeneSymbolsThatAreNotFound = null;
		BufferedWriter bufferedWriterGeneSymbolsThatAreNotFound = null;


		
		int promoterStart = -1;
		int promoterEnd = -1;
		
		try {
			fileWriter = FileOperations.createFileWriter(LGMDRelatedGenesFolder + LGMDRelatedGenesIntervalsFile);
			bufferedWriter = new BufferedWriter(fileWriter);
			
			fileWriterGeneSymbolsThatAreFound = FileOperations.createFileWriter(LGMDRelatedGenesFolder + LGMDRelatedGenesSymbolsThatAreFoundFile);
			bufferedWriterGeneSymbolsThatAreFound = new BufferedWriter(fileWriterGeneSymbolsThatAreFound);

			fileWriterGeneSymbolsThatAreNotFound = FileOperations.createFileWriter(LGMDRelatedGenesFolder + LGMDRelatedGenesSymbolsThatAreNotFoundFile);
			bufferedWriterGeneSymbolsThatAreNotFound = new BufferedWriter(fileWriterGeneSymbolsThatAreNotFound);
			
			for(Iterator<String> itr=LGMDAssociatedGeneSymbolsList.iterator();itr.hasNext();){
				
				geneSymbol = itr.next();
				
				geneNumberInternal = geneSymbol2GeneInternalNumberMap.get(geneSymbol);
				geneInformation = geneInternalNumber2GeneInformationMap.get(geneNumberInternal);
			
				//Found
				if (geneInformation!=null){
					
					geneStrand = geneInformation.getStrand();
					
					//Always from 5' prime to 3' prime
					//At '+' strand gene starts at TxStart (5'prime txStart txEnd 3'prime)
					//At '-' strand gene starts at TxEnd   (3'prime txStart txEnd 5'prime)
					if (geneStrand=='+'){
						promoterStart = geneInformation.getTxStart()-numberofBasesBeforeTSS;
						promoterEnd = geneInformation.getTxStart()+numberofBasedAfterTSS;
						
					}else if (geneStrand=='-'){
						promoterEnd = geneInformation.getTxEnd()+numberofBasesBeforeTSS;
						promoterStart = geneInformation.getTxEnd()-numberofBasedAfterTSS;						
					}
					
					bufferedWriterGeneSymbolsThatAreFound.write(geneInformation.getGeneSymbol() + System.getProperty("line.separator"));
					bufferedWriter.write(geneInformation.getChromName().convertEnumtoString() + "\t" + promoterStart + "\t" + promoterEnd + System.getProperty("line.separator"));				
				}
				//Not Found
				else{
					bufferedWriterGeneSymbolsThatAreNotFound.write(geneSymbol + System.getProperty("line.separator"));
					
				}
							
			}//End of FOR			
			
			//Close
			bufferedWriter.close();
			bufferedWriterGeneSymbolsThatAreFound.close();
			bufferedWriterGeneSymbolsThatAreNotFound.close();
			
		
		} catch (IOException e) {
			e.printStackTrace();
		}

		
	}
	
	
	public static void main(String[] args) {
		
		//Read LGMD Associated gene symbols and get their entrez IDs.
		//Overlap all LGMD variants with RefSeq Genes, check whether they include any LGMD Associated gene
		//Overlap all rare LGMD variants with RefSeq Genes, check whether they include any LGMD Associated genes
		
		
		/*********************************************************************************/
		/*****************Get LGMD related genes******************************************/
		/*********************************************************************************/
		String OMIMLGMDDirectory = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\OMIM\\";		
		//String LGMDAssociatedGenesFileName = "OMIM_LGMD_Associated_Genes.txt";
		
		//only added ANPEP for LAP1 using http://www.genenames.org/
		String LGMDAssociatedGenesFileName_Burcak_Updated = "OMIM_LGMD_Associated_Genes_Burcak_Manually_Updated.txt";
		
		String dataFolder = "C:\\Users\\Burçak\\Google Drive\\Data\\";
		
		List<String> LGMDAssociatedGeneSymbolsList = new ArrayList<String>();
		List<Integer> LGMDAssociatedGeneIdsList = new ArrayList<Integer>();
		
		
		//LGMD has two forms: Autosomal dominant (LGMD1) and the autosomal recessive forms (LGMD2).
		
		//Fill LGMDAssociatedGeneSymbolsList
		readLGMDAssociatedGenes(OMIMLGMDDirectory,LGMDAssociatedGenesFileName_Burcak_Updated,LGMDAssociatedGeneSymbolsList);
		
		//Manually add genes from GulsumKayman email
		addGeneSymbols(LGMDAssociatedGeneSymbolsList,"HOPX");
		addGeneSymbols(LGMDAssociatedGeneSymbolsList,"RASL11B");
		addGeneSymbols(LGMDAssociatedGeneSymbolsList,"CLOCK");
		addGeneSymbols(LGMDAssociatedGeneSymbolsList,"SGCG");
		addGeneSymbols(LGMDAssociatedGeneSymbolsList,"SGCA");
		addGeneSymbols(LGMDAssociatedGeneSymbolsList,"SGCB");
		addGeneSymbols(LGMDAssociatedGeneSymbolsList,"SGCD");
		
		
		//Manually add genes from https://ghr.nlm.nih.gov/condition/limb-girdle-muscular-dystrophy#genes
		addGeneSymbols(LGMDAssociatedGeneSymbolsList,"LMNA");
		addGeneSymbols(LGMDAssociatedGeneSymbolsList,"CAV3");
		addGeneSymbols(LGMDAssociatedGeneSymbolsList,"CAPN3");
		addGeneSymbols(LGMDAssociatedGeneSymbolsList,"DYSF");
		addGeneSymbols(LGMDAssociatedGeneSymbolsList,"TTN");
		
		//Manually added genes from http://www.kegg.jp/dbget-bin/www_bget?ds:H00593
		//KEGG Pathway Disease Code for LGMD is H00593
		addGeneSymbols(LGMDAssociatedGeneSymbolsList,"TTID");
		addGeneSymbols(LGMDAssociatedGeneSymbolsList,"LMNA");
		addGeneSymbols(LGMDAssociatedGeneSymbolsList,"CAV3");
		
				
		//Fill LGMDAssociatedGeneIdsList
		//Pay attention I'm using only gene symbols in the rest of the code.
		getLGMDAssociatedGenesEntrezIDs(dataFolder,LGMDAssociatedGeneSymbolsList,LGMDAssociatedGeneIdsList);
		
		
		//Just for information has nothing to do with the rest of the code.
		System.out.println("Number of gene symbols in the list \t" + LGMDAssociatedGeneSymbolsList.size());
		System.out.println("Number of gene ids in the list \t" + LGMDAssociatedGeneIdsList.size());
		/*********************************************************************************/
		/*****************Get LGMD related genes******************************************/
		/*********************************************************************************/

		
		
//		/*********************************************************************************/
//		/*****************Generate Intervals for LGMD related genes***********************/
//		/*********************************************************************************/
//		//Question: where does the TFs bind in genes?
//		//TFs binds close to gene's promoter regions.
//		String UCSC_GENOME_HG38_REFSEQ_GENES_FILE = dataFolder + Commons.UCSCGENOME_HG38_REFSEQ_GENES_DOWNLOADED_2_DEC_2016;
//		TObjectIntMap<String> geneSymbol2GeneInternalNumberMap = new TObjectIntHashMap<String>();		
//		TIntObjectMap<HG38RefSeqGeneInformation> geneInternalNumber2HG38GeneInformationMap = new TIntObjectHashMap<HG38RefSeqGeneInformation>();	
//		
//		//Fill the maps using HG38 RefSeq genes file
//		//This file only contains geneSymbol and geneRefSeqName. 
//		//It doesn't contains gene entrez id
//		//Therefore we have nothing to do with gene entrez id.
//		HG38_RefSeq_Genes.readUCSCGenomeHG38RefSeqGenes(
//				UCSC_GENOME_HG38_REFSEQ_GENES_FILE,
//				geneSymbol2GeneInternalNumberMap,
//				geneInternalNumber2HG38GeneInformationMap);
//		
//		String LGMDRelatedGenesFolder = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\OMIM\\";
//		//String LGMDRelatedGenesIntervalsFile = "LGMD_Related_Genes_Promoters_1KB_Intervals_In_HG38_Coordinates.txt"; 
//		String LGMDRelatedGenesIntervalsFile = "LGMD_Related_Genes_Promoter_2KB_Intervals_In_HG38_Coordinates.txt"; 
//		String LGMDRelatedGenesSymbolsThatAreFoundFile = "LGMD_Related_Genes_Symbols_That_Are_Found.txt"; 
//		String LGMDRelatedGenesSymbolsThatAreNotFoundFile = "LGMD_Related_Genes_Symbols_That_Are_Not_Found.txt"; 
//
//		generateIntervalsForGivenGeneSymbols(
//				LGMDRelatedGenesFolder,
//				LGMDRelatedGenesIntervalsFile,
//				LGMDRelatedGenesSymbolsThatAreFoundFile,
//				LGMDRelatedGenesSymbolsThatAreNotFoundFile,
//				LGMDAssociatedGeneSymbolsList,
//				geneSymbol2GeneInternalNumberMap,
//				geneInternalNumber2HG38GeneInformationMap,
//				2000,
//				200);
//		/*********************************************************************************/
//		/*****************Generate Intervals for LGMD related genes***********************/
//		/*********************************************************************************/
		
//		/*********************************************************************************/
//		/*****************Find SNPs within provided LGMD related genes intervals**********/
//		/*********************************************************************************/
//		//Read intervals data in hg38
//		//Write SNPs data (their rsIDs)
//		//String SNPsForLGMDRelatedGenesIntervalsFile = "SNPs_for_LGMD_Related_Genes_Promoters_1KB_Intervals.txt"; 
//		String SNPsForLGMDRelatedGenesIntervalsFile = "SNPs_for_LGMD_Related_Genes_Promoter_2KB_Intervals.txt"; 
//		
//		readIntervalsAndWriteSNPs(
//				LGMDRelatedGenesFolder,
//				LGMDRelatedGenesIntervalsFile,
//				SNPsForLGMDRelatedGenesIntervalsFile);
//		/*********************************************************************************/
//		/*****************Find SNPs within provided LGMD related genes intervals**********/
//		/*********************************************************************************/

	
		/*********************************************************************************/
		/*****************Find common genes***********************************************/
		/*********************************************************************************/
		String geneOvelapFolder = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\GLANET_Outputs\\LGMD_ENCODE_GENE_TFKEGG\\Annotation\\Hg19RefSeqGene\\Analysis\\";		
		String overlapAnalysisFile = "Overlap_Analysis_File.txt"; 

		String commonLGMDGenes_For_RareLGMDSNPs_GeneAnnotations_File ="CommonLGMDGenes_For_RareLGMDSNPs_GeneAnnotations.txt";
		checkAnyCommonGenesBetweenLGMDAssociatedGenesAndGLANETAnalysis(
				geneOvelapFolder,
				overlapAnalysisFile,
				commonLGMDGenes_For_RareLGMDSNPs_GeneAnnotations_File,
				LGMDAssociatedGeneSymbolsList,
				LGMDAssociatedGeneIdsList);
		/*********************************************************************************/
		/*****************Find common genes***********************************************/
		/*********************************************************************************/
		
		
		/*********************************************************************************/
		/*****************Find common genes for LGMD given 6085 SNPS**********************/
		/*********************************************************************************/
		String postAnalysisRSAFolder = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\GLANET_Outputs\\LGMD_SNPs_RSA\\RegulatorySequenceAnalysis\\";		
		String postAnalysisRSAFile = "PostAnalysisofRegulatorySequenceAnalysisResults_AugmentedWithGeneAnnotations.txt";
		String commonLGMDGenes_For_RareLGMDSNPs_PostAnalysisRSA_File ="CommonLGMDGenes_For_RareLGMDSNPs_PostAnalysisRSA.txt";
		
		checkAnyCommonGenesBetweenLGMDAssociatedGenesAndGLANETAnalysis(
				postAnalysisRSAFolder,
				postAnalysisRSAFile,
				commonLGMDGenes_For_RareLGMDSNPs_PostAnalysisRSA_File,
				LGMDAssociatedGeneSymbolsList,
				LGMDAssociatedGeneIdsList);
		/*********************************************************************************/
		/*****************Find common genes***********************************************/
		/*********************************************************************************/
		
		
		
		/*********************************************************************************/
		/***Find common genes for 619 SNPs in 2KB promoter region of LGMD related genes***/
		/*********************************************************************************/
		postAnalysisRSAFolder = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\GLANET_Outputs\\LGMD_Related_Genes_SNPs_IN_2KB_RSA\\RegulatorySequenceAnalysis\\";		
		postAnalysisRSAFile = "PostAnalysisofRegulatorySequenceAnalysisResults_AugmentedWithGeneAnnotations.txt";
		commonLGMDGenes_For_RareLGMDSNPs_PostAnalysisRSA_File ="CommonLGMDGenes_For_SNPs_In_2KB_Promoters_of_LGMDRelatedGenes_PostAnalysisRSA.txt";
		
		checkAnyCommonGenesBetweenLGMDAssociatedGenesAndGLANETAnalysis(
				postAnalysisRSAFolder,
				postAnalysisRSAFile,
				commonLGMDGenes_For_RareLGMDSNPs_PostAnalysisRSA_File,
				LGMDAssociatedGeneSymbolsList,
				LGMDAssociatedGeneIdsList);
		/*********************************************************************************/
		/*****************Find common genes***********************************************/
		/*********************************************************************************/


	}

}
