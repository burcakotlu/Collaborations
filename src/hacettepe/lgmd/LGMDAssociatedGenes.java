/**
 * 
 */
package hacettepe.lgmd;

import gnu.trove.map.TIntObjectMap;
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

import augmentation.humangenes.HumanGenesAugmentation;
import auxiliary.FileOperations;

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
					gene = genes.substring(0,indexofFormerComma).trim();
					indexofLatterComma = genes.indexOf(",",indexofFormerComma+1);	
					
					if (!LGMDAssociatedGeneSymbolsList.contains(gene)){
						LGMDAssociatedGeneSymbolsList.add(gene);
					}
					
				}else{
					//There is only one gene
					gene = genes.trim();
					if (!LGMDAssociatedGeneSymbolsList.contains(gene)){
						LGMDAssociatedGeneSymbolsList.add(gene);
					}
				}
				
				//Middle gene or genes
				while(indexofLatterComma>0){
					gene = genes.substring(indexofFormerComma+1,indexofLatterComma).trim();
					if (!LGMDAssociatedGeneSymbolsList.contains(gene)){
						LGMDAssociatedGeneSymbolsList.add(gene);
					}

					indexofFormerComma = indexofLatterComma;
					indexofLatterComma = genes.indexOf(",",indexofFormerComma+1);										
				}//End of while
				
				//Last gene
				if (indexofFormerComma>0){
					gene = genes.substring(indexofFormerComma+1).trim();
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
			List<Integer> LGMDAssociatedGeneIdsList,
			Map<String,Integer> LGMDAssociatedGeneName2EntrezIDMap){
		
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
		Integer geneID = -1;
		
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
	

	
	public static void main(String[] args) {
		
		//Read LGMD Associated gene symbols and get their entrez IDs.
		//Overlap all LGMD variants with RefSeq Genes, check whether they include any LGMD Associated gene
		//Overlap all rare LGMD variants with RefSeq Genes, check whether they include any LGMD Associated genes
		
		String OMIMLGMDDirectory = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\OMIM\\";
		
		String LGMDAssociatedGenesFileName = "OMIM_LGMD_Associated_Genes.txt";
		
		//only added ANPEP for LAP1 using http://www.genenames.org/
		String LGMDAssociatedGenesFileName_Burcak_Updated = "OMIM_LGMD_Associated_Genes_Burcak_Manually_Updated.txt";
		
		String dataFolder = "C:\\Users\\Burçak\\Google Drive\\Data\\";
		
		List<String> LGMDAssociatedGeneSymbolsList = new ArrayList<String>();
		List<Integer> LGMDAssociatedGeneIdsList = new ArrayList<Integer>();
		
		Map<String,Integer> LGMDAssociatedGeneSymbol2GeneIDMap = new HashMap<String,Integer>();
		
		//Fill LGMDAssociatedGeneSymbolsList
		readLGMDAssociatedGenes(OMIMLGMDDirectory,LGMDAssociatedGenesFileName_Burcak_Updated,LGMDAssociatedGeneSymbolsList);
		
		//Pay attention
		//LGMDAssociatedGeneSymbol2GeneIDMap is empty
		getLGMDAssociatedGenesEntrezIDs(dataFolder,LGMDAssociatedGeneSymbolsList,LGMDAssociatedGeneIdsList,LGMDAssociatedGeneSymbol2GeneIDMap);
		
		//LGMDAssociatedGeneSymbols
		//LGMDAssociatedGeneIds
		
		System.out.println("Number of gene symbols in the list \t" + LGMDAssociatedGeneSymbolsList.size());
		System.out.println("Number of gene ids in the list \t" + LGMDAssociatedGeneIdsList.size());
		
		String geneOvelapFolder = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\GeneOverlaps\\";
		String overlapAnalysisFile = "Overlap_Analysis_File.txt"; 
		String commonLGMDGenes_For_RareLGMDSNPs_GeneAnnotations_File ="commonLGMDGenes_For_RareLGMDSNPs_GeneAnnotations.txt";
		checkAnyCommonGenesBetweenLGMDAssociatedGenesAndGLANETAnalysis(
				geneOvelapFolder,
				overlapAnalysisFile,
				commonLGMDGenes_For_RareLGMDSNPs_GeneAnnotations_File,
				LGMDAssociatedGeneSymbolsList,
				LGMDAssociatedGeneIdsList);
		
		
		String postAnalysisRSAFolder = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\PostAnalysisRSA\\";
		String postAnalysisRSAFile = "PostAnalysisofRegulatorySequenceAnalysisResults_AugmentedWithGeneAnnotations.txt";
		String commonLGMDGenes_For_RareLGMDSNPs_PostAnalysisRSA_File ="commonLGMDGenes_For_RareLGMDSNPs_PostAnalysisRSA.txt";
		
		checkAnyCommonGenesBetweenLGMDAssociatedGenesAndGLANETAnalysis(
				postAnalysisRSAFolder,
				postAnalysisRSAFile,
				commonLGMDGenes_For_RareLGMDSNPs_PostAnalysisRSA_File,
				LGMDAssociatedGeneSymbolsList,
				LGMDAssociatedGeneIdsList);
		
		
		

	}

}
