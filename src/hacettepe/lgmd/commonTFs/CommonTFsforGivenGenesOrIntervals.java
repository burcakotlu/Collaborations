/**
 * 
 */
package hacettepe.lgmd.commonTFs;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jaxbxjctool.NCBIEutils;
import remap.Remap;
import rsat.GenerationofSequencesandMatricesforSNPs;
import rsat.PositionFrequencyAndLogoMatrices;
import auxiliary.FileOperations;
import auxiliary.Interval;
import common.Commons;
import enumtypes.ChromosomeName;

/**
 * @author Burçak Otlu
 * @date Jan 2, 2017
 * @project Collaborations 
 *
 */
public class CommonTFsforGivenGenesOrIntervals {

	

	public static void getDNASequenceCallRSATWriteResults(
			Map<String,String> ENCODE_TFName2PFMMap,
			Map<String,Interval> intervalName2IntervalMap){
		
		String intervalName = null;
		Interval interval = null;
		
		ChromosomeName chromosomeName = null;
		String chrName = null;
		String fastaFile = null;
		String referenceSequence = null;
		
		String dataFolder = "C:\\Users\\Burçak\\Google Drive\\Data\\";
		Map<String, String> chrName2RefSeqIdforLatestAssemblyReturnedByNCBIEutilsMap = new HashMap<String, String>();	
		
		
		
		/***************************************************************************************/
		/***************************************************************************************/
		/***************************************************************************************/
		String latestAssemblyNameReturnedByNCBIEutils = NCBIEutils.getLatestAssemblyNameReturnedByNCBIEutils();
		/***************************************************************************************/
		/***************************************************************************************/
		/***************************************************************************************/

		
		/***************************************************************************************/
		/***************************************************************************************/
		/***************************************************************************************/
		Map<String, String> assemblyName2RefSeqAssemblyIDMap = new HashMap<String, String>();
		
		Remap.remap_show_batches(dataFolder, Commons.NCBI_REMAP_API_SUPPORTED_ASSEMBLIES_FILE);
		
		Remap.fillAssemblyName2RefSeqAssemblyIDMap(
				dataFolder, 
				Commons.NCBI_REMAP_API_SUPPORTED_ASSEMBLIES_FILE,
				assemblyName2RefSeqAssemblyIDMap);
		/***************************************************************************************/
		/***************************************************************************************/
		/***************************************************************************************/


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

		
		for(Map.Entry<String, Interval> entry:intervalName2IntervalMap.entrySet()){
			intervalName = entry.getKey();
			interval = entry.getValue();
			
			chromosomeName = interval.getChromosomeName();
			chrName = chromosomeName.convertEnumtoString();
			
			
			//TODO
			//Think about it.
			//Can we declare strand for genes on negative strand? Yes for positive strand, we give strand=1
			//Most probably for negative strand, we give strand=2, check it.
			//get DNA Sequence
			//Check it: whether gene on - strand get the best matches on Reverse
			fastaFile = GenerationofSequencesandMatricesforSNPs.getDNASequence(
					chrName.substring(3),
					interval.getLow(),
					interval.getHigh(),
					chrName2RefSeqIdforLatestAssemblyReturnedByNCBIEutilsMap);

			referenceSequence = GenerationofSequencesandMatricesforSNPs.getDNASequenceFromFastaFile(fastaFile);
			
			//TODO
			//Call RSAT for each ENCODE TF PFM and reference sequence
			//look at matrixScan of RegulatorySequenceAnalysisUsingRSATMatrixScan

			
		}//End of FOR each interval
		
	}

	
	public static void fillENCODETFName2PFMMap(
			List<String> ENCODE_TFs_List,
			Map<String, String> tfName2TFPFMMap,
			Map<String,String> ENCODE_TFName2PFMMap){
		
		
		for(String tfName : ENCODE_TFs_List){
			
			for(Map.Entry<String, String> entry:tfName2TFPFMMap.entrySet()){
				if (entry.getKey().contains(tfName)){
					
					if (ENCODE_TFName2PFMMap.get(tfName)==null){
						ENCODE_TFName2PFMMap.put(tfName, entry.getValue());
						
					}else{
						//append 
						ENCODE_TFName2PFMMap.put(tfName,ENCODE_TFName2PFMMap.get(tfName) + entry.getValue());
					}
					
					System.out.println("EntryKey\t" + entry.getKey() +  "\ttfName\t" +  tfName);
				}
			}//End of for each TF in Jaspar Core			
			
		}//End of FOR each ENCODE TF
		
		System.out.println("xxx " + ENCODE_TFName2PFMMap.size());
		
	}
	
	
	public static void getENCODETFs(List<String> ENCODE_TFs_List){
		
		
		// ID	Accession	Assay Type	Assay Nickname	Target label	Target gene	Biosample summary	Biosample	Description	Lab	Project	Status	Biological replicate	Technical replicate	Linked Antibody	Species	Life stage	Age	Age Units	Treatment	Term ID	Concentration	Concentration units	Duration	Duration units	Synchronization	Post-synchronization time	Post-synchronization time units	Replicates
		// /experiments/ENCSR000ECT/	ENCSR000ECT	ChIP-seq	ChIP-seq	POLR2AphosphoS2	POLR2A	HeLa-S3	HeLa-S3	POLR2AphosphoS2 ChIP-seq on human HeLa-S3	Michael Snyder, Stanford	ENCODE	released	2,1	1	ENCAB000AOB	Homo sapiens	adult	31	year										/replicates/5254b72d-f123-4888-b268-2a4c89a6145d/,/replicates/0832cd43-ca61-4e28-8473-f902df66b109/
		// /experiments/ENCSR000DMY/	ENCSR000DMY	ChIP-seq	ChIP-seq	CTCF	CTCF	medulloblastoma	medulloblastoma	CTCF ChIP-seq on human Medullo	Vishwanath Iyer, UTA	ENCODE	released	2,1	1	ENCAB000AXY	Homo sapiens	child	2	year										/replicates/14bf8d50-3785-4e99-a05d-370930ea5dd1/,/replicates/5f25328f-abfb-41ec-b927-ba14bdfe27a1/

		
		FileReader fileReader = null;
		BufferedReader bufferedReader = null;
		
		String strLine = null;
		
		int indexofFirstTab = 0;
		int indexofSecondTab = 0;
		int indexofThirdTab = 0;
		int indexofFourthTab = 0;
		int indexofFifthTab = 0;
		
		int index = -1;
		int index_underscore = -1;
		
		String ENCODE_TFs_FileName = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\ENCODE_TFs\\ENCODE_Human_ChIP_Seq_report.tsv";		
		String tfName = null;
		
		try {

			fileReader = FileOperations.createFileReader(ENCODE_TFs_FileName);
			bufferedReader = new BufferedReader(fileReader);
			
			strLine = bufferedReader.readLine();
			
			while((strLine = bufferedReader.readLine())!=null){
				
				indexofFirstTab = strLine.indexOf( '\t');
				indexofSecondTab = ( indexofFirstTab > 0)?strLine.indexOf( '\t', indexofFirstTab + 1):-1;
				indexofThirdTab = ( indexofSecondTab > 0)?strLine.indexOf( '\t', indexofSecondTab + 1):-1;
				indexofFourthTab = ( indexofThirdTab > 0)?strLine.indexOf( '\t', indexofThirdTab + 1):-1;
				indexofFifthTab = ( indexofFourthTab > 0)?strLine.indexOf( '\t', indexofFourthTab + 1):-1;
				
				tfName = strLine.substring(indexofFourthTab+1, indexofFifthTab);
				
				index = tfName.indexOf('-');
				if (index!=-1){
					tfName = tfName.substring(index+1);
				}
				
				
				index_underscore = tfName.indexOf('_');
				if (index_underscore!=-1){
					tfName = tfName.substring(0,index_underscore);
				}
						
				
				if (!ENCODE_TFs_List.contains(tfName)){
					ENCODE_TFs_List.add(tfName);
				}
				
			}//End of WHILE
			
			
			bufferedReader.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}

	public static void createIntervals(
			List<String> lgmdGeneSymbolList,
			Map<String,Interval> geneSymbol2IntervalMap){
		
		String dataFolder = "C:\\Users\\Burçak\\Google Drive\\Data\\";
		
		// Downloaded from UCSC Genome Table Browser contains RNA_NUCLEOTIDE_ACCESSION and GENE_SYMBOL
		String UCSC_GENOME_HG38_REFSEQ_GENES_FILE = dataFolder + Commons.UCSCGENOME_HG38_REFSEQ_GENES_DOWNLOADED_2_DEC_2016;
		
		FileReader fileReader = null;
		BufferedReader bufferedReader = null;
		
		String strLine = null;
		
		int indexofFirstTab = 0;
		int indexofSecondTab = 0;
		int indexofThirdTab = 0;
		int indexofFourthTab = 0;
		int indexofFifthTab = 0;
		int indexofSixthTab = 0;
		int indexofSeventhTab = 0;
		int indexofEigthTab = 0;
		int indexofNinethTab = 0;
		int indexofTenthTab = 0;
		int indexofEleventhTab = 0;
		int indexofTwelfthTab = 0;
		int indexofThirteenthTab = 0;

		//String refSeqGeneName;
		ChromosomeName chromName;
		char strand;
		int txStart;
		int txEnd;		
		String alternateGeneName;
		
		
		int start = -1;
		int end = -1;
		
		try {
			
			fileReader = FileOperations.createFileReader(UCSC_GENOME_HG38_REFSEQ_GENES_FILE);
			bufferedReader = new BufferedReader(fileReader);
			
			//Skip header line
			strLine = bufferedReader.readLine();
			//#bin	name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	score	name2	cdsStartStat	cdsEndStat	exonFrames

			
			while( ( strLine = bufferedReader.readLine()) != null){
				
				indexofFirstTab = strLine.indexOf( '\t');
				indexofSecondTab = ( indexofFirstTab > 0)?strLine.indexOf( '\t', indexofFirstTab + 1):-1;
				indexofThirdTab = ( indexofSecondTab > 0)?strLine.indexOf( '\t', indexofSecondTab + 1):-1;
				indexofFourthTab = ( indexofThirdTab > 0)?strLine.indexOf( '\t', indexofThirdTab + 1):-1;
				indexofFifthTab = ( indexofFourthTab > 0)?strLine.indexOf( '\t', indexofFourthTab + 1):-1;
				indexofSixthTab = ( indexofFifthTab > 0)?strLine.indexOf( '\t', indexofFifthTab + 1):-1;
				indexofSeventhTab = ( indexofSixthTab > 0)?strLine.indexOf( '\t', indexofSixthTab + 1):-1;
				indexofEigthTab = ( indexofSeventhTab > 0)?strLine.indexOf( '\t', indexofSeventhTab + 1):-1;
				indexofNinethTab = ( indexofEigthTab > 0)?strLine.indexOf( '\t', indexofEigthTab + 1):-1;
				indexofTenthTab = ( indexofNinethTab > 0)?strLine.indexOf( '\t', indexofNinethTab + 1):-1;
				indexofEleventhTab = ( indexofTenthTab > 0)?strLine.indexOf( '\t', indexofTenthTab + 1):-1;
				indexofTwelfthTab = ( indexofEleventhTab > 0)?strLine.indexOf( '\t', indexofEleventhTab + 1):-1;
				indexofThirteenthTab = ( indexofTwelfthTab > 0)?strLine.indexOf( '\t', indexofTwelfthTab + 1):-1;

				//refSeqGeneName = strLine.substring( indexofFirstTab + 1, indexofSecondTab);
				chromName = ChromosomeName.convertStringtoEnum( strLine.substring( indexofSecondTab + 1,
						indexofThirdTab));

				strand = strLine.substring( indexofThirdTab + 1, indexofFourthTab).trim().charAt( 0);

				txStart = Integer.parseInt( strLine.substring( indexofFourthTab + 1, indexofFifthTab));
				// Convert one based end to zero based end
				txEnd = Integer.parseInt( strLine.substring( indexofFifthTab + 1, indexofSixthTab)) - 1;

				alternateGeneName = strLine.substring( indexofTwelfthTab + 1, indexofThirteenthTab);
				
				
				//Left here
				for(String geneSymbol: lgmdGeneSymbolList){
					
					if (alternateGeneName.equalsIgnoreCase(geneSymbol)){
						
						System.out.println(strLine);
						System.out.println(alternateGeneName + " gene length: " + (txEnd-txStart));
						
						switch(strand){
						
							case '+': 
								start = txStart-2000;
								end = txStart +500;
								geneSymbol2IntervalMap.put(geneSymbol + "_promoter" , new Interval(chromName,start,end));
								break;
								
							case '-':
								start = txEnd -500;
								end = txEnd+2000;
								geneSymbol2IntervalMap.put(geneSymbol + "_promoter" , new Interval(chromName,start,end));
								break;
						
						}//End of SWITCH
						
					}//End of IF
					
				}//End of FOR								

			}//End of WHILE reading UCSC_GENOME_HG38_REFSEQ_GENES_FILE
			
			bufferedReader.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
	}

	public static void main(String[] args) {
		
		//Task1
		//Create Intervals for given genes		
		List<String> lgmdGeneSymbolList = new  ArrayList<String>();
		lgmdGeneSymbolList.add("SGCA");
		lgmdGeneSymbolList.add("SGCB");
		lgmdGeneSymbolList.add("SGCG");
		lgmdGeneSymbolList.add("SGCD");
		Map<String,Interval> intervalName2IntervalMap = new HashMap<String,Interval>();
		createIntervals(lgmdGeneSymbolList,intervalName2IntervalMap);
		

		//Task2
		//Consider ENCODE TFs
		//https://www.encodeproject.org/matrix/?type=Experiment
		List<String> ENCODE_TFs_List = new ArrayList<String>();
		getENCODETFs(ENCODE_TFs_List);
		
		//Task3
		// Construct position frequency matrices from Jaspar Core
		// Construct logo matrices from Jaspar Core
		Map<String,String> tfName2TFPFMMap =new HashMap<String,String>();
		Map<String,String> tfName2TFLogoMatricesMap =new HashMap<String,String>();		
		Map<String,String> ENCODE_TFName2PFMMap = new HashMap<String,String>();
			
		String dataFolder = "C:\\Users\\Burçak\\Google Drive\\Data\\";
		String jasparCoreInputFileName = Commons.JASPAR_CORE;
		PositionFrequencyAndLogoMatrices.constructPfmMatricesandLogoMatricesfromJasparCore(
				dataFolder, 
				jasparCoreInputFileName, 
				tfName2TFPFMMap,
				tfName2TFLogoMatricesMap);
		
		fillENCODETFName2PFMMap(ENCODE_TFs_List,tfName2TFPFMMap,ENCODE_TFName2PFMMap);
		
		
		
		//Task4
		//Get DNA sequences for these intervals for the latest assembly
		//Call RSAT web service
		//Write the best TF matches for each gene 
		getDNASequenceCallRSATWriteResults(ENCODE_TFName2PFMMap,intervalName2IntervalMap);
		
		
	}

}
