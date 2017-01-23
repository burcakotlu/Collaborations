/**
 * 
 */
package hacettepe.lgmd.commonTFs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.StringReader;
import java.rmi.RemoteException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.xml.rpc.ServiceException;

import jaxbxjctool.NCBIEutils;
import remap.Remap;
import rsat.GenerationofSequencesandMatricesforSNPs;
import rsat.PositionFrequencyAndLogoMatrices;
import rsat.RegulatorySequenceAnalysisUsingRSATMatrixScan;
import RSATWS.MatrixScanRequest;
import RSATWS.MatrixScanResponse;
import RSATWS.RSATWSPortType;
import RSATWS.RSATWebServicesLocator;
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

	
	public static void fillUsingTheFirstResult( 
			RSATResult rsatResult,
			String matrixScanResult) {

		String resultLine = null;

		BufferedReader bufferedReader = new BufferedReader( new StringReader(matrixScanResult));

		try{

			while( ( resultLine = bufferedReader.readLine()) != null){

				if( resultLine.charAt( 0) != ';' && resultLine.charAt( 0) != '#'){
					fillRSATResult(resultLine,rsatResult);
					break;
				}// End of IF

			}// End of WHILE

			// Close bufferedReader
			bufferedReader.close();

		}catch( IOException e){
			System.out.println(e.toString());
		}

	}

	public static void fillRSATResult(
			String resultLine,
			RSATResult rsatResult) {

		/*********************** SET DECIMAL FORMAT SEPARATORS *****************************/
		//DecimalFormat df = GlanetDecimalFormat.getGLANETDecimalFormat( "0.######E0");
		/*********************** SET DECIMAL FORMAT SEPARATORS *****************************/

		String matrixName = null;
		int indexofDot;
		int matrixNumber = 0;

		char direction = ' ';
		int start = 0;
		int end = 0;
		String sequence = null;
		double pValue = 0;

		// example result line
		// gi|568815587:47291327-47291355 site matrix-scan-matrix_2015-01-27.7 D
		// 17 26 TTCACTGGAC 2.8 1.0e-02 -4.562 1.981 1 1

		int indexofFirstTab = resultLine.indexOf( '\t');
		int indexofSecondTab = ( indexofFirstTab > 0)?resultLine.indexOf( '\t', indexofFirstTab + 1):-1;
		int indexofThirdTab = ( indexofSecondTab > 0)?resultLine.indexOf( '\t', indexofSecondTab + 1):-1;
		int indexofFourthTab = ( indexofThirdTab > 0)?resultLine.indexOf( '\t', indexofThirdTab + 1):-1;
		int indexofFifthTab = ( indexofFourthTab > 0)?resultLine.indexOf( '\t', indexofFourthTab + 1):-1;
		int indexofSixthTab = ( indexofFifthTab > 0)?resultLine.indexOf( '\t', indexofFifthTab + 1):-1;
		int indexofSeventhTab = ( indexofSixthTab > 0)?resultLine.indexOf( '\t', indexofSixthTab + 1):-1;
		int indexofEigthTab = ( indexofSeventhTab > 0)?resultLine.indexOf( '\t', indexofSeventhTab + 1):-1;
		int indexofNinethTab = ( indexofEigthTab > 0)?resultLine.indexOf( '\t', indexofEigthTab + 1):-1;

		// RSAT convention MatrixNumber = LastNumber -1
		// matrix name
		// 1 matrix-scan_2014-08-29.2
		// 2 matrix-scan_2014-08-29.3
		// 3 matrix-scan_2014-08-29.4
		// 4 matrix-scan_2014-08-29.5
		// can be matrix-scan_2014-08-29

		// matrix-scan-matrix_2014-08-29.14 means 13rd matrix in the file
		matrixName = resultLine.substring( indexofSecondTab + 1, indexofThirdTab);
		indexofDot = matrixName.indexOf( '.');
		if( indexofDot >= 0){
			matrixNumber = Integer.parseInt( matrixName.substring( indexofDot + 1)) - 1;
		}else{
			// if there is no '.' this means that there is only one matrix which
			// is numbered 1.
			matrixNumber = 1;
		}

		direction = resultLine.substring( indexofThirdTab + 1, indexofFourthTab).charAt( 0);
		start = Integer.parseInt( resultLine.substring( indexofFourthTab + 1, indexofFifthTab));
		end = Integer.parseInt( resultLine.substring( indexofFifthTab + 1, indexofSixthTab));
		sequence = resultLine.substring( indexofSixthTab + 1, indexofSeventhTab);
		pValue = Double.parseDouble( resultLine.substring( indexofEigthTab + 1, indexofNinethTab));

		rsatResult.setMatrixName(matrixName);
		rsatResult.setMatrixNumber(matrixNumber);
		
		rsatResult.setDirection(direction);
		rsatResult.setStart(start);
		rsatResult.setEnd(end);
		rsatResult.setSequence(sequence);
		rsatResult.setpValue(pValue);

	}

	public static String matrixScan(
			String dnaSequenceInFasta, 
			String pfmMatrices,
			MatrixScanRequest matrixScanRequest, 
			RSATWSPortType proxy) {

		String result = null;

		try{

			matrixScanRequest.setSequence(dnaSequenceInFasta);
			matrixScanRequest.setMatrix(pfmMatrices);

			/***************************************************/
			/************ Old Code starts ************************/
			/***************************************************/
			/* Call the service */
			MatrixScanResponse response;
			response = proxy.matrix_scan(matrixScanRequest);

			/* Get the result */
			result = response.getClient();
			/***************************************************/
			/************ Old Code ends **************************/
			/***************************************************/
	
		}catch( RemoteException e){
			System.out.println(e.toString());
		}

		return result;
	}

	public static void getDNASequenceCallRSATWriteResults(
			Map<String,String> ENCODE_TFName2PFMMap,
			Map<String,Interval> intervalName2IntervalMap){
		
		
		//First way
		Map<String,List<RSATResult>>  intervalName2RSATResultListMap = new HashMap<String,List<RSATResult>>();
		Map<String,List<RSATResult>>  tfName2RSATResultListMap = new HashMap<String,List<RSATResult>>();
		List<RSATResult> rsatResultList = null;
		List<RSATResult> tfBased_RSATResultList = null;
		RSATResult rsatResult = null;

//		//Second way
//						
//		List<TFBasedResult> TFBasedResultList = new ArrayList<TFBasedResult>();
//		TFBasedResult tfBasedResult = null;
//			
//		Double avgPValue = 0d;
//		Double smallestPValue = 0d;
//		Double accumulatedLogRatio = 0d;		
		
		String intervalName = null;
		Interval interval = null;
		
		String tfName = null;
		String tfPFM = null;
				
		ChromosomeName chromosomeName = null;
		String chrName = null;
		String fastaFile = null;
		
		String dataFolder = "C:\\Users\\Burçak\\Google Drive\\Data\\";
		Map<String, String> chrName2RefSeqIdforLatestAssemblyReturnedByNCBIEutilsMap = new HashMap<String, String>();	
		
		FileWriter fileWriter = null;
		BufferedWriter bufferedWriter = null;

//		FileWriter fileWriter2 = null;
//		BufferedWriter bufferedWriter2 = null;

		String RSATMatrixScanResult = null;
			

		try {
			fileWriter = FileOperations.createFileWriter("C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\CommonTFs\\RSATResults.txt");
			bufferedWriter = new BufferedWriter( fileWriter);
	
//			fileWriter2 = FileOperations.createFileWriter("C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\CommonTFs\\RSATResults_sorted_wrt_avg_pvalue.txt");
//			bufferedWriter2 = new BufferedWriter(fileWriter2);
	
		
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
			
			//Remap.remap_show_batches(dataFolder, Commons.NCBI_REMAP_API_SUPPORTED_ASSEMBLIES_FILE);
			
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

			
			/***************************************************************************************/
			/***********************************RSAT starts*****************************************/
			/***************************************************************************************/
			RSATWebServicesLocator service = new RSATWebServicesLocator();

			//User the server address for RSAT Metazoa
			service.setRSATWSPortTypeEndpointAddress("http://rsat.sb-roscoff.fr/web_services/RSATWS.cgi");

			RSATWSPortType proxy = null;

			try{
				proxy = service.getRSATWSPortType();
			}catch( ServiceException e){
			}
			/***************************************************************************************/
			/***********************************RSAT ends*******************************************/
			/***************************************************************************************/
			
			
			/**************************************************************************************************************************/
			/************************* Initialize Matrix Scan Request Parameters starts ***********************************************/
			/**************************************************************************************************************************/
			MatrixScanRequest matrixScanRequest = new MatrixScanRequest();
			
			RegulatorySequenceAnalysisUsingRSATMatrixScan.initializeMatrixScanParameters(matrixScanRequest);
			/**************************************************************************************************************************/
			/************************* Initialize Matrix Scan Request Parameters ends *************************************************/
			/**************************************************************************************************************************/
		
			
//			//The Second Way
//			//Call RSAT
//			//For each TF
//			//And then for each interval
//			for(Map.Entry<String, String> tfEntry : ENCODE_TFName2PFMMap.entrySet()){
//				
//				tfName = tfEntry.getKey();
//				tfPFM = tfEntry.getValue();
//				
//				for(Map.Entry<String, Interval> intervalEntry:intervalName2IntervalMap.entrySet()){
//					intervalName = intervalEntry.getKey();
//					interval = intervalEntry.getValue();
//					
//					chromosomeName = interval.getChromosomeName();
//					chrName = chromosomeName.convertEnumtoString();
//					
//					fastaFile = GenerationofSequencesandMatricesforSNPs.getDNASequence(
//							chrName.substring(3),
//							interval.getLow(),
//							interval.getHigh(),
//							chrName2RefSeqIdforLatestAssemblyReturnedByNCBIEutilsMap);
//					
//					rsatResult = new RSATResult();
//					
//					rsatResult.setIntervalName(intervalName);
//					rsatResult.setTfName(tfName);
//					
//					RSATMatrixScanResult = matrixScan( 
//							fastaFile, 
//							tfPFM,  
//							matrixScanRequest, 
//							proxy);
//					
//					if( RSATMatrixScanResult != null){
//						
//						fillUsingTheFirstResult(rsatResult,RSATMatrixScanResult);
//						
//						rsatResultList = tfName2RSATResultListMap.get(tfName);
//						
//						if (rsatResultList==null){
//							rsatResultList = new ArrayList<RSATResult>();
//							rsatResultList.add(rsatResult);
//							tfName2RSATResultListMap.put(tfName, rsatResultList);
//						}else{
//							rsatResultList.add(rsatResult);
//						}
//					}//End of IF RSATResult is not null
//					
//					
//				}//End of for each interval
//				
//			}//End of for each Encode TF
				
			
			//First way
			for(Map.Entry<String, Interval> intervalEntry:intervalName2IntervalMap.entrySet()){
				
				intervalName = intervalEntry.getKey();
				interval = intervalEntry.getValue();
				
				chromosomeName = interval.getChromosomeName();
				chrName = chromosomeName.convertEnumtoString();				
				
				//TODO
				//Think about it.
				//Can we declare strand for genes on negative strand? Yes for positive strand, we give strand=1
				//Most probably for negative strand, we give strand=2, check it.
				//get DNA Sequence
				//Check it: whether gene on - strand get the best matches with TFs on Reverse side
				//RSAT looks for both strands so getting DNA sequence on + or - strand does not matter.
				fastaFile = GenerationofSequencesandMatricesforSNPs.getDNASequence(
						chrName.substring(3),
						interval.getLow(),
						interval.getHigh(),
						chrName2RefSeqIdforLatestAssemblyReturnedByNCBIEutilsMap);
	
				
				for(Map.Entry<String, String> tfEntry : ENCODE_TFName2PFMMap.entrySet()){
					
					tfName = tfEntry.getKey();
					tfPFM = tfEntry.getValue();
					
					rsatResult = new RSATResult();
					
					rsatResult.setIntervalName(intervalName);
					rsatResult.setTfName(tfName);
					
					RSATMatrixScanResult = matrixScan( 
							fastaFile, 
							tfPFM,  
							matrixScanRequest, 
							proxy);
					
					if( RSATMatrixScanResult != null){
						
						fillUsingTheFirstResult(rsatResult,RSATMatrixScanResult);
						
						rsatResultList = intervalName2RSATResultListMap.get(intervalName);
						
						if (rsatResultList==null){
							rsatResultList = new ArrayList<RSATResult>();
							rsatResultList.add(rsatResult);
							intervalName2RSATResultListMap.put(intervalName, rsatResultList);
						}else{
							rsatResultList.add(rsatResult);
						}
						
						//Fill other map
						tfBased_RSATResultList = tfName2RSATResultListMap.get(tfName);
						if (tfBased_RSATResultList==null){
							tfBased_RSATResultList = new ArrayList<RSATResult>();
							tfBased_RSATResultList.add(rsatResult);
							tfName2RSATResultListMap.put(tfName, tfBased_RSATResultList);
						}else{
							tfBased_RSATResultList.add(rsatResult);
						}
						
					}
						
				}//End of FOR each TF PFM
								
			}//End of FOR each interval
			
			/***********************************************************************/
			/***********************Output for first way starts*********************/
			/***********************************************************************/
			//Now sort the results with ascending p-Value
			for (Map.Entry<String, List<RSATResult>> entry: intervalName2RSATResultListMap.entrySet()){
				rsatResultList = entry.getValue();
				Collections.sort(rsatResultList, RSATResult.P_VALUE);
			}
			
			//Write header line
			bufferedWriter.write(
					"IntervalName"+ "\t" + 
					"TfName"+ "\t" + 					
					"MatrixName" + "\t" + 
					"MatrixNumber"+ "\t" + 
					"Direction"+ "\t" + 
					"Start"+ "\t" + 
					"End"+ "\t" + 
					"pValue"+ "\t" +
					System.getProperty("line.separator"));
			
			
			//Now write the results with ascending p-Value
			for (Map.Entry<String, List<RSATResult>> entry: intervalName2RSATResultListMap.entrySet()){
				
				intervalName = entry.getKey();
				rsatResultList = entry.getValue();
				
				for(int i=0;i<rsatResultList.size();i++){		
										
					bufferedWriter.write(
							rsatResultList.get(i).getIntervalName()+ "\t" + 
							rsatResultList.get(i).getTfName()+ "\t" + 
							rsatResultList.get(i).getMatrixName()+ "\t" + 
							rsatResultList.get(i).getMatrixNumber()+ "\t" + 
							rsatResultList.get(i).getDirection()+ "\t" + 
							rsatResultList.get(i).getStart()+ "\t" + 
							rsatResultList.get(i).getEnd()+ "\t" + 
							rsatResultList.get(i).getpValue()+ "\t" +
							System.getProperty("line.separator"));
										
				}//End of for each RSATResult												
				
			}//End of for each RSATResultlist			
			/***********************************************************************/
			/***********************Output for first way ends***********************/
			/***********************************************************************/
			
			findTheCommonTFs(intervalName2RSATResultListMap,tfName2RSATResultListMap);

			
//			/***********************************************************************/
//			/***********************Output for second way starts********************/
//			/***********************************************************************/
//			//Now sort the results with ascending p-Value
//			for (Map.Entry<String, List<RSATResult>> entry: tfName2RSATResultListMap.entrySet()){
//				rsatResultList = entry.getValue();
//				Collections.sort(rsatResultList, RSATResult.P_VALUE);
//			}
//			
//			
//			
//			//Write header line
//			bufferedWriter.write(
//					"TfName"+ "\t" + 
//					"IntervalName"+ "\t" + 
//					"MatrixName" + "\t" + 
//					"MatrixNumber"+ "\t" + 
//					"Direction"+ "\t" + 
//					"Start"+ "\t" + 
//					"End"+ "\t" + 
//					"pValue"+ "\t" +
//					System.getProperty("line.separator"));
//			
//			//Now write the results with ascending p-Value
//			for (Map.Entry<String, List<RSATResult>> entry: tfName2RSATResultListMap.entrySet()){
//				
//				tfName = entry.getKey();
//				rsatResultList = entry.getValue();
//				
//				avgPValue = 0d;
//				accumulatedLogRatio = 0d;
//				
//				for(int i=0;i<rsatResultList.size();i++){		
//					if (i==0){
//						smallestPValue = rsatResultList.get(i).getpValue();
//					}else{
//						accumulatedLogRatio += Math.log(smallestPValue/rsatResultList.get(i).getpValue());
//					}
//					
//					bufferedWriter.write(
//							rsatResultList.get(i).getTfName()+ "\t" + 
//							rsatResultList.get(i).getIntervalName()+ "\t" + 
//							rsatResultList.get(i).getMatrixName()+ "\t" + 
//							rsatResultList.get(i).getMatrixNumber()+ "\t" + 
//							rsatResultList.get(i).getDirection()+ "\t" + 
//							rsatResultList.get(i).getStart()+ "\t" + 
//							rsatResultList.get(i).getEnd()+ "\t" + 
//							rsatResultList.get(i).getpValue()+ "\t" +
//							System.getProperty("line.separator"));
//					
//					avgPValue += rsatResultList.get(i).getpValue();
//					
//				}//End of for each RSATResult
//				
//				avgPValue = avgPValue/rsatResultList.size();
//				
//				tfBasedResult = new TFBasedResult(tfName,avgPValue, accumulatedLogRatio);
//				TFBasedResultList.add(tfBasedResult);
//				
//			}//End of for each RSATResultlist
//			
//			//Sort TFBasedResultList w.r.t. avgPValue in ascending order
//			Collections.sort(TFBasedResultList, TFBasedResult.AVG_P_VALUE);
//			
//			
//			//Write header
//			bufferedWriter2.write("TfName" + "\t" +
//					"AvgPValue" + "\t" +
//					"AccumulatedLogRatio" + "\t" +
//					System.getProperty("line.separator"));
//			for(int i=0;i<TFBasedResultList.size();i++){
//				tfBasedResult = TFBasedResultList.get(i);
//				bufferedWriter2.write(tfBasedResult.getTfName() + "\t" +
//						tfBasedResult.getAvgPValue() + "\t" +
//						tfBasedResult.getAccumulatedLogRatio() + "\t" +
//						System.getProperty("line.separator"));
//			}
//			/***********************************************************************/
//			/***********************Output for second way ends**********************/
//			/***********************************************************************/

			
			//Close bufferedWriter
			bufferedWriter.close();
//			bufferedWriter2.close();
		
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}

		
	}

	
	
	public  static void findTheCommonTFs(
			Map<String,List<RSATResult>> intervalName2RSATResultListMap,
			Map<String,List<RSATResult>> tfName2RSATResultListMap){
		
		FileWriter fileWriter = null;
		BufferedWriter bufferedWriter = null;
		
		String intervalName = null;
		List<RSATResult> RSATResultList = null;
		
		RSATResult rsatResult = null;
		double threshold = 0.0001d;
		
		Set<String> intersection = null;
				
		try {
			fileWriter = FileOperations.createFileWriter("C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\CommonTFs\\CommonTFs.txt");
			bufferedWriter = new BufferedWriter( fileWriter);

			Map<String,Set<String>> intervalName2TFSetMap = new HashMap<String,Set<String>>();
			Set<String> tfSet = null;
			
			//Initialize intervalName2TFSetMap
			for(Map.Entry<String,List<RSATResult>> entry:intervalName2RSATResultListMap.entrySet()){
				intervalName = entry.getKey();
				intervalName2TFSetMap.put(intervalName, new HashSet<String>());
			}
			
			for(Map.Entry<String,List<RSATResult>> entry:intervalName2RSATResultListMap.entrySet()){
				
				intervalName = entry.getKey();
				RSATResultList = entry.getValue();
				
				for(int i=0; i<RSATResultList.size();i++){
					
					rsatResult = RSATResultList.get(i);
					
					//Fill intervalName2TFSetMap				
					if (rsatResult.getpValue()<threshold){
						
						tfSet = intervalName2TFSetMap.get(intervalName);
						
						if (tfSet==null){
							
							tfSet = new HashSet<String>();
							tfSet.add(rsatResult.getTfName());
							
							intervalName2TFSetMap.put(intervalName, tfSet);
							
						}else{
							tfSet.add(rsatResult.getTfName());						
						}
						
						
					}//End of IF
					
				}//End of for TF RSATResult
				
			}//End of for each intervalName
			
			
			//Initialize intersection with any tfSet
			for(Map.Entry<String, Set<String>>  entry: intervalName2TFSetMap.entrySet()) {
				tfSet = entry.getValue();
				intersection = new HashSet<String>(tfSet);
				break;			
			}//End of for
			
			//Find common TFs in sets
			for(Map.Entry<String, Set<String>>  entry: intervalName2TFSetMap.entrySet()) {
				tfSet = entry.getValue();
				intersection.retainAll(tfSet);
			
			}//End of for
			
			for(String commonTF: intersection) {
				bufferedWriter.write(commonTF + System.getProperty("line.separator") );	
				RSATResultList =  tfName2RSATResultListMap.get(commonTF);
				for(int i=0; i<RSATResultList.size();i++){
					bufferedWriter.write(
							RSATResultList.get(i).getIntervalName()+ "\t" + 
							RSATResultList.get(i).getTfName()+ "\t" + 
							RSATResultList.get(i).getMatrixName()+ "\t" + 
							RSATResultList.get(i).getMatrixNumber()+ "\t" + 
							RSATResultList.get(i).getDirection()+ "\t" + 
							RSATResultList.get(i).getStart()+ "\t" + 
							RSATResultList.get(i).getEnd()+ "\t" + 
							RSATResultList.get(i).getpValue()+ "\t" +
							System.getProperty("line.separator"));
				}
				
			}
			
			//Close
			bufferedWriter.close();
		
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		
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
		
		/***************************************************/
		/*********************Task1 starts******************/
		/***************************************************/
		//Create Intervals for given genes		
		List<String> lgmdGeneSymbolList = new  ArrayList<String>();
		lgmdGeneSymbolList.add("SGCA");
		lgmdGeneSymbolList.add("SGCB");
		lgmdGeneSymbolList.add("SGCG");
		lgmdGeneSymbolList.add("SGCD");
		Map<String,Interval> intervalName2IntervalMap = new HashMap<String,Interval>();
		createIntervals(lgmdGeneSymbolList,intervalName2IntervalMap);
		/***************************************************/
		/*********************Task1 ends********************/
		/***************************************************/
		

		/***************************************************/
		/*********************Task2 starts******************/
		/***************************************************/
		//Consider only ENCODE TFs
		//Bookkeeping
		//https://www.encodeproject.org/matrix/?type=Experiment
		List<String> ENCODE_TFs_List = new ArrayList<String>();
		getENCODETFs(ENCODE_TFs_List);
		
//		//for testing purposes
//		List<String> small_ENCODE_TFs_List = new ArrayList<String>();
//		for(int i=0; i<5;i++){
//			small_ENCODE_TFs_List.add(ENCODE_TFs_List.get(i));
//		}
		/***************************************************/
		/*********************Task2 ends********************/
		/***************************************************/

		
		/***************************************************/
		/*********************Task3 starts******************/
		/***************************************************/
		// Construct position frequency matrices from Jaspar Core
		// Construct logo matrices from Jaspar Core
		Map<String,String> tfName2TFPFMMap =new HashMap<String,String>();
		Map<String,String> tfName2TFLogoMatricesMap =new HashMap<String,String>();		
			
		String dataFolder = "C:\\Users\\Burçak\\Google Drive\\Data\\";
		String jasparCoreInputFileName = Commons.JASPAR_CORE;
		PositionFrequencyAndLogoMatrices.constructPfmMatricesandLogoMatricesfromJasparCore(
				dataFolder, 
				jasparCoreInputFileName, 
				tfName2TFPFMMap,
				tfName2TFLogoMatricesMap);
		
		Map<String,String> ENCODE_TFName2PFMMap = new HashMap<String,String>();
		fillENCODETFName2PFMMap(ENCODE_TFs_List,tfName2TFPFMMap,ENCODE_TFName2PFMMap);
		/***************************************************/
		/*********************Task3 ends********************/
		/***************************************************/

		
	
		/***************************************************/
		/*********************Task4 starts******************/
		/***************************************************/
		//Get DNA sequences for these intervals for the latest assembly
		//Call RSAT web service for each interval-TF pair
		//Sort RSAT result for each gene in ascending order (the lower the p-value, the better the macth is)
		//Find the common best matching TFs in top ten for each gene
		//Write the RSAT outputs
		getDNASequenceCallRSATWriteResults(ENCODE_TFName2PFMMap,intervalName2IntervalMap);
		/***************************************************/
		/*********************Task4 ends********************/
		/***************************************************/

		
	}

}
