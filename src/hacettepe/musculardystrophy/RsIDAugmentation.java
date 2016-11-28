/**
 * 
 */
package hacettepe.musculardystrophy;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;

import jaxbxjctool.AugmentationofGivenIntervalwithRsIds;
import jaxbxjctool.NCBIEutils;
import remap.Remap;
import auxiliary.FileOperations;

/**
 * @author Burçak Otlu
 * @date Nov 25, 2016
 * @project Collaboration 
 * 
 * This class prepares input file for GLANET
 *
 */
public class RsIDAugmentation {

	/**
	 * 
	 */
	public RsIDAugmentation() {
		// TODO Auto-generated constructor stub
	}

	
	public static void readRareVariantsAugmentWithRsIDs(
			String dataFolder,
			String RareVariants_ChrName_GLANET1BasedCoordinateLatestAssembly_rsID_InputFileGLANET,
			String RareVariants_ChrName_GLANET1BasedCoordinateLatestAssembly_rsID_Augmented_InputFileGLANET){
		
		//Read 
		FileReader fileReader = null;
		BufferedReader bufferedReader = null;
		
		//FirstTask
		FileWriter fileWriter= null;
		BufferedWriter bufferedWriter = null;
		
		String strLine = null;
		
		String chrName = null;
		
		int _1BasedPosition_LatestAssembly = -1;
		int _1BasedPosition_Plus1_LatestAssembly = -1;
		int _1BasedPosition_Minus1_LatestAssembly = -1;
			
		String rsID = null;
		String rsID_Plus1 = null;
		String rsID_Minus1 = null;

		int indexofFirstTab;
		int indexofSecondTab;
		
		List<Integer> rsIdList = null;
		List<Integer> rsIdList_plus1 = null;
		List<Integer> rsIdList_minus1 = null;
		
		String chrNamewithoutPreceedingChr = null;
		
		int rsIDInt = -1;
		
		int numberofLinesWithNoRsId = 0;
		
		AugmentationofGivenIntervalwithRsIds augmentationOfAGivenIntervalWithRsIDs = null;
		
		try {
			augmentationOfAGivenIntervalWithRsIDs = new AugmentationofGivenIntervalwithRsIds();
		} catch (Exception e1) {
			e1.printStackTrace();
		}

		try {
			fileReader = FileOperations.createFileReader(dataFolder + RareVariants_ChrName_GLANET1BasedCoordinateLatestAssembly_rsID_InputFileGLANET);
			bufferedReader = new BufferedReader(fileReader);
			
			fileWriter = FileOperations.createFileWriter(dataFolder + RareVariants_ChrName_GLANET1BasedCoordinateLatestAssembly_rsID_Augmented_InputFileGLANET);
			bufferedWriter = new BufferedWriter(fileWriter);
			
			while((strLine = bufferedReader.readLine())!=null){
			
				indexofFirstTab = strLine.indexOf('\t');
				indexofSecondTab  = strLine.indexOf('\t',indexofFirstTab+1);
				
				chrName = strLine.substring(0, indexofFirstTab);
				_1BasedPosition_LatestAssembly = Integer.parseInt(strLine.substring(indexofFirstTab+1, indexofSecondTab));
				rsID =  strLine.substring(indexofSecondTab+1);
								
				_1BasedPosition_Plus1_LatestAssembly = _1BasedPosition_LatestAssembly+1; 
				_1BasedPosition_Minus1_LatestAssembly = _1BasedPosition_LatestAssembly-1; 
				
				chrNamewithoutPreceedingChr = chrName.substring(3);
				
				//Initialize
				rsIdList = null;
				rsIdList_plus1 = null;
				rsIdList_minus1 = null;
				
				//Version2
				rsIdList = augmentationOfAGivenIntervalWithRsIDs.getRsIdsInAGivenInterval(
					chrNamewithoutPreceedingChr, 
					_1BasedPosition_LatestAssembly, 
					_1BasedPosition_LatestAssembly);					
									
				if (!rsIdList.isEmpty()){
					rsID = "";
					for(Iterator<Integer> itr=rsIdList.iterator();itr.hasNext();) {
						rsIDInt = itr.next();							
						rsID = rsID + "rs" + rsIDInt + "\t";
					}
					
					bufferedWriter.write(chrName + "\t" + _1BasedPosition_LatestAssembly + "\t" + rsID + "\t" + "Augmented Original Position" +System.getProperty("line.separator"));
				}
			
				else if (rsIdList.isEmpty()){
					
					rsIdList_plus1 = augmentationOfAGivenIntervalWithRsIDs.getRsIdsInAGivenInterval(
							chrNamewithoutPreceedingChr, 
							_1BasedPosition_LatestAssembly+1, 
							_1BasedPosition_LatestAssembly+1);
					
					rsIdList_minus1 = augmentationOfAGivenIntervalWithRsIDs.getRsIdsInAGivenInterval(
							chrNamewithoutPreceedingChr, 
							_1BasedPosition_LatestAssembly-1, 
							_1BasedPosition_LatestAssembly-1);
					
					if (!rsIdList_plus1.isEmpty()){	
						rsID_Plus1 ="";
						for(Iterator<Integer> itr=rsIdList_plus1.iterator();itr.hasNext();) {
							rsIDInt = itr.next();							
							rsID_Plus1 = rsID_Plus1 +  "rs" + rsIDInt + "\t";
						}
						
						bufferedWriter.write(chrName + "\t" + _1BasedPosition_Plus1_LatestAssembly + "\t" + rsID_Plus1 + "\t" + "Augmented Plus1" + System.getProperty("line.separator"));
					}
					
					if (!rsIdList_minus1.isEmpty()){	
						rsID_Minus1="";
						for(Iterator<Integer> itr=rsIdList_minus1.iterator();itr.hasNext();) {
							rsIDInt = itr.next();							
							rsID_Minus1 = rsID_Minus1 + "rs" + rsIDInt + "\t";
						}
						
						bufferedWriter.write(chrName + "\t" + _1BasedPosition_Minus1_LatestAssembly + "\t" + rsID_Minus1 + "\t" + "Augmented Minus1" +System.getProperty("line.separator"));
					}
		
				}//End of if rsIdList is empty

				if (rsIdList.isEmpty() && rsIdList_plus1.isEmpty() && rsIdList_minus1.isEmpty()) {
					numberofLinesWithNoRsId++;
				}
								
//				//Version1
//				if(rsID.equalsIgnoreCase("-")){
//					
//					rsIdList = augmentationOfAGivenIntervalWithRsIDs.getRsIdsInAGivenInterval(
//							chrNamewithoutPreceedingChr, 
//							_1BasedPosition_LatestAssembly, 
//							_1BasedPosition_LatestAssembly);					
//											
//					if (!rsIdList.isEmpty()){						
//							for(Iterator<Integer> itr=rsIdList.iterator();itr.hasNext();) {
//							rsIDInt = itr.next();							
//							rsID = "rs" + rsIDInt + "\t";
//						}
//						
//						bufferedWriter.write(chrName + "\t" + _1BasedPosition_LatestAssembly + "\t" + rsID + "\t" + "Augmented" +System.getProperty("line.separator"));
//					}
//						
//	
//					
//					else if (rsIdList.isEmpty()){
//						
//						rsIdList_plus1 = augmentationOfAGivenIntervalWithRsIDs.getRsIdsInAGivenInterval(
//								chrNamewithoutPreceedingChr, 
//								_1BasedPosition_LatestAssembly+1, 
//								_1BasedPosition_LatestAssembly+1);
//						
//						rsIdList_minus1 = augmentationOfAGivenIntervalWithRsIDs.getRsIdsInAGivenInterval(
//								chrNamewithoutPreceedingChr, 
//								_1BasedPosition_LatestAssembly-1, 
//								_1BasedPosition_LatestAssembly-1);
//						
//						if (!rsIdList_plus1.isEmpty()){						
//							for(Iterator<Integer> itr=rsIdList_plus1.iterator();itr.hasNext();) {
//								rsIDInt = itr.next();							
//								rsID_Plus1 = "rs" + rsIDInt + "\t";
//							}
//							
//							bufferedWriter.write(chrName + "\t" + _1BasedPosition_Plus1_LatestAssembly + "\t" + rsID_Plus1 + "\t" + "Augmented Plus1" + System.getProperty("line.separator"));
//						}
//						
//						if (!rsIdList_minus1.isEmpty()){						
//							for(Iterator<Integer> itr=rsIdList_minus1.iterator();itr.hasNext();) {
//								rsIDInt = itr.next();							
//								rsID_Minus1 = "rs" + rsIDInt + "\t";
//							}
//							
//							bufferedWriter.write(chrName + "\t" + _1BasedPosition_Minus1_LatestAssembly + "\t" + rsID_Minus1 + "\t" + "Augmented Minus1" +System.getProperty("line.separator"));
//
//						}
//		
//					}//End of if rsIdList is empty
//					
//				}else {
//
//					bufferedWriter.write(chrName + "\t" + _1BasedPosition_LatestAssembly + "\t" + rsID + "\t" + "Original" + System.getProperty("line.separator"));
//					
//				}													
			
		}//End of while
			
		System.out.println("FYI numberofLinesWithNoRsId: " + numberofLinesWithNoRsId);
			
		//close
		bufferedReader.close();
		bufferedWriter.close();

		
	} catch (IOException e) {
		// TODO Auto-generated catch block
		e.printStackTrace();
	}

		
	}
	
	
	public static void readAugmentWrite(){
		
		String dataFolder = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\";
		
		//Given InputFile from eliminating common variants
		String RareVariants_ChrName_GLANET1BasedCoordinate_GRCh37p13_rsID_InputFileGLANET = "RareVariants_ChrName_GLANET1BasedCoordinateGRCh37p13_rsID_LGMD-FamB-WES-All_chr_result.tep.txt";
		
		//After Task1 Uplift
		String RareVariants_ChrName_GLANET1BasedCoordinate_LatestAssembly_rsID_InputFileGLANET = "RareVariants_ChrName_GLANET1BasedCoordinateLatestAssembly_rsID_LGMD-FamB-WES-All_chr_result.tep.txt";
		
		//After Task2 Augment
		String RareVariants_ChrName_GLANET1BasedCoordinate_LatestAssembly_rsID_Augmented_InputFileGLANET = "RareVariants_ChrName_GLANET1BasedCoordinateLatestAssembly_rsID_Augmented_LGMD-FamB-WES-All_chr_result.tep.txt";
		
		//After Task3 Downlift
		String RareVariants_ChrName_GLANET1BasedCoordinate_Ch37p13_rsID_Augmented_InputFileGLANET = "RareVariants_ChrName_GLANET1BasedCoordinateGRCh37p13_rsID_Augmented_LGMD-FamB-WES-All_chr_result.tep.txt";

		//Get Latest Assembly right now it is GRCh38.p7
		String latestAssembyNameReturnedByNCBIEutils = NCBIEutils.getLatestAssemblyNameReturnedByNCBIEutils();
		
		String sourceAssemblyName = "GRCh37.p13"; 
		String targetAssemblyName = latestAssembyNameReturnedByNCBIEutils;

		//First Task: Uplift from GRCh37p13 to latest assembly
		Remap.convertGivenInputCoordinatesFromSourceAssemblytoTargetAssemblyUsingRemap(
				dataFolder, 		
				RareVariants_ChrName_GLANET1BasedCoordinate_GRCh37p13_rsID_InputFileGLANET, 
				RareVariants_ChrName_GLANET1BasedCoordinate_LatestAssembly_rsID_InputFileGLANET, 
				sourceAssemblyName, 
				targetAssemblyName,
				true);
		
		//Second Task: Augment with rsID
		readRareVariantsAugmentWithRsIDs(
				dataFolder,
				RareVariants_ChrName_GLANET1BasedCoordinate_LatestAssembly_rsID_InputFileGLANET,
				RareVariants_ChrName_GLANET1BasedCoordinate_LatestAssembly_rsID_Augmented_InputFileGLANET);
		

		//Third Task
		Remap.convertGivenInputCoordinatesFromSourceAssemblytoTargetAssemblyUsingRemap(
				dataFolder, 		
				RareVariants_ChrName_GLANET1BasedCoordinate_LatestAssembly_rsID_Augmented_InputFileGLANET, 
				RareVariants_ChrName_GLANET1BasedCoordinate_Ch37p13_rsID_Augmented_InputFileGLANET, 
				latestAssembyNameReturnedByNCBIEutils, 
				"GRCh37.p13",
				false);
		
		
	}
	
	
	
	
	
	public static void main(String[] args) {
	
		//Some of the given chromosome positions don't have corresponding rsIDs
		//For those cases I will augment each chromosome position with its corresponding rsID		
		readAugmentWrite();
		
	}

}
