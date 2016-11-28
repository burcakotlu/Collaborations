/**
 * 
 */
package hacettepe.musculardystrophy;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import auxiliary.FileOperations;

/**
 * @author Burçak Otlu
 * @date Nov 21, 2016
 * @project Collaboration 
 *
 */
public class CommonVariantsElimination {

	
	public static void readMuscularDystrophyDataFilterWriteFile(float selectionCriteriaForCommonVariant, float selectionCriteriaForRareVariant){
		
		//Input
		String directoryNameandfileName = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\LGMD-FamB-WES-All_chr_result.tep.txt";

		//Output
		String directoryNameandRareVariantsWithAllColumnsfileName = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\RareVariants_All_Columns_LGMD-FamB-WES-All_chr_result.tep.txt";
		String directoryNameandRareVariants_ChrName_GLANET1BasedCoordinate_InputFileGLANET = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\RareVariants_ChrName_GLANET1BasedCoordinateGRCh37p13_LGMD-FamB-WES-All_chr_result.tep.txt";
		String directoryNameandRareVariants_ChrName_GLANET1BasedCoordinate_rsID_InputFileGLANET = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\RareVariants_ChrName_GLANET1BasedCoordinateGRCh37p13_rsID_LGMD-FamB-WES-All_chr_result.tep.txt";

		
		FileReader fileReader = null;
		BufferedReader bufferedReader = null;
		
		FileWriter fileWriter_RareVariants_WithAllColumns = null;
		BufferedWriter bufferedWriter_RareVariants_WithAllColumns = null;
		
		FileWriter fileWriter_RareVariants_ChrName_GLANET1BasedCoordinate = null;
		BufferedWriter bufferedWriter_RareVariants_ChrName_GLANET1BasedCoordinate = null;
		
		FileWriter fileWriter_RareVariants_ChrName_GLANET1BasedCoordinate_rsID = null;
		BufferedWriter bufferedWriter_RareVariants_ChrName_GLANET1BasedCoordinate_rsID = null;
		
		String strLine = null;
		
		String chrName = null;
		int _1BasedPosition = -1;
		
		//Chromosome	Position	Reference	GeneName	OMIM	Inheritance	Function	HGVS	ShareNumber	Control_13D0201099_mut	Ration	Control_13D0201100_mut	Ration	Case_13D0201103_mut	Ration	dbSNP_fre	1000human_fre	Hapmap_fre	Agilent_38M_fre	Agilent_46M_fre	Agilent_50M_fre	Nimblegen_44M_fre	Prediction from SIFT	Score from SIFT	RS-ID	NM-ID	Sub-region	Strand	ResidueChange	Gene description	GO_BP	GO_MF	GO_CC	KEGG_Pathway

		//15 dbSNP_fre 16
		// 1000human_fre
		// Hapmap_fre 
		// Agilent_38M_fre 
		// Agilent_46M_fre 20	
		// Agilent_50M_fre 21	
		// Nimblegen_44M_fre 22
		int indexofSixteenthTab;
		int indexofSeventeenthTab;
		int indexofEighteenthTab;
		int indexofNineteenthTab;
		int indexofTwenythTab;
		int indexofTwentyFirstTab;
		int indexofTwentySecondTab;
		int indexofTwentyThirdTab;
		int indexofTwentyFourthTab;
		int indexofTwentyFifthTab;
				
		int indexofTab = -1;
		int indexofFormerTab = -1;
		int count = 0;
		
		String function = null;
		String rsID = null;
		
		float dbSNP_fre = 0f;
		float _1000human_fre = 0f;
		float Hapmap_fre = 0f;
		float Agilent_38M_fre = 0f;
		float Agilent_46M_fre = 0f;
		float Agilent_50M_fre = 0f;
		float Nimblegen_44M_fre = 0f;
				
		try {
			
			//Input
			fileReader = FileOperations.createFileReader(directoryNameandfileName);
			bufferedReader = new BufferedReader(fileReader);
			
			//Outputs
			fileWriter_RareVariants_WithAllColumns = FileOperations.createFileWriter(directoryNameandRareVariantsWithAllColumnsfileName);
			bufferedWriter_RareVariants_WithAllColumns = new BufferedWriter(fileWriter_RareVariants_WithAllColumns);
	
			fileWriter_RareVariants_ChrName_GLANET1BasedCoordinate = FileOperations.createFileWriter(directoryNameandRareVariants_ChrName_GLANET1BasedCoordinate_InputFileGLANET);
			bufferedWriter_RareVariants_ChrName_GLANET1BasedCoordinate = new BufferedWriter(fileWriter_RareVariants_ChrName_GLANET1BasedCoordinate);
			
			fileWriter_RareVariants_ChrName_GLANET1BasedCoordinate_rsID = FileOperations.createFileWriter(directoryNameandRareVariants_ChrName_GLANET1BasedCoordinate_rsID_InputFileGLANET);
			bufferedWriter_RareVariants_ChrName_GLANET1BasedCoordinate_rsID = new BufferedWriter(fileWriter_RareVariants_ChrName_GLANET1BasedCoordinate_rsID);

			//Skip header line
			strLine = bufferedReader.readLine();
			
			//Write header line
			bufferedWriter_RareVariants_WithAllColumns.write(strLine + System.getProperty("line.separator")); 
			
			while((strLine = bufferedReader.readLine())!=null){
				
				//Initialize
				function = null;
				rsID = null;
				
				dbSNP_fre = 0f;
				_1000human_fre = 0f;
				Hapmap_fre = 0f;
				Agilent_38M_fre = 0f;
				Agilent_46M_fre = 0f;
				Agilent_50M_fre = 0f;
				Nimblegen_44M_fre = 0f;
				count = 0;
				
				indexofTab = strLine.indexOf('\t');
								
				while (indexofTab>0 && count < 14){
					
					count++;
					
					if (count==1){
						chrName = strLine.substring(0, indexofTab);
						indexofFormerTab = indexofTab;
					}else if (count==2){
						_1BasedPosition= Integer.parseInt(strLine.substring(indexofFormerTab+1, indexofTab));
						indexofFormerTab = indexofTab;
					}else if (count==7){
						function = strLine.substring(indexofFormerTab+1, indexofTab);
					}else{
						indexofFormerTab = indexofTab;
					}
					
					indexofTab = strLine.indexOf('\t',indexofTab+1);
				}//End of WHILE
				
				indexofSixteenthTab = strLine.indexOf('\t',indexofTab+1);
				indexofSeventeenthTab  = strLine.indexOf('\t',indexofSixteenthTab+1);
				indexofEighteenthTab  = strLine.indexOf('\t',indexofSeventeenthTab+1);
				indexofNineteenthTab  = strLine.indexOf('\t',indexofEighteenthTab+1);
				indexofTwenythTab  = strLine.indexOf('\t',indexofNineteenthTab+1);
				indexofTwentyFirstTab  = strLine.indexOf('\t',indexofTwenythTab+1);
				indexofTwentySecondTab  = strLine.indexOf('\t',indexofTwentyFirstTab+1);
				indexofTwentyThirdTab  = strLine.indexOf('\t',indexofTwentySecondTab+1);
				indexofTwentyFourthTab  = strLine.indexOf('\t',indexofTwentyThirdTab+1);
				indexofTwentyFifthTab  = strLine.indexOf('\t',indexofTwentyFourthTab+1);
				
				//dbSNP_fre
				if (!strLine.substring(indexofTab+1, indexofSixteenthTab).equalsIgnoreCase("-")){
					dbSNP_fre = Float.parseFloat(strLine.substring(indexofTab+1, indexofSixteenthTab));					
				}
				
				//_1000human_fre
				if (!strLine.substring(indexofSixteenthTab+1, indexofSeventeenthTab).equalsIgnoreCase("-")){
					_1000human_fre = Float.parseFloat(strLine.substring(indexofSixteenthTab+1, indexofSeventeenthTab));
				}
				
				//Hapmap_fre
				if (!strLine.substring(indexofSeventeenthTab+1, indexofEighteenthTab).equalsIgnoreCase("-")){
					Hapmap_fre = Float.parseFloat(strLine.substring(indexofSeventeenthTab+1, indexofEighteenthTab));
				}
				
				//Agilent_38M_fre
				if (!strLine.substring(indexofEighteenthTab+1, indexofNineteenthTab).equalsIgnoreCase("-")){
					Agilent_38M_fre = Float.parseFloat(strLine.substring(indexofEighteenthTab+1, indexofNineteenthTab));
				}
				
				//Agilent_46M_fre
				if (!strLine.substring(indexofNineteenthTab+1, indexofTwenythTab).equalsIgnoreCase("-")){
					Agilent_46M_fre = Float.parseFloat(strLine.substring(indexofNineteenthTab+1, indexofTwenythTab));
				}
				
				//Agilent_50M_fre
				if (!strLine.substring(indexofTwenythTab+1, indexofTwentyFirstTab).equalsIgnoreCase("-")){
					Agilent_50M_fre = Float.parseFloat(strLine.substring(indexofTwenythTab+1, indexofTwentyFirstTab));
				}
				
				//Nimblegen_44M_fre
				if (!strLine.substring(indexofTwentyFirstTab+1, indexofTwentySecondTab).equalsIgnoreCase("-")){
					Nimblegen_44M_fre = Float.parseFloat(strLine.substring(indexofTwentyFirstTab+1, indexofTwentySecondTab));
				}
				
				//RS-ID between 24 and 25
				rsID = strLine.substring(indexofTwentyFourthTab+1, indexofTwentyFifthTab);
				
				
				//Kent and colleagues defined rare variants by MAF < 0.01, but they defined common variants by MAF > 0.1.
				//Choose only Rare Variants
				if ( 	(dbSNP_fre<selectionCriteriaForRareVariant) && (_1000human_fre<selectionCriteriaForRareVariant)  && (Hapmap_fre<selectionCriteriaForRareVariant) &&
						(Agilent_38M_fre<selectionCriteriaForRareVariant) && (Agilent_46M_fre<selectionCriteriaForRareVariant) && (Agilent_50M_fre<selectionCriteriaForRareVariant) &&
						(Nimblegen_44M_fre<selectionCriteriaForRareVariant)){
					
					//Filter Synonymous SNPs
					if (!function.startsWith("Synonymous")){
						bufferedWriter_RareVariants_WithAllColumns.write(strLine + System.getProperty("line.separator")); 
						bufferedWriter_RareVariants_ChrName_GLANET1BasedCoordinate.write(chrName +  "\t" + _1BasedPosition + System.getProperty("line.separator"));
						bufferedWriter_RareVariants_ChrName_GLANET1BasedCoordinate_rsID.write(chrName +  "\t" + _1BasedPosition + "\t" + rsID +   System.getProperty("line.separator"));						
					}
					
				}
				
				
			}//End of while reading input file
			
			//Close
			bufferedReader.close();
			bufferedWriter_RareVariants_WithAllColumns.close();
			bufferedWriter_RareVariants_ChrName_GLANET1BasedCoordinate.close();
			bufferedWriter_RareVariants_ChrName_GLANET1BasedCoordinate_rsID.close();
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// Read the musculardystrophy data
		// Filter the common variants
		// Write the remaining rare variants
		
		float selectionCriteriaForCommonVariant = 0.1f;
		float selectionCriteriaForRareVariant = 0.01f;
		readMuscularDystrophyDataFilterWriteFile(selectionCriteriaForCommonVariant,selectionCriteriaForRareVariant);
		

	}

}
