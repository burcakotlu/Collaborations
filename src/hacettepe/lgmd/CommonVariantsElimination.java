/**
 * 
 */
package hacettepe.lgmd;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;

import enumtypes.ChromosomeName;
import auxiliary.FileOperations;

/**
 * @author Burçak Otlu
 * @date Nov 21, 2016
 * @project Collaboration 
 * 
 * LGMD stands for Limb-girdle muscular dystrophy.
 * LGMD is a term for a group of diseases that cause weakness and wasting of the muscles in the arms and legs.
 *
 */
public class CommonVariantsElimination {

	
	public static void readMuscularDystrophyDataFilterWriteFile(
			float selectionCriteriaForCommonVariant, 
			float selectionCriteriaForRareVariant){
		
		//Input
		String directoryNameandfileName = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\LGMD-FamB-WES-All_chr_result.tep.txt";

		//Output
		String directoryNameandRareVariantsWithAllColumnsfileName = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\RareVariants_All_Columns_LGMD-FamB-WES-All_chr_result.tep.txt";
		String directoryNameandRareVariants_ChrName_GLANET1BasedCoordinate_InputFileGLANET = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\RareVariants_ChrName_GLANET1BasedCoordinateGRCh37p13_LGMD-FamB-WES-All_chr_result.tep.txt";
		String directoryNameandRareVariants_ChrName_GLANET1BasedCoordinate_rsID_InputFileGLANET = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\RareVariants_ChrName_GLANET1BasedCoordinateGRCh37p13_rsID_LGMD-FamB-WES-All_chr_result.tep.txt";
		String directoryNameandRareVariants_ChrName_GLANET1BasedCoordinate_father_mother_case_rsID_InputFileGLANET = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\RareVariants_ChrName_GLANET1BasedCoordinateGRCh37p13_father_mother_case_rsID_LGMD-FamB-WES-All_chr_result.tep.txt";
		String directoryNameandRareVariants_ChrName_GLANET1BasedCoordinate_rsID__chr4_related_haplotype_InputFileGLANET = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\RareVariants_ChrName_GLANET1BasedCoordinateGRCh37p13_rsID_chr4_related_haplotype__LGMD-FamB-WES-All_chr_result.tep.txt";
		
		FileReader fileReader = null;
		BufferedReader bufferedReader = null;
		
		FileWriter fileWriter_RareVariants_WithAllColumns = null;
		BufferedWriter bufferedWriter_RareVariants_WithAllColumns = null;
		
		FileWriter fileWriter_RareVariants_ChrName_GLANET1BasedCoordinate = null;
		BufferedWriter bufferedWriter_RareVariants_ChrName_GLANET1BasedCoordinate = null;
		
		FileWriter fileWriter_RareVariants_ChrName_GLANET1BasedCoordinate_rsID = null;
		BufferedWriter bufferedWriter_RareVariants_ChrName_GLANET1BasedCoordinate_rsID = null;
		
		FileWriter fileWriter_RareVariants_ChrName_GLANET1BasedCoordinate_father_mother_case_rsID = null;
		BufferedWriter bufferedWriter_RareVariants_ChrName_GLANET1BasedCoordinate_father_mother_case_rsID = null;

		
		FileWriter fileWriter_RareVariants_ChrName_GLANET1BasedCoordinate_rsID_on_chr4_related_haplotype = null;
		BufferedWriter bufferedWriter_RareVariants_ChrName_GLANET1BasedCoordinate_rsID_on_chr4_related_haplotype = null;

		
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
		String Case_13D0201103_mut = null;
		String Control_13D0201099_mut_father = null;
		String Control_13D0201100_mut_mother = null;		
		String rsID = null;
		
		float dbSNP_fre = 0f;
		float _1000human_fre = 0f;
		float Hapmap_fre = 0f;
		float Agilent_38M_fre = 0f;
		float Agilent_46M_fre = 0f;
		float Agilent_50M_fre = 0f;
		float Nimblegen_44M_fre = 0f;
				
		
		Map<String,Integer> functionMap = new HashMap<String,Integer>();
		int numberofALLSNPS = 0;

		Map<String,Integer> subsetFunctionMap = new HashMap<String,Integer>();
		int numberofSubSetofSNPS = 0;
		
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

			fileWriter_RareVariants_ChrName_GLANET1BasedCoordinate_father_mother_case_rsID = FileOperations.createFileWriter(directoryNameandRareVariants_ChrName_GLANET1BasedCoordinate_father_mother_case_rsID_InputFileGLANET);
			bufferedWriter_RareVariants_ChrName_GLANET1BasedCoordinate_father_mother_case_rsID = new BufferedWriter(fileWriter_RareVariants_ChrName_GLANET1BasedCoordinate_father_mother_case_rsID);
			
			fileWriter_RareVariants_ChrName_GLANET1BasedCoordinate_rsID_on_chr4_related_haplotype = FileOperations.createFileWriter(directoryNameandRareVariants_ChrName_GLANET1BasedCoordinate_rsID__chr4_related_haplotype_InputFileGLANET);;
			bufferedWriter_RareVariants_ChrName_GLANET1BasedCoordinate_rsID_on_chr4_related_haplotype = new BufferedWriter(fileWriter_RareVariants_ChrName_GLANET1BasedCoordinate_rsID_on_chr4_related_haplotype);

			
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
					}else if (count==2){
						_1BasedPosition= Integer.parseInt(strLine.substring(indexofFormerTab+1, indexofTab));
					}else if (count==7){
						function = strLine.substring(indexofFormerTab+1, indexofTab);
					}else if(count==10){
						Control_13D0201099_mut_father = strLine.substring(indexofFormerTab+1, indexofTab);
					}else if(count==12){
						Control_13D0201100_mut_mother = strLine.substring(indexofFormerTab+1, indexofTab);
					}else if(count==14){
						Case_13D0201103_mut = strLine.substring(indexofFormerTab+1, indexofTab);
					}
					
					indexofFormerTab = indexofTab;
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
				
				//Initialize for each line
				//If there is "-" then accept frequency as 0.
				dbSNP_fre = 0f;
				_1000human_fre = 0f;
				Hapmap_fre = 0f;
				Agilent_38M_fre = 0f;
				Agilent_46M_fre = 0f;
				Agilent_50M_fre = 0f;
				Nimblegen_44M_fre = 0f;
				
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
	
						if ( (Case_13D0201103_mut.contains("Hom") || Case_13D0201103_mut.contains("Het")) && Control_13D0201099_mut_father.contains("Het") && Control_13D0201100_mut_mother.contains("Het")  ){
							bufferedWriter_RareVariants_ChrName_GLANET1BasedCoordinate_father_mother_case_rsID.write(chrName +  "\t" + _1BasedPosition + "\t" + Control_13D0201099_mut_father + "\t" + Control_13D0201100_mut_mother + "\t" + Case_13D0201103_mut + "\t" + rsID +   System.getProperty("line.separator"));
						}
						
						if (subsetFunctionMap.get(function)==null){
							subsetFunctionMap.put(function,1);					
						}else{
							subsetFunctionMap.put(function,subsetFunctionMap.get(function)+1);					
						}
						
						
						if (chrName.equalsIgnoreCase(ChromosomeName.CHROMOSOME4.convertEnumtoString()) && _1BasedPosition < 70325709  &&  _1BasedPosition> 43884368){
							bufferedWriter_RareVariants_ChrName_GLANET1BasedCoordinate_rsID_on_chr4_related_haplotype.write(chrName +  "\t" + _1BasedPosition + "\t" + rsID +  "\t"  + function+ System.getProperty("line.separator"));
						}//End of if on chr4 haplotype
					
					}//End of if synonmyous SNPs
					
				}//End of rare variants
				
				
				if (functionMap.get(function)==null){
					functionMap.put(function,1);					
				}else{
					functionMap.put(function,functionMap.get(function)+1);					
				}
					
					
				
			}//End of while reading input file
			
			
			//All SNPs
			System.out.println("******************************************************************");
			for(Entry<String, Integer> entry: functionMap.entrySet()){
				System.out.println(entry.getKey() + "\t" + entry.getValue());
				numberofALLSNPS += entry.getValue();
			}
			System.out.println("Number of all SNPS" + "\t" + numberofALLSNPS);
			System.out.println("******************************************************************");
			
			//Rare not synnonmous SNPs
			System.out.println("******************************************************************");
			for(Entry<String, Integer> entry: subsetFunctionMap.entrySet()){
				System.out.println(entry.getKey() + "\t" + entry.getValue());
				numberofSubSetofSNPS += entry.getValue();
			}
			System.out.println("Number of rare not synnonmous SNPS" + "\t" + numberofSubSetofSNPS);
			System.out.println("******************************************************************");

			
			
			//Close
			bufferedReader.close();
			bufferedWriter_RareVariants_WithAllColumns.close();
			bufferedWriter_RareVariants_ChrName_GLANET1BasedCoordinate.close();
			bufferedWriter_RareVariants_ChrName_GLANET1BasedCoordinate_rsID.close();
			bufferedWriter_RareVariants_ChrName_GLANET1BasedCoordinate_rsID_on_chr4_related_haplotype.close();
			bufferedWriter_RareVariants_ChrName_GLANET1BasedCoordinate_father_mother_case_rsID.close();
			
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
