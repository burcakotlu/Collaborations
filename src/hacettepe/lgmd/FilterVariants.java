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
public class FilterVariants {

	
	//15 DEC 2016
	public static void filterThatAreNotHomozygotOrCompoundHeterozygotInCase(){
		
		//Input
		String directoryNameandInputFileName = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\LGMD-FamB-WES-All_chr_result.tep.txt";

		//Output
		String directoryNameandFilterVariantsThatAreNotHomozygotOrCompoundHeterozygotInCaseOutpuyFileName = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\FilterVariantsThatAreNotHomozygotOrCompoundHeterozygotInCase.tep.txt";

		FileReader fileReader = null;
		BufferedReader bufferedReader = null;
		
		FileWriter fileWriter_WithAllColumns = null;
		BufferedWriter bufferedWriter_WithAllColumns = null;
		
		String strLine = null;
		
		try {
			
			//Input
			fileReader = FileOperations.createFileReader(directoryNameandInputFileName);		
			bufferedReader = new BufferedReader(fileReader);
		
			//Outputs
			fileWriter_WithAllColumns = FileOperations.createFileWriter(directoryNameandFilterVariantsThatAreNotHomozygotOrCompoundHeterozygotInCaseOutpuyFileName);
			bufferedWriter_WithAllColumns = new BufferedWriter(fileWriter_WithAllColumns);
			
			int indexofTab = -1;
			int indexofFormerTab = -1;
			int count = 0;
			
			String Case_13D0201103_mut = null;
			String Control_13D0201099_mut_father = null;
			String Control_13D0201100_mut_mother = null;		

			
			//Skip header line
			strLine = bufferedReader.readLine();
			
			//Write header line
			bufferedWriter_WithAllColumns.write(strLine + System.getProperty("line.separator")); 
			
			while((strLine = bufferedReader.readLine())!=null){
				
				indexofTab = strLine.indexOf('\t');
				count = 0;
					
				while (indexofTab>0 && count < 14){
					
					count++;
					
					if(count==10){
						Control_13D0201099_mut_father = strLine.substring(indexofFormerTab+1, indexofTab);
					}else if(count==12){
						Control_13D0201100_mut_mother = strLine.substring(indexofFormerTab+1, indexofTab);
					}else if(count==14){
						Case_13D0201103_mut = strLine.substring(indexofFormerTab+1, indexofTab);
					}
					
					indexofFormerTab = indexofTab;
					indexofTab = strLine.indexOf('\t',indexofTab+1);
					
				}//End of WHILE
				
				if ( (Case_13D0201103_mut.contains("Hom") || Case_13D0201103_mut.contains("Het")) && Control_13D0201099_mut_father.contains("Het") && Control_13D0201100_mut_mother.contains("Het")){
					bufferedWriter_WithAllColumns.write(strLine + System.getProperty("line.separator"));
				}
				
			}//End of while read WES data
				
			//Close
			bufferedReader.close();
			bufferedWriter_WithAllColumns.close();

		
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	public static void readMuscularDystrophyDataFilterWriteFile(
			float selectionCriteriaForCommonVariant, 
			float selectionCriteriaForRareVariant){
		
		//Input
		String directoryNameandfileName = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\LGMD-FamB-WES-All_chr_result.tep.txt";
		
		//Output This will be used for RSA for user defined alleles for ENCODE annotated TFs.
		//This one will be used
		String RareNonSynonmousVariants_withAutosomalRecessiceModel_SomeColumns_InputFileGLANET = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\FilterVariants\\RareNonSynonmousVariants_withAutosomalRecessiceModel_SomeColumns_LGMD-FamB-WES-All_chr_result.tep.txt";		
		String RareNonSynonmousVariants_withAutosomalRecessiceModel_AllColumns_InputFileGLANET = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\FilterVariants\\RareNonSynonmousVariants_withAutosomalRecessiceModel_AllColumns_LGMD-FamB-WES-All_chr_result.tep.txt";

		String RareNonSynonmousVariants_SomeColumns_InputFileGLANET = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\FilterVariants\\RareNonSynonmousVariants_SomeColumns_LGMD-FamB-WES-All_chr_result.tep.txt";		
		String RareNonSynonmousVariants_AllColumns_InputFileGLANET = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\FilterVariants\\RareNonSynonmousVariants_AllColumns_LGMD-FamB-WES-All_chr_result.tep.txt";

		FileReader fileReader = null;
		BufferedReader bufferedReader = null;
		
		// with AutosomalRecessiveModel
		FileWriter fileWriter_RareNonSynonmousVariants_withAutosomalRecessiceModel_SomeColumns = null;
		BufferedWriter bufferedWriter_RareNonSynonmousVariants_withAutosomalRecessiceModel_SomeColumns = null;
		
		FileWriter fileWriter_RareNonSynonmousVariants_withAutosomalRecessiceModel_AllColumns = null;
		BufferedWriter bufferedWriter_RareNonSynonmousVariants_withAutosomalRecessiceModel_AllColumns = null;
		
		// without AutosomalRecessiveModel
		FileWriter fileWriter_RareNonSynonmousVariants_SomeColumns = null;
		BufferedWriter bufferedWriter_RareNonSynonmousVariants_SomeColumns = null;
		
		FileWriter fileWriter_RareNonSynonmousVariants_AllColumns = null;
		BufferedWriter bufferedWriter_RareNonSynonmousVariants_AllColumns = null;

		
		String strLine = null;
		
		String chrName = null;
		int _1BasedPosition = -1;
		String reference = null;
		String geneName=null;
		String HGVS = null;
		
		//Chromosome	Position	Reference	GeneName	OMIM	Inheritance	Function	HGVS	ShareNumber	Control_13D0201099_mut	Ration	Control_13D0201100_mut	Ration	Case_13D0201103_mut	Ration	dbSNP_fre	1000human_fre	Hapmap_fre	Agilent_38M_fre	Agilent_46M_fre	Agilent_50M_fre	Nimblegen_44M_fre	Prediction from SIFT	Score from SIFT	RS-ID	NM-ID	Sub-region	Strand	ResidueChange	Gene description	GO_BP	GO_MF	GO_CC	KEGG_Pathway

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
		
		String slashSeparatedObservedAlleles = null;
		
		try {
			
			//Input
			fileReader = FileOperations.createFileReader(directoryNameandfileName);
			bufferedReader = new BufferedReader(fileReader);
			
			//Outputs with AutosomalRecessiveModel
			fileWriter_RareNonSynonmousVariants_withAutosomalRecessiceModel_SomeColumns = FileOperations.createFileWriter(RareNonSynonmousVariants_withAutosomalRecessiceModel_SomeColumns_InputFileGLANET);
			bufferedWriter_RareNonSynonmousVariants_withAutosomalRecessiceModel_SomeColumns = new BufferedWriter(fileWriter_RareNonSynonmousVariants_withAutosomalRecessiceModel_SomeColumns);
			
			fileWriter_RareNonSynonmousVariants_withAutosomalRecessiceModel_AllColumns = FileOperations.createFileWriter(RareNonSynonmousVariants_withAutosomalRecessiceModel_AllColumns_InputFileGLANET);
			bufferedWriter_RareNonSynonmousVariants_withAutosomalRecessiceModel_AllColumns = new BufferedWriter(fileWriter_RareNonSynonmousVariants_withAutosomalRecessiceModel_AllColumns);
			
			//Outputs without AutosomalRecessiveModel
			fileWriter_RareNonSynonmousVariants_SomeColumns = FileOperations.createFileWriter(RareNonSynonmousVariants_SomeColumns_InputFileGLANET);
			bufferedWriter_RareNonSynonmousVariants_SomeColumns = new BufferedWriter(fileWriter_RareNonSynonmousVariants_SomeColumns);
			
			fileWriter_RareNonSynonmousVariants_AllColumns = FileOperations.createFileWriter(RareNonSynonmousVariants_AllColumns_InputFileGLANET);
			bufferedWriter_RareNonSynonmousVariants_AllColumns = new BufferedWriter(fileWriter_RareNonSynonmousVariants_AllColumns);

		
			//Skip header line
			strLine = bufferedReader.readLine();
			
			//Write header line
			bufferedWriter_RareNonSynonmousVariants_withAutosomalRecessiceModel_SomeColumns.write("#chrName" +  "\t" + "_1BasedStartPositionGRCh37.p13" + "\t" + "_1BasedEndPositionGRCh37.p13" + "\t" + "SlashSeparatedObservedAlleles" + "\t" + "Reference" + "\t" + "GeneName" + "\t" + "Function" + "\t" + "HGVS" + "\t" + "Control_13D0201099_mut_father" + "\t" + "Control_13D0201100_mut_mother" + "\t" + "Case_13D0201103_mut" + "\t" + "rsID" +   System.getProperty("line.separator"));
			bufferedWriter_RareNonSynonmousVariants_withAutosomalRecessiceModel_AllColumns.write("#" + strLine + System.getProperty("line.separator")); 
			
			bufferedWriter_RareNonSynonmousVariants_SomeColumns.write("#chrName" +  "\t" + "_1BasedStartPositionGRCh37.p13" + "\t" + "_1BasedEndPositionGRCh37.p13" + "\t" + "SlashSeparatedObservedAlleles" + "\t" + "Reference" + "\t" + "GeneName" + "\t" + "Function" + "\t" + "HGVS" + "\t" + "Control_13D0201099_mut_father" + "\t" + "Control_13D0201100_mut_mother" + "\t" + "Case_13D0201103_mut" + "\t" + "rsID" +   System.getProperty("line.separator"));
			bufferedWriter_RareNonSynonmousVariants_AllColumns.write("#" + strLine + System.getProperty("line.separator")); 

			while((strLine = bufferedReader.readLine())!=null){
				
				//Initialize
				function = null;
				rsID = null;
								
				//Initialize for each line
				//If there is "-" then accept frequency as 0.
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
					}else if(count==3){
						reference = strLine.substring(indexofFormerTab+1, indexofTab);
					}else if(count==4){
						geneName = strLine.substring(indexofFormerTab+1, indexofTab);
					}else if (count==7){
						function = strLine.substring(indexofFormerTab+1, indexofTab);
					}else if (count==8){
						HGVS = strLine.substring(indexofFormerTab+1, indexofTab);
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
						
						//prepareDataforUserDefinedSNPs
						slashSeparatedObservedAlleles = getObservedAlleles(Case_13D0201103_mut);

						if (slashSeparatedObservedAlleles!=null){
							
							bufferedWriter_RareNonSynonmousVariants_SomeColumns.write(chrName + "\t" + _1BasedPosition + "\t" + _1BasedPosition + "\t" + slashSeparatedObservedAlleles + "\t" + reference + "\t" + geneName + "\t" + function + "\t" + HGVS + "\t" + Control_13D0201099_mut_father + "\t" + Control_13D0201100_mut_mother + "\t" + Case_13D0201103_mut + "\t" + rsID +   System.getProperty("line.separator"));
							bufferedWriter_RareNonSynonmousVariants_AllColumns.write(strLine +   System.getProperty("line.separator"));
						
							//Autosomal Recessive Model
							if ( (Case_13D0201103_mut.contains("Hom") || Case_13D0201103_mut.contains("Het")) && Control_13D0201099_mut_father.contains("Het") && Control_13D0201100_mut_mother.contains("Het")){
								

								bufferedWriter_RareNonSynonmousVariants_withAutosomalRecessiceModel_SomeColumns.write(chrName + "\t" + _1BasedPosition + "\t" + _1BasedPosition + "\t" + slashSeparatedObservedAlleles + "\t" + reference + "\t" + geneName + "\t" + function + "\t" + HGVS + "\t" + Control_13D0201099_mut_father + "\t" + Control_13D0201100_mut_mother + "\t" + Case_13D0201103_mut + "\t" + rsID +   System.getProperty("line.separator"));
								bufferedWriter_RareNonSynonmousVariants_withAutosomalRecessiceModel_AllColumns.write(strLine +   System.getProperty("line.separator"));
								
							}

						}//End of IF slashSeparatedObservedAlleles is not NULL
						
						if (subsetFunctionMap.get(function)==null){
							subsetFunctionMap.put(function,1);					
						}else{
							subsetFunctionMap.put(function,subsetFunctionMap.get(function)+1);					
						}
						
						
					
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
			bufferedWriter_RareNonSynonmousVariants_withAutosomalRecessiceModel_SomeColumns.close();
			bufferedWriter_RareNonSynonmousVariants_withAutosomalRecessiceModel_AllColumns.close();
			
			bufferedWriter_RareNonSynonmousVariants_SomeColumns.close();
			bufferedWriter_RareNonSynonmousVariants_AllColumns.close();
			
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	public static String getObservedAlleles(String Case_13D0201103_mut){
		
		int indexOfFirstSemiColon = Case_13D0201103_mut.indexOf(';');
		int indexOfSecondSemiColon = Case_13D0201103_mut.indexOf(';',indexOfFirstSemiColon+1);
		
		String slashSeparatedObservedAlleles = null;
		
		String observedAlleles = null;
		
		if (indexOfSecondSemiColon >= 0){
			//+2AA;V33/W9;Hom
			observedAlleles = Case_13D0201103_mut.substring(0,indexOfFirstSemiColon);
			slashSeparatedObservedAlleles = observedAlleles;
		}else if (indexOfFirstSemiColon >= 0){
			//G90G44A0;Hom
			System.out.println(Case_13D0201103_mut);
			observedAlleles = Case_13D0201103_mut.substring(0,indexOfFirstSemiColon);
			
			int i=0;
			int positionofSecondChar = 0;
			int positionofFirstDigitAfterSecondChar = 0;
			int positionofThirdChar = 0;
			int positionofFirstDigitAfterThirdChar = 0;
			
			//Find the first character
			do{
				
				
				if (Character.isLetter(observedAlleles.charAt(i))){
					//Skip that one
					i++;
					break;
				}				
			}while(i<observedAlleles.length());
			
			//Find the second character
			do{
				if (Character.isLetter(observedAlleles.charAt(i))){
					//Second character is found
					//Keep that one
					positionofSecondChar = i;
					break;
				}
				i++;
			}while(i<observedAlleles.length());
			
			//Find the first digit char after the second character
			do{
				if (Character.isDigit(observedAlleles.charAt(i))){
					positionofFirstDigitAfterSecondChar = i;
					break;
				}
				i++;
			}while(i<observedAlleles.length());
			
			//Find the third character
			do{
				if (Character.isLetter(observedAlleles.charAt(i))){
					//Third character is found
					//Keep that one
					positionofThirdChar = i;
					break;
				}
				i++;
			}while(i<observedAlleles.length());
			
			//Find the first digit char after the third character
			do{
				if (Character.isDigit(observedAlleles.charAt(i))){
					positionofFirstDigitAfterThirdChar = i;
					break;
				}
				i++;
			}while(i<observedAlleles.length());
			
			
			slashSeparatedObservedAlleles = observedAlleles.substring(positionofSecondChar, positionofFirstDigitAfterSecondChar) +"/"  + observedAlleles.substring(positionofThirdChar, positionofFirstDigitAfterThirdChar);

		}
		
		
		return slashSeparatedObservedAlleles;
		
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
		
		//Do I get the same number of SNPs after filter?
		//filterThatAreNotHomozygotOrCompoundHeterozygotInCase();
		
	}

}
