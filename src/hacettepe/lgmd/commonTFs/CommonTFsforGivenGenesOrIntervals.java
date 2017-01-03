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

import auxiliary.FileOperations;
import collaboration.Interval;

import common.Commons;

import enumtypes.ChromosomeName;

/**
 * @author Burçak Otlu
 * @date Jan 2, 2017
 * @project Collaborations 
 *
 */
public class CommonTFsforGivenGenesOrIntervals {

	/**
	 * 
	 */
	public CommonTFsforGivenGenesOrIntervals() {
		// TODO Auto-generated constructor stub
	}

	
	public static void createIntervals(List<String> lgmdGeneSymbolList){
		
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

		String refSeqGeneName;
		ChromosomeName chromName;
		char strand;
		int txStart;
		int txEnd;		
		String alternateGeneName;
		
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

				refSeqGeneName = strLine.substring( indexofFirstTab + 1, indexofSecondTab);
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
					}
				}
				
				


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
		Map<String,Interval> geneSymbol2IntervalMap = new HashMap<String,Interval>();
		createIntervals(lgmdGeneSymbolList);
		
		//Task2
		//Get DNA sequences for these intervals for the latest assembly
		
		//Get TF matrices
		
		//Call RSAT web service
		
		//Write the best TF matches for each gene 
		
	}

}
