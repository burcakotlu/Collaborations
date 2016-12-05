package auxiliary;

import enumtypes.ChromosomeName;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class HG38_RefSeq_Genes {

	

	public static void readUCSCGenomeHG38RefSeqGenes(
			String UCSC_GENOME_HG38_REFSEQ_GENES_FILE,
			TObjectIntMap<String> geneSymbol2GeneInternalNumberMap,
			TIntObjectMap<HG38RefSeqGeneInformation> geneInternalNumber2GeneInformationMap){


		// GLANET convention is 0_based_start and 0 _based_end
		// UCSC GENOME table browser convention: Our internal database
		// representations of coordinates always have a zero-based start and a one-based end.
		// Therefore convert 1_based_end to 0_based_end

		FileReader fileReader = null;
		BufferedReader bufferedReader = null;

		String strLine;
		//#bin	name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	score	name2	cdsStartStat	cdsEndStat	exonFrames

		String refSeqGeneName; //name
		ChromosomeName chromName;
		char strand;
		int txStart;
		int txEnd;
		int exonCounts;
		String exonStarts;
		String exonEnds;
		String geneSymbol; //name2

		TIntList exonStartList = null;
		TIntList exonEndList = null;

		int indexofFormerComma = -1;
		int indexofLatterComma = -1;

		int indexofFirstTab = -1;
		int indexofSecondTab = -1;
		int indexofThirdTab = -1;
		int indexofFourthTab = -1;
		int indexofFifthTab = -1;
		int indexofSixthTab = -1;
		int indexofSeventhTab  = -1;
		int indexofEigthTab = -1;
		int indexofNinethTab = -1;
		int indexofTenthTab = -1;
		int indexofEleventhTab = -1;
		int indexofTwelfthTab = -1;
		int indexofThirteenthTab = -1;
		
		int geneNumber= 1;

		HG38RefSeqGeneInformation gene = null;
		
		try{

			fileReader = FileOperations.createFileReader(UCSC_GENOME_HG38_REFSEQ_GENES_FILE);
			bufferedReader = new BufferedReader( fileReader);

			// skip headerLine
			//#bin	name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	score	name2	cdsStartStat	cdsEndStat	exonFrames

			strLine = bufferedReader.readLine();

			while( ( strLine = bufferedReader.readLine()) != null){

				// example strLine
				//125	NM_004898	chr4	-	55427900	55546909	55435414	55482785	23	55427900,55438281,55442431,55443686,55444632,55448778,55449395,55450090,55453053,55453676,55455896,55456217,55458891,55459147,55463684,55470716,55475962,55478814,55479639,55482738,55489373,55509911,55546781,	55435594,55438537,55442634,55443896,55444785,55448868,55449496,55450232,55453129,55453824,55456003,55456300,55459010,55459261,55463805,55470806,55476054,55478963,55479699,55482828,55489465,55510065,55546909,	0	CLOCK	cmpl	cmpl	0,2,0,0,0,0,1,0,2,1,2,0,1,1,0,0,1,2,2,0,-1,-1,-1,
				//122	NM_000232	chr4	-	52020694	52038319	52023956	52038259	6	52020694,52027967,52028729,52029677,52033430,52038226,	52024160,52028099,52028921,52029863,52033640,52038319,	0	SGCB	cmpl	cmpl	0,0,0,0,0,0,

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
				chromName = ChromosomeName.convertStringtoEnum( strLine.substring( indexofSecondTab + 1,indexofThirdTab));

				strand = strLine.substring( indexofThirdTab + 1, indexofFourthTab).trim().charAt(0);

				txStart = Integer.parseInt( strLine.substring( indexofFourthTab + 1, indexofFifthTab));
				// Convert one based end to zero based end
				txEnd = Integer.parseInt( strLine.substring( indexofFifthTab + 1, indexofSixthTab)) - 1;

				// cdsStart =Integer.parseInt(strLine.substring(indexofSixthTab+1,indexofSeventhTab));
				// cdsEnd =Integer.parseInt(strLine.substring(indexofSeventhTab+1,indexofEigthTab))-1;

				// 28FEB2014
				// Convert one based end to zero based end
				// Or 0Based EndExclusive to End
				// How do we know that ends are OBased exclusive?
				// Is it written in a readme file?

				exonCounts = Integer.parseInt( strLine.substring( indexofEigthTab + 1, indexofNinethTab));
				exonStarts = strLine.substring( indexofNinethTab + 1, indexofTenthTab);
				exonEnds = strLine.substring( indexofTenthTab + 1, indexofEleventhTab);
				geneSymbol = strLine.substring( indexofTwelfthTab + 1, indexofThirteenthTab);

				// For each strLine
				exonStartList = new TIntArrayList( exonCounts);
				exonEndList = new TIntArrayList( exonCounts);

				//Start filling exonStartList
				// Initialize before for loop
				indexofFormerComma = 0;
				indexofLatterComma = 0;
				
				// Fill exonStartList
				for( int i = 0; i < exonCounts; i++){
					indexofFormerComma = indexofLatterComma;
					indexofLatterComma = exonStarts.indexOf(',', indexofFormerComma+1);
					
					//First exonStart
					if (i==0){
						exonStartList.add(Integer.parseInt(exonStarts.substring(0,indexofLatterComma)));
					}else{
						exonStartList.add(Integer.parseInt(exonStarts.substring(indexofFormerComma+1,indexofLatterComma)));
					}
	
				}// End of for: filling exonStartList
				
				
				//Then fill exonEndList
				// Initialize before for loop
				indexofFormerComma = 0;
				indexofLatterComma = 0;
				
				// Fill exonEndList
				for( int i = 0; i < exonCounts; i++){
					indexofFormerComma = indexofLatterComma;
					indexofLatterComma = exonEnds.indexOf(',', indexofFormerComma+1);
					
					//First exonEnd
					if (i==0){
						exonEndList.add(Integer.parseInt(exonEnds.substring(0,indexofLatterComma)) - 1);
					}else{						
						// 28FEB2014
						// Convert 1_based_end to 0_based_end
						exonEndList.add(Integer.parseInt(exonEnds.substring(indexofFormerComma+1,indexofLatterComma)) - 1);
					}
	
				}// End of for: filling exonStartList and exonEndList
				
				geneSymbol2GeneInternalNumberMap.put(geneSymbol,geneNumber);
				
				gene = new HG38RefSeqGeneInformation();
				
				gene.setRefSeqGeneName(refSeqGeneName);
				gene.setGeneSymbol(geneSymbol);
				gene.setChromName(chromName);
				gene.setTxStart(txStart);
				gene.setTxEnd(txEnd);
				gene.setExonCounts(exonCounts);
				gene.setExonStartList(exonStartList);
				gene.setExonEndList(exonEndList);
				gene.setStrand(strand);
				
				geneInternalNumber2GeneInformationMap.put(geneNumber,gene);
				
				geneNumber++;
	
			}// End of while : each read strLine

			bufferedReader.close();			

		}catch( IOException e){
			e.printStackTrace();
		}

	}
	
}
