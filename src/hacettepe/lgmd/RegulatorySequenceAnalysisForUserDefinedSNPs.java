/**
 * 
 */
package hacettepe.lgmd;

import giveninputdata.InputDataProcess;
import giveninputdata.InputDataRemoveOverlaps;
import giveninputdata.Preparation;

import org.apache.log4j.Logger;

import rsat.GeneAnnotationForPostAnalysisRSAResults;
import rsat.GenerationofAllTFAnnotationsFileInGRCh37p13AndInLatestAssemblySimplified;
import rsat.GenerationofSequencesandMatricesforSNPs;
import rsat.RegulatorySequenceAnalysisPostAnalysis;
import rsat.RegulatorySequenceAnalysisUsingRSATMatrixScan;
import annotation.Annotation;

/**
 * @author Burçak Otlu
 * @date Dec 16, 2016
 * @project Collaborations 
 *
 */
public class RegulatorySequenceAnalysisForUserDefinedSNPs {
	
	final static Logger logger = Logger.getLogger(RegulatorySequenceAnalysisForUserDefinedSNPs.class);


	public static void main(String[] args) {
		
		//This run uses input file "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\RareVariants_SomeColumns_LGMD-FamB-WES-All_chr_result.tep.txt"
		//#chrName	_1BasedStartPositionGRCh37.p13	_1BasedEndPositionGRCh37.p13	SlashSeparatedObservedAlleles	Reference	GeneName	Function	HGVS	Control_13D0201099_mut_father	Control_13D0201100_mut_mother	Case_13D0201103_mut	rsID
		//First four columns are required for this run.
		//Existence of fifth or more columns are optional. They may exists or may not. Does not matter.
		
		/************************ Preparation starts ********************************************/
		Preparation.main(args);
		/************************ Preparation ends **********************************************/
		
		/************************ InputDataProcess starts ***************************************/
		InputDataProcess.main(args);
		/************************ InputDataProcess ends *****************************************/

		/************************ RemoveOverlaps starts ******************************************/
		InputDataRemoveOverlaps.main(args);
		/************************ RemoveOverlaps ends ********************************************/
		
		/************************ Annotation starts ***********************************************/
		Annotation.main(args);
		/************************ Annotation ends *************************************************/

		/************* Regulatory Sequence Analysis starts ****************************************/					
		//part1
		//DONE Simplify the conversion
		GenerationofAllTFAnnotationsFileInGRCh37p13AndInLatestAssemblySimplified.main(args);

		//part2
		//DONE Generate altered sequences by using user defined observed alleles.
		GenerationofSequencesandMatricesforSNPs.main(args);

		RegulatorySequenceAnalysisUsingRSATMatrixScan.main(args);
		
		RegulatorySequenceAnalysisPostAnalysis.main(args);	
		
		GeneAnnotationForPostAnalysisRSAResults.main(args);			
		/************* Regulatory Sequence Analysis ends ******************************************/
		
	}

}
