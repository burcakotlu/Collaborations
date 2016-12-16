/**
 * 
 */
package hacettepe.lgmd;

import giveninputdata.InputDataProcess;
import giveninputdata.InputDataRemoveOverlaps;
import giveninputdata.Preparation;

import org.apache.log4j.Logger;

import rsat.GeneAnnotationForPostAnalysisRSAResults;
import rsat.GenerationofAllTFAnnotationsFileInGRCh37p13AndInLatestAssembly;
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
		
		/************************ Preparation starts ********************************************/
		Preparation.main(args);
		/************************ Preparation ends **********************************************/
		
		/************************ InputDataProcess starts ***************************************/
		InputDataProcess.main(args);
		/************************ InputDataProcess ends *****************************************/

		/************************ RemoveOverlaps starts ******************************************/
		InputDataRemoveOverlaps.main( args);
		/************************ RemoveOverlaps ends ********************************************/
		
		/************************ Annotation starts ***********************************************/
		Annotation.main(args);
		/************************ Annotation ends *************************************************/

		/************* Regulatory Sequence Analysis starts ****************************************/					
		//part1
		//I might change this one, simplify the conversion
		GenerationofAllTFAnnotationsFileInGRCh37p13AndInLatestAssembly.main(args);

		//I have to change the way I generate altered sequences by using user defined observed alleles.
		///This might require the way I handle part1
		GenerationofSequencesandMatricesforSNPs.main(args);

		RegulatorySequenceAnalysisUsingRSATMatrixScan.main(args);
		
		RegulatorySequenceAnalysisPostAnalysis.main(args);
		
		GeneAnnotationForPostAnalysisRSAResults.main(args);			
		/************* Regulatory Sequence Analysis ends ******************************************/
		
	}

}
