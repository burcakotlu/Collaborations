/**
 * 
 */
package hacettepe.lgmd;

import giveninputdata.InputDataProcess;
import giveninputdata.InputDataRemoveOverlaps;
import giveninputdata.Preparation;
import rsat.GeneAnnotationForPostAnalysisRSAResults;
import rsat.RegulatorySequenceAnalysisPostAnalysis;

/**
 * @author Burçak Otlu
 * @date Dec 23, 2016
 * @project Collaborations 
 *
 */
public class RegulayorySequenceAnalysisForAllTFsForUserDefinedSNPs {

	/**
	 * 
	 */
	public RegulayorySequenceAnalysisForAllTFsForUserDefinedSNPs() {
		// TODO Auto-generated constructor stub
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		
//		/************************ Preparation starts ********************************************/
//		Preparation.main(args);
//		/************************ Preparation ends **********************************************/
//		
//		/************************ InputDataProcess starts ***************************************/
//		InputDataProcess.main(args);
//		/************************ InputDataProcess ends *****************************************/
//
//		/************************ RemoveOverlaps starts ******************************************/
//		InputDataRemoveOverlaps.main(args);
//		/************************ RemoveOverlaps ends ********************************************/
//		
//		
//		/************* Regulatory Sequence Analysis starts ****************************************/					
//		//Generate altered sequences by using user defined observed alleles.
//		//And PFM matrices for all TFs
//		GenerationofSequencesForSNPsandMatricesforAllTFs.main(args);

		RegulatorySequenceAnalysisforAllTFsUsingRSATMatrixScan.main(args);
		
		RegulatorySequenceAnalysisPostAnalysis.main(args);	
		
		GeneAnnotationForPostAnalysisRSAResults.main(args);			
		/************* Regulatory Sequence Analysis ends ******************************************/


	}

}
