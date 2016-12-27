package hacettepe.lgmd;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.rpc.ServiceException;

import org.apache.log4j.Logger;

import rsat.RegulatorySequenceAnalysisUsingRSATMatrixScan;
import ui.GlanetRunner;
import RSATWS.MatrixScanRequest;
import RSATWS.RSATWSPortType;
import RSATWS.RSATWebServicesLocator;
import auxiliary.FileOperations;

import common.Commons;

import enumtypes.CommandLineArguments;




/**
 * @author Burçak Otlu
 * @date Dec 26, 2016
 * @project Collaborations 
 *
 */
public class RegulatorySequenceAnalysisforAllTFsUsingRSATMatrixScan {
	
	final static Logger logger = Logger.getLogger(RegulatorySequenceAnalysisforAllTFsUsingRSATMatrixScan.class);


	public static void matrixScan( 
			String outputFolder, 
			String forRSASNPTFSequencesMatricesDirectory,
			BufferedWriter bufferedWriter) {

		int matrixScanNumber = 1;

		Map<String, String> snpReferenceSequence2RSATResultMap = new HashMap<String, String>();
		Map<String, String> snpAlteredSequence2RSATResultMap = new HashMap<String, String>();
		Map<String, String> tfExtendedPeakSequence2RSATResultMap = new HashMap<String, String>();

		File mainSNPsDirectory = new File(outputFolder + forRSASNPTFSequencesMatricesDirectory + Commons.SNPs);
		File mainTFPFMAndLogoMatricesDirectory = new File(outputFolder + forRSASNPTFSequencesMatricesDirectory + Commons.TF_PFM_AND_LOGO_Matrices);
		
		//26 DEC 2016 stars
		//Another way of doing it
		//TODO read this consideredTFsFile
		//File consideredTFsFile = new File(outputFolder + forRSASNPTFSequencesMatricesDirectory + Commons.CONSIDERED_TFs + Commons.TEXT_FILE_EXTENSION);
		//Then Fill tfNames
		//Fill tfName2TFPfmMatricesFileMap, key is tfName, value is FPfmMatricesFileName (its absolute path)
		//Considered TFs: They can be all TF or all ENCODE TFs		
		List<String> tfNames = new ArrayList<String>();
		Map<String, String> tfName2TFPfmMatricesFileMap = new HashMap<String, String>();
		
		String tfName = null;
		String absolutePath = null;
		
		if(mainTFPFMAndLogoMatricesDirectory.exists() && mainTFPFMAndLogoMatricesDirectory.isDirectory()){			
			for( File eachTFDirectory : mainTFPFMAndLogoMatricesDirectory.listFiles()){
				
				tfName = eachTFDirectory.getName();
				absolutePath = eachTFDirectory.getAbsolutePath() + System.getProperty("file.separator") + Commons.PFM_MATRICES  + Commons.UNDERSCORE + tfName + Commons.TEXT_FILE_EXTENSION;	
				
				tfName2TFPfmMatricesFileMap.put(tfName, absolutePath);
				tfNames.add(tfName);
				
			}//End of for each TF

		}//End of If mainTFPFMAndLogoMatricesDirectory exists 
		//26 DEC 2016 ends

		String snpReferenceSequenceFile = null;
		String snpAlteredSequenceFile = null;
		String tfExtendedPeakSequenceFile = null;
		
		List<String> snpAlteredSequenceFileList = null;		

		String fileName = null;
		String fileAbsolutePath = null;

		RSATWebServicesLocator service = new RSATWebServicesLocator();

		// User the server address for RSAT Metazoa
		service.setRSATWSPortTypeEndpointAddress("http://rsat.sb-roscoff.fr/web_services/RSATWS.cgi");

		RSATWSPortType proxy = null;

		// mainSNPsDirectory is Commons.SNPs
		if(mainSNPsDirectory.exists() && mainSNPsDirectory.isDirectory()){

			try{
				proxy = service.getRSATWSPortType();
			}catch( ServiceException e){
				if( GlanetRunner.shouldLog())
				logger.error( e.toString());
			}

			MatrixScanRequest matrixScanRequest = new MatrixScanRequest();

			/**************************************************************************************************************************/
			/************************* Initialize Matrix Scan Request Parameters starts ***********************************************/
			/**************************************************************************************************************************/
			RegulatorySequenceAnalysisUsingRSATMatrixScan.initializeMatrixScanParameters(matrixScanRequest);
			/**************************************************************************************************************************/
			/************************* Initialize Matrix Scan Request Parameters ends *************************************************/
			/**************************************************************************************************************************/

			// example eachSNPDirectory is chr1_11802721_rs17367504
			for( File eachSNPDirectory : mainSNPsDirectory.listFiles()){

				// Initialize input files
				snpReferenceSequenceFile = null;
				snpAlteredSequenceFile = null;
				tfExtendedPeakSequenceFile = null;

				snpAlteredSequenceFileList = new ArrayList<String>();
				
				fileName = null;
				fileAbsolutePath = null;

				// example eachSNPDirectory chr21_42416281_rs9976767
				if( eachSNPDirectory.isDirectory()){

					// Now get the SNP specific files
					for( File eachSNPFile : eachSNPDirectory.listFiles()){

						fileName = eachSNPFile.getName();
						fileAbsolutePath = eachSNPFile.getAbsolutePath();

						if( fileName.startsWith( Commons.SNP_REFERENCE_SEQUENCE)){
							snpReferenceSequenceFile = fileAbsolutePath;
						}else if( fileName.startsWith( Commons.SNP_ALTERED_SEQUENCE)){
							snpAlteredSequenceFile = fileAbsolutePath;
							snpAlteredSequenceFileList.add( snpAlteredSequenceFile);
						}else if( fileName.startsWith( Commons.TF_EXTENDED_PEAK_SEQUENCE)){							
							//There is only one TF extended peak sequence file for all TFS that overlap with that SNP
							tfExtendedPeakSequenceFile = fileAbsolutePath;
						}

					}// End of FOR eachSNPFile under eachSNPDirectory

					// If all necessary files are not null
					//if( snpReferenceSequenceFile != null && snpAlteredSequenceFileList.size() > 0 && tfName2TfExtendedPeakSequenceFileMap.size() > 0){
					if( snpReferenceSequenceFile != null && snpAlteredSequenceFileList.size() > 0 && tfExtendedPeakSequenceFile!=null){
						
						// Matrix Scan Call
						// what is enrichedElement
						// what is given interval name
						// what is snp
						if( GlanetRunner.shouldLog())
						logger.info("RSAT MatrixScan " + matrixScanNumber++ + " for " + eachSNPDirectory.getPath());
						RegulatorySequenceAnalysisUsingRSATMatrixScan.matrixScan( eachSNPDirectory.getName(), 
								tfName2TFPfmMatricesFileMap, 
								snpReferenceSequenceFile,
								snpAlteredSequenceFileList, 
								tfExtendedPeakSequenceFile,
								tfNames,
								proxy,
								matrixScanRequest, 
								bufferedWriter, 
								snpReferenceSequence2RSATResultMap,
								snpAlteredSequence2RSATResultMap, 
								tfExtendedPeakSequence2RSATResultMap);

					}// End of IF RSAT matrix scan

				}// End of IF eachSNPDirectory is a directory
			}// End of FOR eachSNPDirectory under mainSNPsDirectory

		}// End of IF mainSNPsDirectory exists and mainSNPsDirectory is a directory

	}


	
	public static void main(String[] args) {
		
		// jobName starts
		String jobName = args[CommandLineArguments.JobName.value()].trim();
		if( jobName.isEmpty()){
			jobName = Commons.NO_NAME;
		}
		// jobName ends

		String outputFolder = args[CommandLineArguments.OutputFolder.value()];

		String forRSASNPTFSequencesMatricesDirectory = Commons.FOR_RSA_SNP_TF_SEQUENCES_MATRICES_DIRECTORY;

		FileWriter fileWriter;
		BufferedWriter bufferedWriter;

		try{
			fileWriter = FileOperations.createFileWriter( outputFolder + Commons.REGULATORY_SEQUENCE_ANALYSIS_DIRECTORY + Commons.RSA_RESULTS);
			bufferedWriter = new BufferedWriter( fileWriter);

			/***************************************************************************************************/
			/***************** Regulatory Sequence Analysis for All TF Annotations starts ************************/
			/***************************************************************************************************/
			if( GlanetRunner.shouldLog())
				logger.info( "RSAT starts for TF");
			
				matrixScan(outputFolder, forRSASNPTFSequencesMatricesDirectory, bufferedWriter);
				
			if( GlanetRunner.shouldLog())
				logger.info( "RSAT ends for TF");
			/***************************************************************************************************/
			/***************** Regulatory Sequence Analysis for All TF Annotations ends **************************/
			/***************************************************************************************************/

			// Close bufferedWriter
			bufferedWriter.close();

		}catch( IOException e){

			if( GlanetRunner.shouldLog())
			logger.error( e.toString());
		}
		
	}

}
