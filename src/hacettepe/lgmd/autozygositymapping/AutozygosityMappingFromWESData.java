/**
 * 
 */
package hacettepe.lgmd.autozygositymapping;

import intervaltree.IntervalTree;
import intervaltree.IntervalTreeNode;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import auxiliary.FileOperations;
import enumtypes.ChromosomeName;

/**
 * @author Burçak Otlu
 * @date Dec 26, 2016
 * @project Collaborations 
 *
 */
public class AutozygosityMappingFromWESData {
	
	public static List<IntervalTreeNode> createIntervals( 
			IntervalTreeNode intervalToBeChanged,
			IntervalTreeNode overlappingInterval) {

		List<IntervalTreeNode> createdIntervals = new ArrayList<IntervalTreeNode>();

		/**************************************************************************************/
		// Left Interval
		IntervalTreeNode leftInterval = null;

		int xLow = intervalToBeChanged.getLow();
		int yLow = overlappingInterval.getLow();

		if( (xLow < yLow)){
			leftInterval = new IntervalTreeNode(intervalToBeChanged.getChromName(), xLow, yLow - 1);
		}
		/**************************************************************************************/

		/**************************************************************************************/
		// Right Interval
		IntervalTreeNode rightInterval = null;

		int xHigh = intervalToBeChanged.getHigh();
		int yHigh = overlappingInterval.getHigh();

		if( (xHigh > yHigh)){
			rightInterval = new IntervalTreeNode(intervalToBeChanged.getChromName(), yHigh + 1, xHigh);
		}
		/**************************************************************************************/

		if( leftInterval != null){
			createdIntervals.add(leftInterval);
		}

		if( rightInterval != null){
			createdIntervals.add(rightInterval);
		}

		return createdIntervals;

	}
	
	public static void writeCaseIntervalTree(Map<ChromosomeName,IntervalTree> case_chromosomeName2IntervalTreeMap,
			BufferedWriter case_AfterOverlapsRemoved_bufferedWriter_LROHs){
		
		IntervalTree intervalTree = null;
		
		for(Map.Entry<ChromosomeName, IntervalTree> entry: case_chromosomeName2IntervalTreeMap.entrySet()) {
			
			intervalTree = entry.getValue();
			
			if (intervalTree!=null && intervalTree.getRoot()!=null && intervalTree.getRoot().getNodeName().isNotSentinel()){
				IntervalTree.intervalTreeInfixTraversal(intervalTree.getRoot(),case_AfterOverlapsRemoved_bufferedWriter_LROHs);
			}
			
		}
		
	}
	
	
	public static void removeOverlapsFromCaseUsingControlAutozygotRegions(
			Map<ChromosomeName,IntervalTree> case_chromosomeName2IntervalTreeMap,
			List<IntervalTreeNode> control_IntervalTreeNodeList){
		
		IntervalTreeNode controlIntervalTreeNode  = null;
		IntervalTree caseIntervalTree = null;
		List<IntervalTreeNode> overlappedNodeList = null;
		
		List<IntervalTreeNode> createdIntervals = null;
		
		for(Iterator<IntervalTreeNode> itr= control_IntervalTreeNodeList.iterator();itr.hasNext();){
			
			controlIntervalTreeNode = itr.next();
			
			//get the corresponding case Interval Tree
			caseIntervalTree = case_chromosomeName2IntervalTreeMap.get(controlIntervalTreeNode.getChromName());
			
			overlappedNodeList = new ArrayList<IntervalTreeNode>();
			
			IntervalTreeNode caseOverlappedNode = null;
			
			if(caseIntervalTree!=null){
				
				//Find the case overlapping IntervalTreeNodes with the control IntervalTreeNode
				caseIntervalTree.findAllOverlappingIntervals(overlappedNodeList, caseIntervalTree.getRoot(),controlIntervalTreeNode);

				// there is overlap
				if( overlappedNodeList != null && overlappedNodeList.size() > 0){

					
					IntervalTreeNode splicedoutNode = null;
					IntervalTreeNode nodetoBeDeleted = null;
					// you may try to delete a node which is already spliced out by former deletions
					// therefore you must keep track of the real node to be deleted 
					// in case of trial of deletion of an already spliced out node.
					Map<IntervalTreeNode, IntervalTreeNode> splicedoutNode2CopiedNodeMap = new HashMap<IntervalTreeNode, IntervalTreeNode>();

					for(int i=0; i<overlappedNodeList.size(); i++){

						caseOverlappedNode = overlappedNodeList.get(i);
						
						//TODO
						//Remove the overlap
						//There can be one or two remaining intervals left
						//updateMergedNode(mergedNode, caseOverlappedNode);
						createdIntervals = createIntervals(caseOverlappedNode,controlIntervalTreeNode);
						

						// if the to be deleted, intended interval tree node
						// is an already spliced out node
						// in other words if it is copied into another
						// interval tree node
						// then you have to delete that node
						// not the already spliced out node

						nodetoBeDeleted = IntervalTree.compute(splicedoutNode2CopiedNodeMap, caseOverlappedNode);

						if( nodetoBeDeleted != null){
							// they are the same
							// current overlapped node has been copied to
							// the previously deleted overlapped node
							// current overlapped node is spliced out by the
							// previous delete operation
							// so delete that previously deleted overlapped
							// node in order to delete the current
							// overlapped node
							// since current overlapped node is copied to
							// the previously deleted overlapped node
							// Now we can delete this overlappedNode
							splicedoutNode = caseIntervalTree.intervalTreeDelete(caseIntervalTree, nodetoBeDeleted);

							if( splicedoutNode != nodetoBeDeleted)
								splicedoutNode2CopiedNodeMap.put( splicedoutNode, nodetoBeDeleted);
						}else{
							// Now we can delete this overlappedNode
							splicedoutNode = caseIntervalTree.intervalTreeDelete(caseIntervalTree, caseOverlappedNode);

							if( splicedoutNode != caseOverlappedNode)
								splicedoutNode2CopiedNodeMap.put( splicedoutNode, caseOverlappedNode);

						}

					}//End of FOR: each case overlapping IntervalTreeNode
					
					//Add created nodes
					for(Iterator<IntervalTreeNode> itr2 =createdIntervals.iterator();itr2.hasNext();){						
						caseIntervalTree.intervalTreeInsert(caseIntervalTree, itr2.next());
					}//End of for each created node

				}//End of IF overlappedNodeList is not null

			}//End of caseIntervalTree is not NULL


		}//End of for each control IntervalTreeNode in the list
		
	}

	
	
	public static void insert(
			ChromosomeName chromosomeName,
			IntervalTreeNode intervalTreeNode,
			Map<ChromosomeName,IntervalTree> chromosomeName2IntervalTreeMap){
		
		IntervalTree intervalTree = chromosomeName2IntervalTreeMap.get(chromosomeName);
		
		if( intervalTree == null){
			intervalTree = new IntervalTree();
			intervalTree.intervalTreeInsert(intervalTree, intervalTreeNode);
			chromosomeName2IntervalTreeMap.put(chromosomeName, intervalTree);
		}else{
			intervalTree.intervalTreeInsert(intervalTree, intervalTreeNode);
		}
		
	}
	
	
	public static void isHomozygotRegionFound(
			int numberofConsecutiveHomozygotsFound,
			int numberofConsecutiveHomVariantsRequired,			
			String chrName,
			ChromosomeName chromosomeName,
			int saved1BasedStart,
			int lastSaved1BasedStart,
			BufferedWriter bufferedWriter,
			Map<ChromosomeName,IntervalTree> chrName2IntervalTreeMap) throws IOException{
		
		int lengthofHomozygotRegion = -1;
		IntervalTreeNode intervalTreeNode = null;
		
		if (numberofConsecutiveHomozygotsFound>= numberofConsecutiveHomVariantsRequired){
			
			lengthofHomozygotRegion = lastSaved1BasedStart-saved1BasedStart+1;
			
			//Output the homozygot region
			bufferedWriter.write(chrName + "\t" + saved1BasedStart + "\t" + lastSaved1BasedStart  + "\t" + lengthofHomozygotRegion + System.getProperty("line.separator"));
			
			//Construct Interval Tree
			intervalTreeNode = new IntervalTreeNode(chromosomeName,saved1BasedStart,lastSaved1BasedStart);								
			insert(chromosomeName,intervalTreeNode,chrName2IntervalTreeMap);
			
		}
		
	}
	
	public static void isHomozygotRegionFound(
			int numberofConsecutiveHomozygotsFound,
			int numberofConsecutiveHomVariantsRequired,			
			String chrName,
			ChromosomeName chromosomeName,
			int saved1BasedStart,
			int lastSaved1BasedStart,
			BufferedWriter bufferedWriter,
			List<IntervalTreeNode> intervalTreeNodeList) throws IOException{
		
		int lengthofHomozygotRegion = -1;
		IntervalTreeNode intervalTreeNode = null;
		
		if (numberofConsecutiveHomozygotsFound>= numberofConsecutiveHomVariantsRequired){
			
			lengthofHomozygotRegion = lastSaved1BasedStart-saved1BasedStart+1;
			
			//Output the homozygot region
			bufferedWriter.write(chrName + "\t" + saved1BasedStart + "\t" + lastSaved1BasedStart  + "\t" + lengthofHomozygotRegion + System.getProperty("line.separator"));
			
			//Fill Interval Tree Node List
			intervalTreeNode = new IntervalTreeNode(chromosomeName,saved1BasedStart,lastSaved1BasedStart);	
			intervalTreeNodeList.add(intervalTreeNode);
		
		}
		
	}
	
	
	
	public static void findLROHs(
			Map<ChromosomeName, IntervalTree> case_chromosomeName2IntervalTreeMap,
			int numberofConsecutiveHomVariantsRequired,
			String sortedWESDataInputFileName,
			String case_outputFileName,
			String case_AfterOverlapsRemoved_outputFileName,
			String control_mother_outputFileName,
			String control_father_outputFileName){
		
		
		FileReader fileReader = null;
		BufferedReader bufferedReader = null;
		
		FileWriter case_fileWriter_LROHs = null;
		BufferedWriter case_bufferedWriter_LROHs = null;
		
		FileWriter case_AfterOverlapsRemoved_fileWriter_LROHs = null;
		BufferedWriter case_AfterOverlapsRemoved_bufferedWriter_LROHs = null;
	
		FileWriter control_mother_fileWriter_LROHs = null;
		BufferedWriter control_mother_bufferedWriter_LROHs = null;

		FileWriter control_father_fileWriter_LROHs = null;
		BufferedWriter control_father_bufferedWriter_LROHs = null;

		String strLine = null;
		
		String chrName = null;
		ChromosomeName chromosomeName = null;
		int _1BasedPosition = -1;
				
		int case_Saved1BasedStart = -1;
		int case_LastSaved1BasedStart = -1;
		int case_NumberofConsecutiveHomozygots = 0;
		
		int control_mother_Saved1BasedStart = -1;
		int control_mother_LastSaved1BasedStart = -1;
		int control_mother_NumberofConsecutiveHomozygots = 0;

		int control_father_Saved1BasedStart = -1;
		int control_father_LastSaved1BasedStart = -1;
		int control_father_NumberofConsecutiveHomozygots = 0;

		int indexofTab = -1;
		int indexofFormerTab = -1;
		int count = 0;
		
		String Case_13D0201103_mut = null;
		String Control_13D0201099_mut_father = null;
		String Control_13D0201100_mut_mother = null;	
		
		List<IntervalTreeNode> control_mother_IntervalTreeNodeList = new ArrayList<IntervalTreeNode>();
		List<IntervalTreeNode> control_father_IntervalTreeNodeList = new ArrayList<IntervalTreeNode>();
		
		//In order to do this input file must be sorted in ascending order for each chromosome				
		try {
			
			//Input
			fileReader = FileOperations.createFileReader(sortedWESDataInputFileName);
			bufferedReader = new BufferedReader(fileReader);
			
			//Outputs
			case_fileWriter_LROHs = FileOperations.createFileWriter(case_outputFileName);
			case_bufferedWriter_LROHs = new BufferedWriter(case_fileWriter_LROHs);
			
			case_AfterOverlapsRemoved_fileWriter_LROHs = FileOperations.createFileWriter(case_AfterOverlapsRemoved_outputFileName);
			case_AfterOverlapsRemoved_bufferedWriter_LROHs = new BufferedWriter(case_AfterOverlapsRemoved_fileWriter_LROHs);

			control_mother_fileWriter_LROHs = FileOperations.createFileWriter(control_mother_outputFileName);
			control_mother_bufferedWriter_LROHs = new BufferedWriter(control_mother_fileWriter_LROHs);
			
			control_father_fileWriter_LROHs = FileOperations.createFileWriter(control_father_outputFileName);
			control_father_bufferedWriter_LROHs = new BufferedWriter(control_father_fileWriter_LROHs);
			
			//Fill chromosomeBasedIntervalTrees for case and controls
								
			//Skip header line
			strLine = bufferedReader.readLine();
			
			//Write header line
			case_bufferedWriter_LROHs.write("chrName" + "\t" + "1BasedStart" + "\t" + "1BasedEnd" + "\t" + "lengthofHomozygotRegion" + "\t" + "numberofConsecutiveHomozygotsVariantsInCase" + System.getProperty("line.separator"));
			case_AfterOverlapsRemoved_bufferedWriter_LROHs.write("chrName" + "\t" + "1BasedStart" + "\t" + "1BasedEnd" + "\t" + "lengthofHomozygotRegion" + "\t" + "numberofConsecutiveHomozygotsVariantsInCase" + System.getProperty("line.separator"));
			control_mother_bufferedWriter_LROHs.write("chrName" + "\t" + "1BasedStart" + "\t" + "1BasedEnd" + "\t" + "lengthofHomozygotRegion" + "\t" + "numberofConsecutiveHomozygotsVariantsInCase" + System.getProperty("line.separator"));
			control_father_bufferedWriter_LROHs.write("chrName" + "\t" + "1BasedStart" + "\t" + "1BasedEnd" + "\t" + "lengthofHomozygotRegion" + "\t" + "numberofConsecutiveHomozygotsVariantsInCase" + System.getProperty("line.separator"));
		
			String formerChrName = "dummy";
			
			while((strLine = bufferedReader.readLine())!=null){
				
				//Chromosome	Position	Reference	GeneName	OMIM	Inheritance	Function	HGVS	ShareNumber	Control_13D0201099_mut	Ration	Control_13D0201100_mut	Ration	Case_13D0201103_mut	Ration	dbSNP_fre	1000human_fre	Hapmap_fre	Agilent_38M_fre	Agilent_46M_fre	Agilent_50M_fre	Nimblegen_44M_fre	Prediction from SIFT	Score from SIFT	RS-ID	NM-ID	Sub-region	Strand	ResidueChange	Gene description	GO_BP	GO_MF	GO_CC	KEGG_Pathway
				//Initialize				
				count = 0;
				
				indexofTab = strLine.indexOf('\t');
								
				while (indexofTab>0 && count < 14){
					
					count++;
					
					if (count==1){
						chrName = strLine.substring(0, indexofTab);	
						chromosomeName = ChromosomeName.convertStringtoEnum(chrName);
						
						if (!chrName.equalsIgnoreCase(formerChrName)){
							
							//Maybe we have not written the last homozygot region during chrName change							
							//Case starts
							isHomozygotRegionFound(
									case_NumberofConsecutiveHomozygots,
									numberofConsecutiveHomVariantsRequired,
									chrName,
									chromosomeName,
									case_Saved1BasedStart,
									case_LastSaved1BasedStart,
									case_bufferedWriter_LROHs,
									case_chromosomeName2IntervalTreeMap);							
							//Case ends
							
							//Control_mother starts
							isHomozygotRegionFound(
									control_mother_NumberofConsecutiveHomozygots,
									numberofConsecutiveHomVariantsRequired,
									chrName,
									chromosomeName,
									control_mother_Saved1BasedStart,
									control_mother_LastSaved1BasedStart,
									control_mother_bufferedWriter_LROHs,
									control_mother_IntervalTreeNodeList);

							//Control_mother ends
							
							
							//Control_father starts
							isHomozygotRegionFound(
									control_father_NumberofConsecutiveHomozygots,
									numberofConsecutiveHomVariantsRequired,
									chrName,
									chromosomeName,
									control_father_Saved1BasedStart,
									control_father_LastSaved1BasedStart,
									control_father_bufferedWriter_LROHs,
									control_father_IntervalTreeNodeList);
							//Control_father ends
														
							//We have started a new chromosome
							//Chromosome has changed
							//Then initialize
							case_NumberofConsecutiveHomozygots = 0;
							control_mother_NumberofConsecutiveHomozygots = 0;
							control_father_NumberofConsecutiveHomozygots = 0;
							formerChrName = chrName;
							
						}
					}else if (count==2){
						_1BasedPosition= Integer.parseInt(strLine.substring(indexofFormerTab+1, indexofTab));
					}else if(count==10){

						//CONTROL_FATHER
						Control_13D0201099_mut_father = strLine.substring(indexofFormerTab+1, indexofTab);
						
						if (Control_13D0201099_mut_father.endsWith("Hom")){
							if (control_father_NumberofConsecutiveHomozygots==0){
								//We have just started 
								control_father_Saved1BasedStart = _1BasedPosition;
							}
							control_father_NumberofConsecutiveHomozygots++;
							control_father_LastSaved1BasedStart = _1BasedPosition;
						}else if (Control_13D0201099_mut_father.endsWith("Het") || Control_13D0201099_mut_father.endsWith("ref")){
							
							isHomozygotRegionFound(
									control_father_NumberofConsecutiveHomozygots,
									numberofConsecutiveHomVariantsRequired,
									chrName,
									chromosomeName,
									control_father_Saved1BasedStart,
									control_father_LastSaved1BasedStart,
									control_father_bufferedWriter_LROHs,
									control_father_IntervalTreeNodeList);
							
							//Homozygot region is lost
							//Initialize
							control_father_NumberofConsecutiveHomozygots = 0;
						}
						
						
					}else if(count==12){
						
						//CONTROL_MOTHER
						Control_13D0201100_mut_mother = strLine.substring(indexofFormerTab+1, indexofTab);
						
						if (Control_13D0201100_mut_mother.endsWith("Hom")){
							if (control_mother_NumberofConsecutiveHomozygots==0){
								//We have just started 
								control_mother_Saved1BasedStart = _1BasedPosition;
							}
							control_mother_NumberofConsecutiveHomozygots++;
							control_mother_LastSaved1BasedStart = _1BasedPosition;
						}else if (Control_13D0201100_mut_mother.endsWith("Het") || Control_13D0201100_mut_mother.endsWith("ref")){
							
							isHomozygotRegionFound(
									control_mother_NumberofConsecutiveHomozygots,
									numberofConsecutiveHomVariantsRequired,
									chrName,
									chromosomeName,
									control_mother_Saved1BasedStart,
									control_mother_LastSaved1BasedStart,
									control_mother_bufferedWriter_LROHs,
									control_mother_IntervalTreeNodeList);
							
							//Homozygot region is lost
							//Initialize
							control_mother_NumberofConsecutiveHomozygots = 0;
						}
						
						
					}else if(count==14){
						
						//CASE PROBAND
						Case_13D0201103_mut = strLine.substring(indexofFormerTab+1, indexofTab);
						
						if (Case_13D0201103_mut.endsWith("Hom")){
							if (case_NumberofConsecutiveHomozygots==0){
								//We have just started 
								case_Saved1BasedStart = _1BasedPosition;
							}
							case_NumberofConsecutiveHomozygots++;
							case_LastSaved1BasedStart = _1BasedPosition;
						}else if (Case_13D0201103_mut.endsWith("Het") || Case_13D0201103_mut.endsWith("ref")){
							
							isHomozygotRegionFound(
									case_NumberofConsecutiveHomozygots,
									numberofConsecutiveHomVariantsRequired,
									chrName,
									chromosomeName,
									case_Saved1BasedStart,
									case_LastSaved1BasedStart,
									case_bufferedWriter_LROHs,
									case_chromosomeName2IntervalTreeMap);		
							
							//Homozygot region is lost
							//Initialize
							case_NumberofConsecutiveHomozygots = 0;
						}
						
						//Please notice that we don't care - cases for Case_13D0201103_mut
						//We skip them
						
					}
					
					indexofFormerTab = indexofTab;
					indexofTab = strLine.indexOf('\t',indexofTab+1);
					
				}//End of WHILE
								
							
			}//End of while reading input file
			
			
			
			removeOverlapsFromCaseUsingControlAutozygotRegions(
					case_chromosomeName2IntervalTreeMap,
					control_mother_IntervalTreeNodeList);
			
			removeOverlapsFromCaseUsingControlAutozygotRegions(
					case_chromosomeName2IntervalTreeMap,
					control_father_IntervalTreeNodeList);
			
			writeCaseIntervalTree(case_chromosomeName2IntervalTreeMap,case_AfterOverlapsRemoved_bufferedWriter_LROHs);
			
			//Close
			bufferedReader.close();
			case_bufferedWriter_LROHs.close();
			case_AfterOverlapsRemoved_bufferedWriter_LROHs.close();
			control_mother_bufferedWriter_LROHs.close();
			control_father_bufferedWriter_LROHs.close();

			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		

		
	}
	
	public static void augmentWithVariants(
			Map<ChromosomeName, IntervalTree> case_chromosomeName2IntervalTreeMap,
			String sortedWESDataInputFileName,
			String case_AfterOverlapsRemoved_AugmentedWithVariants_outputFileName,
			String case_AfterOverlapsRemoved_AugmentedWithNotCommonandSynonmousVariants_outputFileName,
			float selectionCriteriaForRareVariant){
		

		FileReader fileReader = null;
		BufferedReader bufferedReader = null;
		
		FileWriter case_fileWriter_AfterOverlapsRemoved_AugmentedWithVariants = null;
		BufferedWriter case_bufferedWriter_AfterOverlapsRemoved_AugmentedWithVariants = null;
		
		FileWriter case_fileWriter_AfterOverlapsRemoved_AugmentedWithNotCommonandSynonmousVariants = null;
		BufferedWriter case_bufferedWriter_AfterOverlapsRemoved_AugmentedWithNotCommonandSynonmousVariants = null;
		
		String strLine = null;
		
		String chrName = null;
		ChromosomeName chromosomeName = null;
		int _1BasedPosition = -1;
//		String reference = null;
//		String geneName=null;
//		String HGVS = null;
		
		//Chromosome	Position	Reference	GeneName	OMIM	Inheritance	Function	HGVS	ShareNumber	Control_13D0201099_mut	Ration	Control_13D0201100_mut	Ration	Case_13D0201103_mut	Ration	dbSNP_fre	1000human_fre	Hapmap_fre	Agilent_38M_fre	Agilent_46M_fre	Agilent_50M_fre	Nimblegen_44M_fre	Prediction from SIFT	Score from SIFT	RS-ID	NM-ID	Sub-region	Strand	ResidueChange	Gene description	GO_BP	GO_MF	GO_CC	KEGG_Pathway

		int indexofSixteenthTab;
		int indexofSeventeenthTab;
		int indexofEighteenthTab;
		int indexofNineteenthTab;
		int indexofTwenythTab;
		int indexofTwentyFirstTab;
		int indexofTwentySecondTab;
				
		int indexofTab = -1;
		int indexofFormerTab = -1;
		int count = 0;
		
		String function = null;
//		String Case_13D0201103_mut = null;
//		String Control_13D0201099_mut_father = null;
//		String Control_13D0201100_mut_mother = null;		
//		String rsID = null;
		
		float dbSNP_fre = 0f;
		float _1000human_fre = 0f;
		float Hapmap_fre = 0f;
		float Agilent_38M_fre = 0f;
		float Agilent_46M_fre = 0f;
		float Agilent_50M_fre = 0f;
		float Nimblegen_44M_fre = 0f;

		IntervalTreeNode intervalTreeNode = null;
		IntervalTree intervalTree = null;
		
		String LROH_chrName = null;
		int LROH_start = -1;
		int LROH_end = -1;
		int LROH_length = -1;
		
		
		try {
			
			//Input
			fileReader = FileOperations.createFileReader(sortedWESDataInputFileName);
			bufferedReader = new BufferedReader(fileReader);
			
			//Outputs
			case_fileWriter_AfterOverlapsRemoved_AugmentedWithVariants = FileOperations.createFileWriter(case_AfterOverlapsRemoved_AugmentedWithVariants_outputFileName);
			case_bufferedWriter_AfterOverlapsRemoved_AugmentedWithVariants = new BufferedWriter(case_fileWriter_AfterOverlapsRemoved_AugmentedWithVariants);
			
			case_fileWriter_AfterOverlapsRemoved_AugmentedWithNotCommonandSynonmousVariants = FileOperations.createFileWriter(case_AfterOverlapsRemoved_AugmentedWithNotCommonandSynonmousVariants_outputFileName);
			case_bufferedWriter_AfterOverlapsRemoved_AugmentedWithNotCommonandSynonmousVariants = new BufferedWriter(case_fileWriter_AfterOverlapsRemoved_AugmentedWithNotCommonandSynonmousVariants);
			
			
			//Skip header line
			strLine = bufferedReader.readLine();
			
			//Write header line
			case_bufferedWriter_AfterOverlapsRemoved_AugmentedWithVariants.write(strLine + "\t" + "LROH_chrName" + "\t" + "LROH_start" + "\t" + "LROH_end" + "\t" + "LROH_length" + System.getProperty("line.separator"));
			case_bufferedWriter_AfterOverlapsRemoved_AugmentedWithNotCommonandSynonmousVariants.write(strLine + "\t" + "LROH_chrName" + "\t" + "LROH_start" + "\t" + "LROH_end" + "\t" + "LROH_length" + System.getProperty("line.separator"));
	
			
			//We are reading Sorted WES data
			while((strLine = bufferedReader.readLine())!=null){
				
				//Initialize
				function = null;
				
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
					}else if (count==7){
						function = strLine.substring(indexofFormerTab+1, indexofTab);
					}
//					else if(count==10){
//						Control_13D0201099_mut_father = strLine.substring(indexofFormerTab+1, indexofTab);
//					}else if(count==12){
//						Control_13D0201100_mut_mother = strLine.substring(indexofFormerTab+1, indexofTab);
//					}else if(count==14){
//						Case_13D0201103_mut = strLine.substring(indexofFormerTab+1, indexofTab);
//					}
					
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
				
				
				chromosomeName = ChromosomeName.convertStringtoEnum(chrName);

				//create the interval
				intervalTreeNode = new IntervalTreeNode(chromosomeName,_1BasedPosition,_1BasedPosition);
				
				//get the intervalTree
				intervalTree = case_chromosomeName2IntervalTreeMap.get(chromosomeName);
				
				if(intervalTree!=null && intervalTree.getRoot()!=null & intervalTree.getRoot().getNodeName().isNotSentinel()){
					
					//Does it overlap with any interval in the case_afterOverlappingLROHsFromControlsRemoved_intervalTree?
					List<IntervalTreeNode> overlappedNodeList = new ArrayList<IntervalTreeNode>();
					intervalTree.findAllOverlappingIntervals(overlappedNodeList, intervalTree.getRoot(),intervalTreeNode);
					
					//If overlaps add it to the chrName2VariantList
					if (overlappedNodeList.size()>0){
						
						//debug delete later
						if (overlappedNodeList.size()>1){
							System.out.println("Can it be?");
						}
						//debug delete later
						
						LROH_chrName = ChromosomeName.convertEnumtoString(overlappedNodeList.get(0).getChromName());
						LROH_start = overlappedNodeList.get(0).getLow();
						LROH_end = overlappedNodeList.get(0).getHigh();
						LROH_length = LROH_end-LROH_start+1;
						
						case_bufferedWriter_AfterOverlapsRemoved_AugmentedWithVariants.write(strLine + "\t" + LROH_chrName + "\t" + LROH_start + "\t" +LROH_end + "\t" + LROH_length + System.getProperty("line.separator"));
						
						if ( 	(dbSNP_fre<selectionCriteriaForRareVariant) && (_1000human_fre<selectionCriteriaForRareVariant)  && (Hapmap_fre<selectionCriteriaForRareVariant) &&
								(Agilent_38M_fre<selectionCriteriaForRareVariant) && (Agilent_46M_fre<selectionCriteriaForRareVariant) && (Agilent_50M_fre<selectionCriteriaForRareVariant) &&
								(Nimblegen_44M_fre<selectionCriteriaForRareVariant)){
							
							//Filter Synonymous SNPs
							if (!function.startsWith("Synonymous")){
								
								case_bufferedWriter_AfterOverlapsRemoved_AugmentedWithNotCommonandSynonmousVariants.write(strLine +  "\t" + LROH_chrName + "\t" + LROH_start + "\t" + LROH_end + "\t" + LROH_length+ System.getProperty("line.separator"));
								
							}//End of IF NOT Synonymous
						
						}//End of IF RARE Variant
						
						
					}//End of if there is overlap
					
					overlappedNodeList = null;
					
				}//End of IF intervalTree is not null
				
				
				
			}//End of WHILE reading sorted WES data in ascending 1BasedSNP position
			
			
	
			
			//Close
			bufferedReader.close();
			case_bufferedWriter_AfterOverlapsRemoved_AugmentedWithVariants.close();
			case_bufferedWriter_AfterOverlapsRemoved_AugmentedWithNotCommonandSynonmousVariants.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}


	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		//Read WES Data containing Case and Controls
		//Input
		String sortedWESDataInputFileName = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\LROHs\\Sorted_LGMD-FamB-WES-All_chr_result.tep.txt";
		
		//Output
		String case_outputFileName = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\LROHs\\Case_LROHs_FROM_WES.txt";
		String case_AfterOverlapsRemoved_outputFileName = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\LROHs\\Case_LROHs_AfterOverlappingControlLROHsRemoved.txt";
		String control_mother_outputFileName = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\LROHs\\Control_mother_LROHs_FROM_WES.txt";
		String control_father_outputFileName = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\LROHs\\Control_father_LROHs_FROM_WES.txt";
		
		//Output Augmented with Variants
		String case_AfterOverlapsRemoved_AugmentedWithVariants_outputFileName = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\LROHs\\Case_LROHs_AfterOverlappingControlLROHsRemoved_AugmentedWithVariants.txt";
		String case_AfterOverlapsRemoved_AugmentedWithNotCommonandSynonmousVariants_outputFileName = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\LROHs\\Case_LROHs_AfterOverlappingControlLROHsRemoved_AugmentedWithNotCommonandSynonmousVariants.txt";
	
		
		int numberofConsecutiveHomVariantsRequired = 10;
		float selectionCriteriaForRareVariant = 0.01f;

		Map<ChromosomeName, IntervalTree> case_chromosomeName2IntervalTreeMap= new HashMap<ChromosomeName, IntervalTree>();
		
		//Find Long Runs of homozygosity (LROH) with at least 20 consecutive homozygot variants
		findLROHs(
				case_chromosomeName2IntervalTreeMap,
				numberofConsecutiveHomVariantsRequired,
				sortedWESDataInputFileName,
				case_outputFileName,
				case_AfterOverlapsRemoved_outputFileName,
				control_mother_outputFileName,
				control_father_outputFileName);
		
		augmentWithVariants(
				case_chromosomeName2IntervalTreeMap,
				sortedWESDataInputFileName,
				case_AfterOverlapsRemoved_AugmentedWithVariants_outputFileName,
				case_AfterOverlapsRemoved_AugmentedWithNotCommonandSynonmousVariants_outputFileName,
				selectionCriteriaForRareVariant);
		
	}

}
