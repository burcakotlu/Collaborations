/**
 * 
 */
package hacettepe.lgmd.lrohs;

import intervaltree.IntervalTree;
import intervaltree.IntervalTreeNode;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
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

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextField;

import auxiliary.FileOperations;
import enumtypes.ChromosomeName;

/**
 * @author Burçak Otlu
 * @date Dec 26, 2016
 * @project Collaborations 
 *
 */
public class LROHsFromWESData extends JPanel {
	
	
	private static final long serialVersionUID = -6778100197471577791L;	
	
	public final static int ONE_MILLION = 1000000;
	
	static Map<ChromosomeName, IntervalTree> caseBefore_chromosomeName2IntervalTreeMap= new HashMap<ChromosomeName, IntervalTree>();
	static Map<ChromosomeName, IntervalTree> caseAfter_chromosomeName2IntervalTreeMap= new HashMap<ChromosomeName, IntervalTree>();
	
	static Map<ChromosomeName, IntervalTree> controlFather_chromosomeName2IntervalTreeMap= new HashMap<ChromosomeName, IntervalTree>();
	static Map<ChromosomeName, IntervalTree> controlMother_chromosomeName2IntervalTreeMap= new HashMap<ChromosomeName, IntervalTree>();
	    
	static List<IntervalTreeNode> candidateLocis = new ArrayList<IntervalTreeNode>();
	
	static Map<ChromosomeName,Integer> hg19ChromosomeName2ChromSize= new HashMap<ChromosomeName,Integer>();  
	
    
	
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
	
	
	//For case
	public static void isHomozygotRegionFound(
			int numberofConsecutiveHomozygotsFound,
			int numberofConsecutiveHomVariantsRequired,			
			String chrName,
			ChromosomeName chromosomeName,
			int saved1BasedStart,
			int lastSaved1BasedStart,
			BufferedWriter bufferedWriter,
			Map<ChromosomeName,IntervalTree> before_chrName2IntervalTreeMap,
			Map<ChromosomeName,IntervalTree> after_chrName2IntervalTreeMap) throws IOException{
		
		int lengthofHomozygotRegion = -1;
		IntervalTreeNode intervalTreeNode = null;
		IntervalTreeNode copyIntervalTreeNode = null;
		
		if (numberofConsecutiveHomozygotsFound>= numberofConsecutiveHomVariantsRequired){
			
			lengthofHomozygotRegion = lastSaved1BasedStart-saved1BasedStart+1;
			
			//Output the homozygot region
			bufferedWriter.write(chrName + "\t" + saved1BasedStart + "\t" + lastSaved1BasedStart  + "\t" + lengthofHomozygotRegion + System.getProperty("line.separator"));
			
			//Construct Interval Tree
			intervalTreeNode = new IntervalTreeNode(chromosomeName,saved1BasedStart,lastSaved1BasedStart);								
			copyIntervalTreeNode = new IntervalTreeNode(chromosomeName,saved1BasedStart,lastSaved1BasedStart);								
			
			insert(chromosomeName,intervalTreeNode,before_chrName2IntervalTreeMap);
			insert(chromosomeName,copyIntervalTreeNode,after_chrName2IntervalTreeMap);
			
		}
		
	}
	
	//For control
	public static void isHomozygotRegionFound(
			int numberofConsecutiveHomozygotsFound,
			int numberofConsecutiveHomVariantsRequired,			
			String chrName,
			ChromosomeName chromosomeName,
			int saved1BasedStart,
			int lastSaved1BasedStart,
			BufferedWriter bufferedWriter,
			List<IntervalTreeNode> intervalTreeNodeList,
			Map<ChromosomeName, IntervalTree> chromosomeName2IntervalTreeMap) throws IOException{
		
		int lengthofHomozygotRegion = -1;
		IntervalTreeNode intervalTreeNode = null;
		
		if (numberofConsecutiveHomozygotsFound>= numberofConsecutiveHomVariantsRequired){
			
			lengthofHomozygotRegion = lastSaved1BasedStart-saved1BasedStart+1;
			
			//Output the homozygot region
			bufferedWriter.write(chrName + "\t" + saved1BasedStart + "\t" + lastSaved1BasedStart  + "\t" + lengthofHomozygotRegion + System.getProperty("line.separator"));
			
			//Fill Interval Tree Node List
			intervalTreeNode = new IntervalTreeNode(chromosomeName,saved1BasedStart,lastSaved1BasedStart);	
			intervalTreeNodeList.add(intervalTreeNode);
			
			//Add to chrName2IntervalTreeMap
			insert(chromosomeName,intervalTreeNode,chromosomeName2IntervalTreeMap);
		
		}
		
	}
	
	
	
	public static void findLROHs(
			Map<ChromosomeName, IntervalTree> caseBefore_chromosomeName2IntervalTreeMap,
			Map<ChromosomeName, IntervalTree> caseAfter_chromosomeName2IntervalTreeMap,
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
									caseBefore_chromosomeName2IntervalTreeMap,
									caseAfter_chromosomeName2IntervalTreeMap);							
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
									control_mother_IntervalTreeNodeList,
									controlMother_chromosomeName2IntervalTreeMap);

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
									control_father_IntervalTreeNodeList,
									controlFather_chromosomeName2IntervalTreeMap);
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
									control_father_IntervalTreeNodeList,
									controlFather_chromosomeName2IntervalTreeMap);
							
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
									control_mother_IntervalTreeNodeList,
									controlMother_chromosomeName2IntervalTreeMap);
							
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
									caseBefore_chromosomeName2IntervalTreeMap,
									caseAfter_chromosomeName2IntervalTreeMap);		
							
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
					caseAfter_chromosomeName2IntervalTreeMap,
					control_mother_IntervalTreeNodeList);
			
			removeOverlapsFromCaseUsingControlAutozygotRegions(
					caseAfter_chromosomeName2IntervalTreeMap,
					control_father_IntervalTreeNodeList);
			
			writeCaseIntervalTree(caseAfter_chromosomeName2IntervalTreeMap,case_AfterOverlapsRemoved_bufferedWriter_LROHs);
			
		
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
			Map<ChromosomeName, IntervalTree> caseAfter_chromosomeName2IntervalTreeMap,
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
		
		String Control_13D0201099_mut_father = null;
		String Control_13D0201100_mut_mother = null;
		String Case_13D0201103_mut = null;
	
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
					else if(count==10){
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
				
				
				if ( 	(dbSNP_fre<selectionCriteriaForRareVariant) && (_1000human_fre<selectionCriteriaForRareVariant)  && (Hapmap_fre<selectionCriteriaForRareVariant) &&
						(Agilent_38M_fre<selectionCriteriaForRareVariant) && (Agilent_46M_fre<selectionCriteriaForRareVariant) && (Agilent_50M_fre<selectionCriteriaForRareVariant) &&
						(Nimblegen_44M_fre<selectionCriteriaForRareVariant)){
					
					//Filter Synonymous SNPs
					if (!function.startsWith("Synonymous")){
						
						//Autosomal Recessive Model
						if ( (Case_13D0201103_mut.contains("Hom") || Case_13D0201103_mut.contains("Het")) && Control_13D0201099_mut_father.contains("Het") && Control_13D0201100_mut_mother.contains("Het")){
							
							chromosomeName = ChromosomeName.convertStringtoEnum(chrName);
	
							//create the interval from SNP position
							intervalTreeNode = new IntervalTreeNode(chromosomeName,_1BasedPosition,_1BasedPosition);
							
							//get the intervalTree
							intervalTree = caseAfter_chromosomeName2IntervalTreeMap.get(chromosomeName);
							
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
									
																			
									//TODO add check whether Case is Hom and Control is Het at candidate Loci
									case_bufferedWriter_AfterOverlapsRemoved_AugmentedWithNotCommonandSynonmousVariants.write(strLine +  "\t" + LROH_chrName + "\t" + LROH_start + "\t" + LROH_end + "\t" + LROH_length+ System.getProperty("line.separator"));								
									candidateLocis.add(intervalTreeNode);
									
								}//End of if there is overlap
								
								overlappedNodeList = null;
								
							}//End of IF intervalTree is not null
						
						}//End of IF Autosomal Recessive Model
						
						
					}//End of IF NOT Synonymous
					
				}//End of IF RARE Variant
				
				
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
	
	
	
	
	
	public void paintComponent(Graphics g) {
		
		
		
		
		char[] chrCharArray = null;
		
		char[] caseCharArray = null;
		char[] beforeCharArray = null;
		char[] afterCharArray = null;
		char[] controlCharArray = null;
		char[] controlMotherCharArray = null;
		char[] controlFatherCharArray = null;
		
		///Enlarge the y position
		int enlargeFactor = 3;
		
		Map<ChromosomeName,Position> caseAfter_chromosomeName2RectangleTopLeftPositionsMap = new HashMap<ChromosomeName,Position> ();
				
		int top_left_x = 30;
		int top_left_y = 150;
		
		int width = 30;
		int widthBetweenRectangles = 30;
		int widthBetweenRectangleAndBorderLine = 20;
		
		int widthForBorderLine = 2;
		int heigthForBorderLine = 900/enlargeFactor;
		
		
		super.paintComponent(g);
	       
		fillHg19ChromosomeSizes();


        for(ChromosomeName chrName :ChromosomeName.values()){
        	
        
        	caseCharArray = "Case".toCharArray();
        	beforeCharArray = "Before".toCharArray();
        	afterCharArray = "After".toCharArray();
        	controlCharArray = "Control".toCharArray();        	
        	controlMotherCharArray = "Mother".toCharArray();
        	controlFatherCharArray = "Father".toCharArray();
        	chrCharArray = chrName.convertEnumtoString().toCharArray();

        	/*************************************************/
        	/*************CASE Before starts******************/
        	/*************************************************/
        	drawLROHs(g,
        			chrName,
        			caseBefore_chromosomeName2IntervalTreeMap,
        			beforeCharArray,
        			caseCharArray,
        			chrCharArray,
        			hg19ChromosomeName2ChromSize,
        			top_left_x,
        			top_left_y,
        			enlargeFactor,
        			width);   
        	
	        //Update top_left_x position for the next chromosome rectangle
        	top_left_x = top_left_x + width + widthBetweenRectangles;

        	/*************************************************/
        	/*************CASE Before ends********************/
        	/*************************************************/
        	
        	
           	/*************************************************/
        	/*************CASE After starts********************/
        	/*************************************************/
        	drawLROHs(g,
        			chrName,
        			caseAfter_chromosomeName2IntervalTreeMap,
        			afterCharArray,
        			caseCharArray,
        			chrCharArray,
        			hg19ChromosomeName2ChromSize,
        			top_left_x,
        			top_left_y,
        			enlargeFactor,
        			width);   
        	
        	caseAfter_chromosomeName2RectangleTopLeftPositionsMap.put(chrName, new Position(top_left_x, top_left_y));
        	
	        //Update top_left_x position for the next chromosome rectangle
        	top_left_x = top_left_x + width + widthBetweenRectangles;
        	/*************************************************/
        	/*************CASE After ends*********************/
        	/*************************************************/
	        
      
       	
        	/*************************************************/
        	/*************Control Mother starts***************/
        	/*************************************************/
        	drawLROHs(g,
        			chrName,
        			controlMother_chromosomeName2IntervalTreeMap,
        			controlMotherCharArray,
        			controlCharArray,
        			chrCharArray,
        			hg19ChromosomeName2ChromSize,
        			top_left_x,
        			top_left_y,
        			enlargeFactor,
        			width);  
        	
          	//Update top_left_x position for the next chromosome rectangle
        	top_left_x = top_left_x + width + widthBetweenRectangles;
        	/*************************************************/
        	/*************Control Mother ends*****************/
        	/*************************************************/
        	
         	
        	
           	/*************************************************/
        	/*************Control Father starts***************/
        	/*************************************************/
        	drawLROHs(g,
        			chrName,
        			controlFather_chromosomeName2IntervalTreeMap,        			
        			controlFatherCharArray,
        			controlCharArray,
        			chrCharArray,
        			hg19ChromosomeName2ChromSize,
        			top_left_x,
        			top_left_y,
        			enlargeFactor,
        			width);     
        	
        	//Update top_left_x position for the next chromosome rectangle
        	top_left_x = top_left_x + width + widthBetweenRectangleAndBorderLine;
        	/*************************************************/
        	/*************Control Father ends*****************/
        	/*************************************************/
        	
        	/*************************************************/
        	/****drawBorderLine Between chromsomes starts*****/
        	/*************************************************/
        	g.setColor(Color.black);
        	g.fillRect(top_left_x, top_left_y, widthForBorderLine, heigthForBorderLine*enlargeFactor);
        	
        	top_left_x = top_left_x + width + widthBetweenRectangleAndBorderLine;
        	/*************************************************/
        	/*****drawBorderLine Between chromsomes ends******/
        	/*************************************************/
        	
        	
        }//End of FOR
        
        showCandidateVariants(g,
        		candidateLocis,
        		width,
        		enlargeFactor,
        		caseAfter_chromosomeName2RectangleTopLeftPositionsMap);
 

    }
	
	public static void fillHg19ChromosomeSizes(){
		hg19ChromosomeName2ChromSize.put(ChromosomeName.CHROMOSOME1,249250621);
		hg19ChromosomeName2ChromSize.put(ChromosomeName.CHROMOSOME2,243199373);
		hg19ChromosomeName2ChromSize.put(ChromosomeName.CHROMOSOME3,198022430);
		hg19ChromosomeName2ChromSize.put(ChromosomeName.CHROMOSOME4,191154276);		
		hg19ChromosomeName2ChromSize.put(ChromosomeName.CHROMOSOME5,180915260);
		hg19ChromosomeName2ChromSize.put(ChromosomeName.CHROMOSOME6,171115067);
		hg19ChromosomeName2ChromSize.put(ChromosomeName.CHROMOSOME7,159138663);
		hg19ChromosomeName2ChromSize.put(ChromosomeName.CHROMOSOME8,146364022);
		hg19ChromosomeName2ChromSize.put(ChromosomeName.CHROMOSOME9,141213431);
		hg19ChromosomeName2ChromSize.put(ChromosomeName.CHROMOSOME10,135534747);
		hg19ChromosomeName2ChromSize.put(ChromosomeName.CHROMOSOME11,135006516);
		hg19ChromosomeName2ChromSize.put(ChromosomeName.CHROMOSOME12,133851895);
		hg19ChromosomeName2ChromSize.put(ChromosomeName.CHROMOSOME13,115169878);
		hg19ChromosomeName2ChromSize.put(ChromosomeName.CHROMOSOME14,107349540);
		hg19ChromosomeName2ChromSize.put(ChromosomeName.CHROMOSOME15,102531392);
		hg19ChromosomeName2ChromSize.put(ChromosomeName.CHROMOSOME16,90354753);
		hg19ChromosomeName2ChromSize.put(ChromosomeName.CHROMOSOME17,81195210);
		hg19ChromosomeName2ChromSize.put(ChromosomeName.CHROMOSOME18,78077248);
		hg19ChromosomeName2ChromSize.put(ChromosomeName.CHROMOSOME19,59128983);
		hg19ChromosomeName2ChromSize.put(ChromosomeName.CHROMOSOME20,63025520);
		hg19ChromosomeName2ChromSize.put(ChromosomeName.CHROMOSOME21,48129895);
		hg19ChromosomeName2ChromSize.put(ChromosomeName.CHROMOSOME22,51304566);
		hg19ChromosomeName2ChromSize.put(ChromosomeName.CHROMOSOMEX,155270560);
		hg19ChromosomeName2ChromSize.put(ChromosomeName.CHROMOSOMEY,59373566);
		
//		chr1	249250621
//		chr2	243199373
//		chr3	198022430
//		chr4	191154276
//		chr5	180915260
//		chr6	171115067
//		chr7	159138663
//		chrX	155270560
//		chr8	146364022
//		chr9	141213431
//		chr10	135534747
//		chr11	135006516
//		chr12	133851895
//		chr13	115169878
//		chr14	107349540
//		chr15	102531392
//		chr16	90354753
//		chr17	81195210
//		chr18	78077248
//		chr20	63025520
//		chrY	59373566
//		chr19	59128983
//		chr22	51304566
//		chr21	48129895

	}
	
	 public void drawLROHs(
	    		Graphics g,
				ChromosomeName chrName,
				Map<ChromosomeName,IntervalTree> chromosomeName2IntervalTreeMap,
				char[] beforeorAfterCharArray,
				char[] caseorControlCharArray,
				char[] chrCharArray,
				Map<ChromosomeName,Integer> hg19ChromosomeName2ChromSize,
				int top_left_x,
				int top_left_y,
				int enlargeFactor,
				int width){
	    	
	    	IntervalTree intervalTree = chromosomeName2IntervalTreeMap.get(chrName);
	    	int chromSize = -1;
	    	int height = -1;
	    	
	    	g.setColor(Color.BLACK);  
	    	
	      	//Write Case or Control at upper y position
	    	g.drawChars(caseorControlCharArray, 0, caseorControlCharArray.length, top_left_x, top_left_y-30);
	    	
	    	//Write Before or After at upper y position
	    	if (beforeorAfterCharArray!=null){
	        	g.drawChars(beforeorAfterCharArray, 0, beforeorAfterCharArray.length, top_left_x, top_left_y-20);    		
	    	}

	     	//Write chrName at lower y position
		     g.drawChars(chrCharArray, 0, chrCharArray.length, top_left_x, top_left_y-10);
		       	
			chromSize = hg19ChromosomeName2ChromSize.get(chrName);
			height =getConvertedNumber(chromSize);
			
			//Draw chromosome whole rectangle
			g.setColor(Color.GRAY);
	        g.fillRect(top_left_x,top_left_y,width,height*enlargeFactor);
	        
	        //Draw position indicator lines for each chromosome
	        g.setColor(Color.black);
	        for(Integer i=0; i< height;i= i+10){
	        	 g.drawChars(i.toString().toCharArray(), 0, i.toString().toCharArray().length, top_left_x-20,top_left_y+i*enlargeFactor);	        	
	        	 g.drawLine(top_left_x-5,top_left_y+i*enlargeFactor,top_left_x, top_left_y+i*enlargeFactor);	  	       
	        }
	        	      
	        //Draw LROHs of chromosome
	        if(intervalTree!= null && intervalTree.getRoot()!=null){
		        fillLROHs(intervalTree.getRoot(),g,top_left_x,top_left_y,width,enlargeFactor);
	        }
	    	
	    }
	    
	 public int getConvertedNumber(int x){
			return Math.round(x*1.0f/ONE_MILLION);
		}
		
	    public void fillLROHs(
	    		IntervalTreeNode node,
	    		Graphics g,
	    		int top_left_x,
	    		int top_left_y,
	    		int width,
	    		int enlargeFactor){
	    	
	    	int low;
	    	int high;
	    	int height;
	    	
	    	//Left Node
			if( node.getLeft()!= null && node.getLeft().getNodeName().isNotSentinel())
				fillLROHs(node.getLeft(),g,top_left_x,top_left_y,width,enlargeFactor);

			//Middle Node
			if( node.getNodeName().isNotSentinel()){
				
				//One way
				low =getConvertedNumber(node.getLow()); 
				high = getConvertedNumber(node.getHigh());
				
				height = high-low+1;
				
				//Fill rectangle for the LROH
				g.setColor(Color.RED);
				
				//Draw if LROH is greater than 10000 bp
				if (node.getHigh()-node.getLow()>10000){
					g.fillRect(top_left_x, top_left_y + (low*enlargeFactor), width, height*enlargeFactor);
				}

				
//				//Second Way didn't work
//				height = getConvertedNumber(node.getHigh()-node.getLow()+1);
//				
//				//Fill rectangle for the LROH
//				g.setColor(Color.RED);
//				g.fillRect(top_left_x, top_left_y + (low*enlargeFactor), width, height*enlargeFactor);

			}

			//Right Node
			if( node.getRight()!= null && node.getRight().getNodeName().isNotSentinel())
				fillLROHs(node.getRight(),g,top_left_x,top_left_y,width,enlargeFactor);

	    }
	    
	    public void showCandidateVariants(
	    		Graphics g,
	    		List<IntervalTreeNode> candidateLocis,
	    		int width,
	    		int enlargeFactor,
	    		Map<ChromosomeName,Position> caseAfter_chromosomeName2RectangleTopLeftPositionsMap){
	    	
//	    	String star = "*";
//	    	char[] starCharArray = star.toCharArray();
	    	
	    	int candidate_y = -1;
	    	ChromosomeName chromosomeName = null; 
	    	
	    	
	    	Position rectangleTopLeftPosition = null;
	    	int x = -1;
	    	int y = -1;
	    	
	    	for(IntervalTreeNode candidate:candidateLocis){
	    		
	    		candidate_y = getConvertedNumber(candidate.getHigh());
	    		
	    		chromosomeName = candidate.getChromName();
	    		
	    	    		
	    		rectangleTopLeftPosition = caseAfter_chromosomeName2RectangleTopLeftPositionsMap.get(chromosomeName);
	    	
	    		x = Math.round(rectangleTopLeftPosition.getX() + (width*1.0f/2));
	    		y = rectangleTopLeftPosition.getY()+ candidate_y*enlargeFactor;
	    		
	    		g.setColor(Color.CYAN);

	    		//Star did not draw at the good position
	    		//g.drawChars(starCharArray, 0, starCharArray.length, x, y);
	    		
	    		//x-4 is for centering it
	    		g.drawRect(x-4, y, 8, 1);
	    		
	    		g.drawRect(x, y-4, 1, 8);
	    
	    	}//End of for each candidate loci
	    	
	    }
	    

	//Implement this method
    public static void findLROHS(
    		String caseColumn, 
    		int[] controlSelectedIndices,
    		String inputFileTextField,
    		String numberOfConsecutiveHomVariantsRequired,
    		String rareVariantThreshold){
    	
    	
    }
        
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
		
		//After Meeting
		//Do not augment with control Hom snps
		String case_AfterOverlapsRemoved_AugmentedWithNotCommonandSynonmousVariants_outputFileName = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\LROHs\\Case_LROHs_AfterOverlappingControlLROHsRemoved_AugmentedWithNotCommonandSynonmousVariants.txt";
	
		int numberofConsecutiveHomVariantsRequired = 10;
		float selectionCriteriaForRareVariant = 0.01f;

		//Find Long Runs of homozygosity (LROH) with at least --numberofConsecutiveHomVariantsRequired-- consecutive homozygot variants
		findLROHs(
				caseBefore_chromosomeName2IntervalTreeMap,
				caseAfter_chromosomeName2IntervalTreeMap,				
				numberofConsecutiveHomVariantsRequired,
				sortedWESDataInputFileName,
				case_outputFileName,
				case_AfterOverlapsRemoved_outputFileName,
				control_mother_outputFileName,
				control_father_outputFileName);
		
	
		augmentWithVariants(
				caseAfter_chromosomeName2IntervalTreeMap,
				sortedWESDataInputFileName,
				case_AfterOverlapsRemoved_AugmentedWithVariants_outputFileName,
				case_AfterOverlapsRemoved_AugmentedWithNotCommonandSynonmousVariants_outputFileName,
				selectionCriteriaForRareVariant);
		
		
		
		//GUI starts
		LROHsFromWESData mainPanel = new LROHsFromWESData();
		//TODO 7000 must be parametric
		mainPanel.setPreferredSize(new Dimension(7000, 1000));
		mainPanel.setLayout(new BorderLayout());
		
		
		GridBagConstraints constraints = new GridBagConstraints();
        constraints.anchor = GridBagConstraints.PAGE_START;
        constraints.insets = new Insets(5, 5, 5, 5);
        constraints.fill = GridBagConstraints.BOTH;

		JPanel userPanel = new JPanel();
		userPanel.setLayout(new GridBagLayout());
		
	
		         
        /******************************************/
        /**********Input File Panel starts*********/
        /******************************************/
		JPanel inputFilePanel = new JPanel();
		inputFilePanel.setLayout(new GridBagLayout());

		JButton browseButton = new JButton("Browse");
		browseButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
            	//Do Something
            }
        });
		
		JTextField inputFileTextField = new JTextField(30);

		constraints.gridx = 0;
        constraints.gridy = 0;     
        inputFilePanel.add(browseButton, constraints);
        
        constraints.gridx = 1;
        constraints.gridy = 0; 
        inputFilePanel.add(inputFileTextField, constraints);
        
		// set border for the panel
		inputFilePanel.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(), "Load Input File"));
	    /******************************************/
        /**********Input File Panel ends***********/
        /******************************************/
		
        /******************************************/
        /**********Threshold Panel starts**********/
        /******************************************/
		JPanel thresholdPanel = new JPanel();
		thresholdPanel.setLayout(new GridBagLayout());
				
		JLabel numberOfConsecutiveHomVariantsRequiredLabel = new JLabel("Number Of Consecutive Homozygot Variants Required for LROH:");
		JTextField numberOfConsecutiveHomVariantsRequired = new JTextField(10);
		numberOfConsecutiveHomVariantsRequired.setText("10");
		
		JLabel rareVariantThresholdLabel = new JLabel("Rare Variant Threshold:");
		JTextField rareVariantThreshold = new JTextField(10);
		rareVariantThreshold.setText("0.05");
	
		constraints.gridx = 0;
	    constraints.gridy = 0; 
	    thresholdPanel.add(numberOfConsecutiveHomVariantsRequiredLabel,constraints);
		
		constraints.gridx = 1;
	    constraints.gridy = 0; 
	    thresholdPanel.add(numberOfConsecutiveHomVariantsRequired,constraints);
		
		constraints.gridx = 0;
	    constraints.gridy = 1; 
	    thresholdPanel.add(rareVariantThresholdLabel,constraints);
		
		constraints.gridx = 1;
	    constraints.gridy = 1; 
	    thresholdPanel.add(rareVariantThreshold,constraints);
		
		// set border for the panel
	    thresholdPanel.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(), "Thresholds"));
        /******************************************/
        /**********Threshold Panel ends************/
        /******************************************/
		
        /******************************************/
        /****Case and control Panel starts*********/
        /******************************************/
		JPanel caseControlPanel =  new JPanel();
		caseControlPanel.setLayout(new GridBagLayout());
		
		String[] options = { "Option1", "Option2", "Option3", "Option4", "Option15" };
		JLabel caseColumnLabel = new JLabel("Select Case Column:");
		JComboBox<String> caseColumn = new JComboBox<String>(options);
		caseColumn.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {
                // Do something when you select a value

            }
        });

		JLabel controlColumnLabel = new JLabel("Select Control Column/s:");
		JList<String> controlList = new JList<String>(options);

				
		constraints.gridx = 0;
	    constraints.gridy = 0; 
		caseControlPanel.add(caseColumnLabel,constraints);
		
		constraints.gridx = 1;
	    constraints.gridy = 0; 
		caseControlPanel.add(caseColumn,constraints);
		
		constraints.gridx = 0;
	    constraints.gridy = 1; 
		caseControlPanel.add(controlColumnLabel,constraints);
		
		constraints.gridx = 1;
	    constraints.gridy = 1; 
		caseControlPanel.add(controlList,constraints);
		
		// set border for the panel
		caseControlPanel.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(), "Select Case Column and Control Column/s"));
        /******************************************/
        /****Case and control Panel ends***********/
        /******************************************/
		
	    /******************************************/
        /**********Run Panel starts****************/
        /******************************************/
		JPanel runPanel =  new JPanel();
		runPanel.setLayout(new GridBagLayout());
		
		// Button submit
        JButton runButton = new JButton("Show LROHs...");
        runButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                findLROHS((String) caseColumn.getSelectedItem(), controlList.getSelectedIndices(),inputFileTextField.getText(),numberOfConsecutiveHomVariantsRequired.getText(),rareVariantThreshold.getText());
            }
        });
        
		constraints.gridx = 0;
	    constraints.gridy = 0; 
	    runPanel.add(runButton,constraints);		
	    /******************************************/
        /**********Run Panel ends******************/
        /******************************************/

		
	    /******************************************/
        /******Add Panels toUser Panel starts******/
        /******************************************/
		constraints.gridx = 0;
        constraints.gridy = 0; 
		userPanel.add(inputFilePanel, constraints);
		
		constraints.gridx = 0;
        constraints.gridy = 1; 
		userPanel.add(thresholdPanel, constraints);
		
		constraints.gridx = 1;
        constraints.gridy = 0; 
       userPanel.add(caseControlPanel, constraints);
        
        constraints.gridx = 1;
        constraints.gridy = 1; 
		userPanel.add(runPanel, constraints);
		
		// set border for the panel
		userPanel.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(), "User Panel"));
	    /******************************************/
        /******Add Panels toUser Panel ends********/
        /******************************************/
        
		
	    JScrollPane scrollPane = new JScrollPane(mainPanel, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);  //Let all scrollPanel has scroll bars
	    scrollPane.setPreferredSize(new Dimension(1000, 900));
		
	  
	    JFrame frame = new JFrame("LROHs from WES");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		//I have added two panels two my frame.
		//How to give upper panel a shorter height and lower panel a higher height?
		frame.add(userPanel);		
		frame.add(scrollPane);
		frame.setSize(1000, 900);
		
		//What does it do?
		//How graphics draws to mainPanel? How is it decided?
		frame.setLayout(new GridLayout(0, 1));;
		
		frame.pack();
		frame.setLocationRelativeTo(null);
		frame.setVisible(true); 
		//GUI ends
		
		
		
	}

}
