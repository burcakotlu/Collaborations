/**
 * 
 */
package candidateVariantIdentification;

import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import intervaltree.IntervalTree;
import intervaltree.IntervalTreeNode;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import javax.swing.BorderFactory;
import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.border.LineBorder;

import auxiliary.FileOperations;
import enumtypes.ChromosomeName;

/**
 * @author Burçak Otlu
 * @date Dec 26, 2016
 * @project Collaborations 
 * 
 * This class finds candidate variants by discarding common LROHs between case and controls from WES data.
 *
 */
public class CandidateVariantIdentificationFromWES extends JPanel {
	
	private static JComboBox<String> caseColumn;
	private static JList<String> controlList;
	

	
	private static final long serialVersionUID = -6778100197471577791L;	
	
	public final static int ONE_MILLION = 1000000;
	
	static Map<ChromosomeName, IntervalTree> caseBefore_chromosomeName2IntervalTreeMap= new HashMap<ChromosomeName, IntervalTree>();
	static Map<ChromosomeName, IntervalTree> caseAfter_chromosomeName2IntervalTreeMap= new HashMap<ChromosomeName, IntervalTree>();
	
	static TIntObjectMap<Control> controlColumnNumber2ControlFeatureMap = new TIntObjectHashMap<Control>();
	    
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
	
	public static void writeCaseIntervalTree(
			Map<ChromosomeName,IntervalTree> case_chromosomeName2IntervalTreeMap,
			BufferedWriter case_bufferedWriter_LROHs_AfterOverlappingLROHsRemoved){
		
		IntervalTree intervalTree = null;
		
		for(Map.Entry<ChromosomeName, IntervalTree> entry: case_chromosomeName2IntervalTreeMap.entrySet()) {
			
			intervalTree = entry.getValue();
			
			if (intervalTree!=null && intervalTree.getRoot()!=null && intervalTree.getRoot().getNodeName().isNotSentinel()){
				IntervalTree.intervalTreeInfixTraversal(intervalTree.getRoot(),case_bufferedWriter_LROHs_AfterOverlappingLROHsRemoved);
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
						
						//Remove the overlap
						//There can be one or two remaining intervals left
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
			bufferedWriter.flush();
			
			//Construct Interval Tree
			intervalTreeNode = new IntervalTreeNode(chromosomeName,saved1BasedStart,lastSaved1BasedStart);								
			copyIntervalTreeNode = new IntervalTreeNode(chromosomeName,saved1BasedStart,lastSaved1BasedStart);								
			
			insert(chromosomeName,intervalTreeNode,before_chrName2IntervalTreeMap);
			insert(chromosomeName,copyIntervalTreeNode,after_chrName2IntervalTreeMap);
			
		}
		
	}
	
	//For control
	public static void isHomozygotRegionFound(
			int numberofConsecutiveHomVariantsRequired,			
			String chrName,
			ChromosomeName chromosomeName,
			Control control) throws IOException{
		
		int lengthofHomozygotRegion = -1;
		IntervalTreeNode intervalTreeNode = null;
		
		if (control.getNumberofConsecutiveHomozygotsFound()>= numberofConsecutiveHomVariantsRequired){
			
			
			lengthofHomozygotRegion = control.getLastSaved1BasedStart() - control.getSaved1BasedStart()+1;
			
			//Output the homozygot region
			control.getBufferedWriter().write(chrName + "\t" + control.getSaved1BasedStart() + "\t" + control.getLastSaved1BasedStart()  + "\t" + lengthofHomozygotRegion + System.getProperty("line.separator"));
			
			//Fill Interval Tree Node List
			intervalTreeNode = new IntervalTreeNode(chromosomeName,control.getSaved1BasedStart(),control.getLastSaved1BasedStart());	
			control.getIntervalTreeNodeList().add(intervalTreeNode);
			
			//Add to chrName2IntervalTreeMap
			insert(chromosomeName,intervalTreeNode,control.getChromosomeName2IntervalTreeMap());
		
		}//End of IF
		
	}
	
	
	public static void createInitializeAndWriteHeaderLines(
			List<String> controlColumnNames,
			int[] controlCounts,
			String outputDirectory,
			TIntObjectMap<Control> controlColumnNumber2ControlFeatureMap) throws IOException{
		
		int controlcolumnNumber = -1;
		String controlcolumnName = null;
		
		FileWriter control_fileWriter_LROHs = null;
		BufferedWriter control_bufferedWriter_LROHs = null;
		
		Control controlFeatures = null;

		
		for(int i=0; i<controlColumnNames.size(); i++){
			
			controlcolumnName = controlColumnNames.get(i);
			controlcolumnNumber = controlCounts[i];
			
			//Create fileWriter and bufferedWriter
			control_fileWriter_LROHs = FileOperations.createFileWriter(outputDirectory +  System.getProperty("file.separator") + controlcolumnName  + "_LROHs_FROM_WES.txt");
			control_bufferedWriter_LROHs = new BufferedWriter(control_fileWriter_LROHs);
			
			//Write header Line
			control_bufferedWriter_LROHs.write("chrName" + "\t" + "1BasedStart" + "\t" + "1BasedEnd" + "\t" + "lengthofHomozygotRegion" + "\t" + "numberofConsecutiveHomozygotsVariantsInCase" + System.getProperty("line.separator"));
			
			
			//Initialize
			controlFeatures = new Control();
			
			controlFeatures.setName(controlcolumnName);
			
			controlFeatures.setIntervalTreeNodeList(new ArrayList<IntervalTreeNode>());
			controlFeatures.setChromosomeName2IntervalTreeMap(new HashMap<ChromosomeName, IntervalTree>());
			controlFeatures.setBufferedWriter(control_bufferedWriter_LROHs);
			
			controlFeatures.setSaved1BasedStart(-1);
			controlFeatures.setLastSaved1BasedStart(-1);
			controlFeatures.setNumberofConsecutiveHomozygotsFound(0);			
			
			//Put into map
			controlColumnNumber2ControlFeatureMap.put(controlcolumnNumber, controlFeatures);
				
		}//End of FOR
		
	}
	
	
	public static int fillCountsforCaseandControls(
			String headerLine,
			String caseColumnName,			
			List<String> controlColumnNames,
			//TIntObjectMap<String>
			int[] controlColumnNumbers){
		
		int caseColumnNumber = -1;
		
		String columnName = null;
		
		int columnNumber = 0;
		
		int  indexofTab = -1;
		int indexofPreviousTab = -1;
		
		indexofTab = headerLine.indexOf('\t');
		
		while(indexofTab>=0){
			
				
			columnName = headerLine.substring(indexofPreviousTab+1, indexofTab); 
			columnNumber++;
			
			if (columnName.equals(caseColumnName)){
				caseColumnNumber = columnNumber;
			}//End of IF
			
			for(int i=0; i<controlColumnNames.size();i++) {
				if (columnName.equals(controlColumnNames.get(i))){
					controlColumnNumbers[i] = columnNumber;
				}//End of IF
			}//End of FOR
			
			//update indexes
			indexofPreviousTab = indexofTab;
			indexofTab = headerLine.indexOf('\t',indexofPreviousTab+1);
			
		}//End of while
		
		
		return caseColumnNumber;
		
	}

	
	public static boolean isCountInControlColumnNumbers(
			int columnNumber,
			int[] controlColumnNumbers){
		
		boolean contains = false;
		
		for(int i=0; i<controlColumnNumbers.length; i++){			
			if (controlColumnNumbers[i]==columnNumber){
				contains = true;
				break;
			}
		}//End of for 
		
		return contains;
		
	}
	
	public static void findLROHs(
			Map<ChromosomeName, IntervalTree> caseBefore_chromosomeName2IntervalTreeMap,
			Map<ChromosomeName, IntervalTree> caseAfter_chromosomeName2IntervalTreeMap,
			String inputFileName,
			String caseColumnName,
			List<String> controlColumnNames,
			int numberofConsecutiveHomVariantsRequired,
			String outputDirectory){
		
		
		FileReader fileReader = null;
		BufferedReader bufferedReader = null;
		
		FileWriter case_fileWriter_LROHs_before = null;
		BufferedWriter case_bufferedWriter_LROHs_before = null;
		
		FileWriter case_fileWriter_LROHs_AfterOverlappingLROHsRemoved = null;
		BufferedWriter case_bufferedWriter_LROHs_AfterOverlappingLROHsRemoved = null;
	
		String strLine = null;
		
		String chrName = null;
		ChromosomeName chromosomeName = null;
		int _1BasedPosition = -1;
				
		int case_Saved1BasedStart = -1;
		int case_LastSaved1BasedStart = -1;
		int case_NumberofConsecutiveHomozygots = 0;
		
		Control control = null;
				
		int indexofTab = -1;
		int indexofFormerTab = -1;
		int count = 0;
		
		String Case_13D0201103_mut = null;
		String Control_Column = null;
		
		int caseColumnNumber = -1;
		int[] controlColumnNumbers = new int[controlColumnNames.size()];
		
		
		//In order to do this input file must be sorted in ascending order for each chromosome				
		try {
			
			//Input
			fileReader = FileOperations.createFileReader(inputFileName);
			bufferedReader = new BufferedReader(fileReader);
			
			//Outputs for Case Before
			case_fileWriter_LROHs_before = FileOperations.createFileWriter(outputDirectory + System.getProperty("file.separator") + caseColumnName + "_LROHs_Before.txt");
			case_bufferedWriter_LROHs_before = new BufferedWriter(case_fileWriter_LROHs_before);
			
			//Outputs for Case After
			case_fileWriter_LROHs_AfterOverlappingLROHsRemoved = FileOperations.createFileWriter(outputDirectory +  System.getProperty("file.separator") + caseColumnName + "_LROHs_AfterOverlappingLROHsRemoved.txt");
			case_bufferedWriter_LROHs_AfterOverlappingLROHsRemoved = new BufferedWriter(case_fileWriter_LROHs_AfterOverlappingLROHsRemoved);
											
			//Skip header line
			strLine = bufferedReader.readLine();
			
			//Get column numbers for case and controls
			caseColumnNumber = fillCountsforCaseandControls(strLine,caseColumnName,controlColumnNames,controlColumnNumbers);
			
			//Does outputDirectory has file separator at the end? No.
			createInitializeAndWriteHeaderLines(
					controlColumnNames,
					controlColumnNumbers,
					outputDirectory,
					controlColumnNumber2ControlFeatureMap);
			
			//For case write header lines
			case_bufferedWriter_LROHs_before.write("chrName" + "\t" + "1BasedStart" + "\t" + "1BasedEnd" + "\t" + "lengthofHomozygotRegion" + "\t" + "numberofConsecutiveHomozygotsVariantsInCase" + System.getProperty("line.separator"));
			case_bufferedWriter_LROHs_AfterOverlappingLROHsRemoved.write("chrName" + "\t" + "1BasedStart" + "\t" + "1BasedEnd" + "\t" + "lengthofHomozygotRegion" + "\t" + "numberofConsecutiveHomozygotsVariantsInCase" + System.getProperty("line.separator"));
		
			String formerChrName = "dummy";
			
			while((strLine = bufferedReader.readLine())!=null){
				
				//Chromosome	Position	Reference	GeneName	OMIM	Inheritance	Function	HGVS	ShareNumber	Control_13D0201099_mut	Ration	Control_13D0201100_mut	Ration	Case_13D0201103_mut	Ration	dbSNP_fre	1000human_fre	Hapmap_fre	Agilent_38M_fre	Agilent_46M_fre	Agilent_50M_fre	Nimblegen_44M_fre	Prediction from SIFT	Score from SIFT	RS-ID	NM-ID	Sub-region	Strand	ResidueChange	Gene description	GO_BP	GO_MF	GO_CC	KEGG_Pathway
				//Initialize				
				count = 0;
				
				indexofTab = strLine.indexOf('\t');
								
				while (indexofTab>0 && count < 14){
					
					count++;
					
					//Same as in vcf format
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
									case_bufferedWriter_LROHs_before,
									caseBefore_chromosomeName2IntervalTreeMap,
									caseAfter_chromosomeName2IntervalTreeMap);							
							//Case ends
							
							//We have started a new chromosome
							//Chromosome has changed
							//Then initialize
							case_NumberofConsecutiveHomozygots = 0;
							
							//For each Control starts
							for (TIntObjectIterator<Control> itr = controlColumnNumber2ControlFeatureMap.iterator();itr.hasNext(); ){
								
								itr.advance();
								
								control = itr.value();
								
								isHomozygotRegionFound(									
										numberofConsecutiveHomVariantsRequired,
										chrName,
										chromosomeName,
										control);
								
								//We have started a new chromosome
								//Chromosome has changed
								//Then initialize
								control.setNumberofConsecutiveHomozygotsFound(0);
								
							}//For each Control ends
														
							formerChrName = chrName;
							
						}
					}//End of count==1 chromosomeColumn
					
					//Same as in vcf format
					else if (count==2){
						_1BasedPosition= Integer.parseInt(strLine.substring(indexofFormerTab+1, indexofTab));
					}
					
					else if(isCountInControlColumnNumbers(count,controlColumnNumbers)){
						
						//Get the control
						control = controlColumnNumber2ControlFeatureMap.get(count);
						
						//Get the control column value
						Control_Column = strLine.substring(indexofFormerTab+1, indexofTab);
						
						if (Control_Column.endsWith("Hom")){
							if (control.getNumberofConsecutiveHomozygotsFound()==0){
								//We have just started 
								control.setSaved1BasedStart(_1BasedPosition);
							}
							control.setNumberofConsecutiveHomozygotsFound(control.getNumberofConsecutiveHomozygotsFound()+1);
							control.setLastSaved1BasedStart(_1BasedPosition);
						}else if (Control_Column.endsWith("Het") || Control_Column.endsWith("ref")){
							
							isHomozygotRegionFound(
									numberofConsecutiveHomVariantsRequired,
									chrName,
									chromosomeName,
									control);
							
							//Homozygot region is lost
							//Initialize
							control.setNumberofConsecutiveHomozygotsFound(0);
						}
						
						
					}//End of IF one of control column
					
					
					else if(count==caseColumnNumber){
						
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
									case_bufferedWriter_LROHs_before,
									caseBefore_chromosomeName2IntervalTreeMap,
									caseAfter_chromosomeName2IntervalTreeMap);		
							
							//Homozygot region is lost
							//Initialize
							case_NumberofConsecutiveHomozygots = 0;
						}
							
					}//End of CASE
					
					indexofFormerTab = indexofTab;
					indexofTab = strLine.indexOf('\t',indexofTab+1);
					
				}//End of WHILE
								
							
			}//End of while reading input file
			
		  
			//Remove overlapping LROHs from case for each control
			for(TIntObjectIterator<Control> itr = controlColumnNumber2ControlFeatureMap.iterator();itr.hasNext();){
				
				itr.advance();			
				control = itr.value();
				
				removeOverlapsFromCaseUsingControlAutozygotRegions(
						caseAfter_chromosomeName2IntervalTreeMap,
						control.getIntervalTreeNodeList());
			}//End of for each control
			
			writeCaseIntervalTree(caseAfter_chromosomeName2IntervalTreeMap,case_bufferedWriter_LROHs_AfterOverlappingLROHsRemoved);
			
		
			//Close BufferedReader and BufferedWriters
			bufferedReader.close();
			case_bufferedWriter_LROHs_before.close();
			case_bufferedWriter_LROHs_AfterOverlappingLROHsRemoved.close();
			
			//Close bufferedWriter for each control
			for (TIntObjectIterator<Control> itr = controlColumnNumber2ControlFeatureMap.iterator();itr.hasNext(); ){
				
				itr.advance();
		
				control = itr.value();
				control.getBufferedWriter().close();
			}
		

			
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
		
		FileWriter case_Hom_or_Het_filter_wrt_aotosomal_recessive_mode_fileWriter = null;
		BufferedWriter case_Hom_or_Het_filter_wrt_aotosomal_recessive_mode_bufferedWriter = null;
		
		FileWriter case_Hom_filter_wrt_aotosomal_recessive_mode_fileWriter = null;
		BufferedWriter case_Hom_filter_wrt_aotosomal_recessive_mode_bufferedWriter = null;

		FileWriter case_fileWriter_AfterOverlapsRemoved_AugmentedWithVariants = null;
		BufferedWriter case_bufferedWriter_AfterOverlapsRemoved_AugmentedWithVariants = null;
		
		FileWriter case_fileWriter_AfterOverlapsRemoved_AugmentedWithNotCommonandSynonmousVariants = null;
		BufferedWriter case_bufferedWriter_AfterOverlapsRemoved_AugmentedWithNotCommonandSynonmousVariants = null;
		
		String strLine = null;
		
		String chrName = null;
		ChromosomeName chromosomeName = null;
		int _1BasedPosition = -1;
		
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
			case_Hom_or_Het_filter_wrt_aotosomal_recessive_mode_fileWriter =  FileOperations.createFileWriter("C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\LROHs\\excel_case_Hom_or_Het_filter.txt");
			case_Hom_or_Het_filter_wrt_aotosomal_recessive_mode_bufferedWriter = new BufferedWriter(case_Hom_or_Het_filter_wrt_aotosomal_recessive_mode_fileWriter);

			case_Hom_filter_wrt_aotosomal_recessive_mode_fileWriter =  FileOperations.createFileWriter("C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\LROHs\\excel_case_Hom_filter.txt");
			case_Hom_filter_wrt_aotosomal_recessive_mode_bufferedWriter = new BufferedWriter(case_Hom_filter_wrt_aotosomal_recessive_mode_fileWriter);

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
						//TODO Can control variant be  "ref"? 
						if ( (Case_13D0201103_mut.contains("Hom") || Case_13D0201103_mut.contains("Het")) && Control_13D0201099_mut_father.contains("Het") && Control_13D0201100_mut_mother.contains("Het")){
							
							//debug for excel filtering versions starts
							case_Hom_or_Het_filter_wrt_aotosomal_recessive_mode_bufferedWriter.write(strLine + System.getProperty("line.separator"));								
							if (Case_13D0201103_mut.contains("Hom")){
								case_Hom_filter_wrt_aotosomal_recessive_mode_bufferedWriter.write(strLine + System.getProperty("line.separator"));																
							}
							//debug for excel filtering versions ends
							
							
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
									
									//TODO for numberofconsecutive homozygot variants required for LROH = 1 debug it
									//It must give the same excel filter output case = Hom and control = Het
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
									
																			
									//This is done in  Autosomal Recessive Model: add check whether Case is Hom and Control is Het at candidate Loci
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
			
			//Close for debugging purposes
			case_Hom_or_Het_filter_wrt_aotosomal_recessive_mode_bufferedWriter.close();
			case_Hom_filter_wrt_aotosomal_recessive_mode_bufferedWriter.close();
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	
	
	
	
	public void paintComponent(Graphics g) {
		
		Control control = null;
		//Integer controlColumnNumber= new Integer(0);
		
		char[] chrCharArray = null;
		
		char[] caseCharArray = null;
		char[] beforeCharArray = null;
		char[] afterCharArray = null;
		char[] controlCharArray = null;
	
		///Enlarge the y position
		int enlargeFactor = 3;
		
		Map<ChromosomeName,Position> caseAfter_chromosomeName2RectangleTopLeftPositionsMap = new HashMap<ChromosomeName,Position> ();
		
		//Shows the top left most point where the drawing starts
		int top_left_x = 30;
		int top_left_y = 60;
		
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
	        
      
			for(TIntObjectIterator<Control> itr = controlColumnNumber2ControlFeatureMap.iterator();itr.hasNext();){
				
				itr.advance();
			
				//controlColumnNumber = itr.key();
				control = itr.value();

				//TODO write the control Column Name not the column number
				//But it is long beware of this
				drawLROHs(g,
	        			chrName,
	        			control.getChromosomeName2IntervalTreeMap(),
	        			control.getName().toCharArray(),
	        			controlCharArray,
	        			chrCharArray,
	        			hg19ChromosomeName2ChromSize,
	        			top_left_x,
	        			top_left_y,
	        			enlargeFactor,
	        			width);  
	        	
	          	//Update top_left_x position for the next chromosome rectangle
	        	top_left_x = top_left_x + width + widthBetweenRectangles;
				
			}//End of for each control
			

        	
 
        	
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
        
        showCandidateVariants(
        		g,
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
    
    public static List<String> readColumnNamesFromHeaderLine(String inputFileName){
    	
    	FileReader fileReader = null;
    	BufferedReader bufferedReader = null;
    	
    	String strLine = null;
    	int indexofTab;
    	int indexofPreviousTab =-1;
    	
    	String columnName;
    	List<String> columnNames = new ArrayList<String>();
    	
		try {
			fileReader = new FileReader(inputFileName);
			bufferedReader = new  BufferedReader(fileReader);
			
			if((strLine = bufferedReader.readLine())!=null){
				
				indexofTab = strLine.indexOf('\t',indexofPreviousTab+1);
				
				while(indexofTab>0){
					columnName = strLine.substring(indexofPreviousTab+1, indexofTab);
					columnNames.add(columnName);
					
					indexofPreviousTab = indexofTab;
					indexofTab = strLine.indexOf('\t',indexofPreviousTab+1);

				}//End of while
				
				//For the last one
				if (indexofPreviousTab>0){
					columnName = strLine.substring(indexofPreviousTab+1);
					columnNames.add(columnName);
					
				}
				
			}//Read First header line
			
			//Close
			bufferedReader.close();
			
			
	    	
		} catch (FileNotFoundException e) {			
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
    	
		return columnNames;
    }
        
	public static void main(String[] args) {
		
		//TODO Simplify the code
		//TODO make GUI parametric
		
		//Input File: Read WES Data containing Case and Controls
		//It has to be sorted w.r.t. chromoseme and start end points in ascending order
		//String sortedWESDataInputFileName = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\LROHs\\Sorted_LGMD-FamB-WES-All_chr_result.tep.txt";
		
		//Output Augmented with Variants
		String case_AfterOverlapsRemoved_AugmentedWithVariants_outputFileName = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\LROHs\\Case_LROHs_AfterOverlappingControlLROHsRemoved_AugmentedWithVariants.txt";
		
		//After Meeting
		//Do not augment with control Hom snps
		String case_AfterOverlapsRemoved_AugmentedWithNotCommonandSynonmousVariants_outputFileName = "C:\\Users\\Burçak\\Google Drive\\Collaborations\\HacettepeUniversity\\LGMD\\LROHs\\Case_LROHs_AfterOverlappingControlLROHsRemoved_AugmentedWithNotCommonandSynonmousVariants.txt";
	
		//GUI starts		
		GridBagConstraints constraints = new GridBagConstraints();
        constraints.anchor = GridBagConstraints.PAGE_START;
        constraints.insets = new Insets(5, 5, 5, 5);
        constraints.fill = GridBagConstraints.BOTH;

		JPanel userPanel = new JPanel();
		userPanel.setLayout(new GridBagLayout());
		//userPanel.setPreferredSize(new Dimension(7000,100));
		userPanel.setMaximumSize(new Dimension(7000,200));
		//userPanel.setMinimumSize(new Dimension(7000,100));
		
		JPanel casePanel =  new JPanel();
		casePanel.setLayout(new GridBagLayout());
	
		JPanel controlPanel =  new JPanel();
		controlPanel.setLayout(new GridBagLayout());


        /******************************************/
        /**********Input File Panel starts*********/
        /******************************************/
		JPanel inputFilePanel = new JPanel();
		inputFilePanel.setLayout(new GridBagLayout());

		String[] options = {};
		caseColumn = new JComboBox<String>(options);
		controlList = new JList<String>(options);

		
		JTextField inputFileTextField = new JTextField(30);
		
		JButton browseButton = new JButton("Browse");
		browseButton.addActionListener(new ActionListener() {
            
			@Override
            public void actionPerformed(ActionEvent e) {
            	
            	JFileChooser fc = new JFileChooser();
            	int returnVal;

            	fc.setFileSelectionMode(JFileChooser.FILES_ONLY);
            	
            	returnVal = fc.showOpenDialog(userPanel);
            	if( returnVal == JFileChooser.APPROVE_OPTION){
            		
            		File file = fc.getSelectedFile();
            		inputFileTextField.setText(file.getPath() + System.getProperty( "file.separator"));
            		
            		//TODO sort coordinates based on start position in ascending order chromosome based
                	List<String> columns = readColumnNamesFromHeaderLine(inputFileTextField.getText());
            		
            		
            		//Populate case columns
                	caseColumn.removeAllItems();
            		for(int i=0;i<columns.size();i++){
            			caseColumn.addItem(columns.get(i));                		
            		}
            		
                  	//Populate control columns
            		controlList.removeAll();
            		DefaultListModel<String> model = new DefaultListModel<String>();
            		//Autoscroll did not work as I expected.
            		controlList.setAutoscrolls(true);
            		controlList.setBorder(new LineBorder(Color.BLACK));
            		for(int i=0;i<columns.size();i++){
            			model.addElement(columns.get(i));            			
            		}
            		controlList.setModel(model);
					
            	}

            }
        });

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
		JLabel caseColumnLabel = new JLabel("Select Case Column:");
		JLabel controlColumnLabel = new JLabel("Select Control Column/s:");

		caseColumn.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {
                // Do something when you select a value
            }
        });

			
		constraints.gridx = 0;
	    constraints.gridy = 0; 
		casePanel.add(caseColumnLabel,constraints);
		
		constraints.gridx = 1;
	    constraints.gridy = 0; 
		casePanel.add(caseColumn,constraints);
		
		constraints.gridx = 0;
	    constraints.gridy = 0; 
		controlPanel.add(controlColumnLabel,constraints);
		
		constraints.gridx = 1;
	    constraints.gridy = 0; 
		controlPanel.add(controlList,constraints);
		
		//Set border for the panel
		casePanel.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(), "Select Case Column"));
		
		//Set border for the panel
		controlPanel.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(), "Select Control Column/s"));
        /******************************************/
        /****Case and control Panel ends***********/
        /******************************************/
		
	    /******************************************/
        /**********Run Panel starts****************/
        /******************************************/
		JPanel runPanel =  new JPanel();
		runPanel.setLayout(new GridBagLayout());
		
		// Button submit
        JButton runButton = new JButton("Show Candidate Variants...");
        runButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
            	
            	int numberofConsecutiveHomVariantsRequired = Integer.parseInt(numberOfConsecutiveHomVariantsRequired.getText());
            	float selectionCriteriaForRareVariant = Float.parseFloat(rareVariantThreshold.getText());
            	
            	String inputFilename = inputFileTextField.getText();
            	String outputDirectory = null;
            			
            	//TODO Minor user can provide output Folder
            	//Right now we get it from inputfile path
            	File file = new File(inputFilename);
            	if (file.isFile()){
            		outputDirectory = file.getParent();              
            	}
            	
            	String caseColumnName = (String) caseColumn.getSelectedItem();
            	
            	int[] controlColumnIndicies = controlList.getSelectedIndices();
            	String controlColumnName = null;
            	List<String> controlColumnNames = new ArrayList<String>();
            	            	
            	for (int i = 0; i < controlColumnIndicies.length; i++) {
            		controlColumnName = controlList.getModel().getElementAt(controlColumnIndicies[i]);
            		controlColumnNames.add(controlColumnName);
            	 }
            	
            	
            	//TODO call the methods below if variables are set correctly            	
            	//This program is more useful in case of control data is provided
            	//On the other hand, it can still information only case data is provided.
            	//Case data is mandatory
            	if(caseColumnName!=null && !caseColumnName.isEmpty() && controlColumnNames!=null){
   
            		//Find Long Runs of homozygosity (LROH) with at least --numberofConsecutiveHomVariantsRequired-- consecutive homozygot variants
            		findLROHs(
            				caseBefore_chromosomeName2IntervalTreeMap,
            				caseAfter_chromosomeName2IntervalTreeMap,				
            				inputFilename,
            				caseColumnName,
            				controlColumnNames,
            				numberofConsecutiveHomVariantsRequired,
            				outputDirectory);
            		
            		//TODO Can we do this part using VCF file?
            		//TODO has to be parametric too.
            		augmentWithVariants(
            				caseAfter_chromosomeName2IntervalTreeMap,
            				inputFilename,
            				case_AfterOverlapsRemoved_AugmentedWithVariants_outputFileName,
            				case_AfterOverlapsRemoved_AugmentedWithNotCommonandSynonmousVariants_outputFileName,
            				selectionCriteriaForRareVariant);

            	}//End of if user provided variables are entered properly 
            	
          		
        		CandidateVariantIdentificationFromWES mainPanel = new CandidateVariantIdentificationFromWES();
        		//TODO 7000 must be parametric
        		//decide on number of case and controls and calculate the real value instead of 7000
        		mainPanel.setPreferredSize(new Dimension(7000, 1000));
        		mainPanel.setLayout(new BorderLayout());
        		
        	    JScrollPane scrollPane = new JScrollPane(mainPanel, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);  //Let all scrollPanel has scroll bars
        	    //TODO make it parametric
        	    scrollPane.setPreferredSize(new Dimension(1000, 900));

        	    JFrame frame2 = new JFrame("Show Candidate Variants by discarding common LROHs between case and controls from WES");
        		frame2.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        		
        		frame2.add(scrollPane);
        		//TODO make it parametric
        	    frame2.setSize(1000, 900);
        		
        		frame2.pack();
        		frame2.setLocationRelativeTo(null);
        		frame2.setVisible(true); 

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
        constraints.gridheight  =1;
        userPanel.add(casePanel, constraints);
        
		constraints.gridx = 1;
        constraints.gridy = 1; 
        constraints.gridheight  =1;
        userPanel.add(controlPanel, constraints);
 
        
        constraints.gridx = 0;
        constraints.gridy = 2; 
        constraints.gridwidth = 2;
		userPanel.add(runPanel, constraints);
		
		// set border for the panel
		userPanel.setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(), "User Panel"));
	    /******************************************/
        /******Add Panels toUser Panel ends********/
        /******************************************/
       
	    JFrame frame = new JFrame("Candidate Variant Identification by discarding common LROHs between case and controls from WES");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		//I have added two panels two my frame.
		//How to give upper panel a shorter height and lower panel a higher height?
		//frame.add(userPanel);		
		//frame.add(scrollPane);
		//frame.setSize(1000, 900);
		
		//What does it do?
		//How graphics draws to mainPanel? How is it decided? Since mainPanel is the current instance that extends JPanel
		//Therefore graphics are drawn on mainPanel
		//frame.setLayout(new GridLayout(0, 1));;
		frame.setLayout(new GridBagLayout());
		
		//I have added two panels two my frame.
		//How to give upper panel a shorter height and lower panel a higher height?
		
		constraints.gridx = 0;
        constraints.gridy = 0;     
        constraints.weighty=1;
       frame.add(userPanel,constraints);		
		
       constraints.gridx = 0;
       constraints.gridy = 1;
       constraints.gridheight=4;
       constraints.weighty=10;
       constraints.fill= GridBagConstraints.VERTICAL;     
       //frame.add(scrollPane,constraints);
       //frame.setSize(1000, 900);
		
		frame.pack();
		frame.setLocationRelativeTo(null);
		frame.setVisible(true); 
		//GUI ends	
		
	}

}
