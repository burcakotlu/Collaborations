/**
 * 
 */
package candidateVariantIdentification;

import intervaltree.IntervalTree;
import intervaltree.IntervalTreeNode;

import java.io.BufferedWriter;
import java.util.List;
import java.util.Map;

import enumtypes.ChromosomeName;

/**
 * @author Burçak Otlu
 * @date Jan 23, 2017
 * @project Collaborations 
 *
 */
public class Control {
	
	String name;
	
	int saved1BasedStart;
	int lastSaved1BasedStart;
	int numberofConsecutiveHomozygotsFound;
	
	BufferedWriter bufferedWriter;
	List<IntervalTreeNode> intervalTreeNodeList;
	Map<ChromosomeName, IntervalTree> chromosomeName2IntervalTreeMap;
	
	
	
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}
	
	public int getSaved1BasedStart() {
		return saved1BasedStart;
	}
	public void setSaved1BasedStart(int saved1BasedStart) {
		this.saved1BasedStart = saved1BasedStart;
	}
	public int getLastSaved1BasedStart() {
		return lastSaved1BasedStart;
	}
	public void setLastSaved1BasedStart(int lastSaved1BasedStart) {
		this.lastSaved1BasedStart = lastSaved1BasedStart;
	}
	public int getNumberofConsecutiveHomozygotsFound() {
		return numberofConsecutiveHomozygotsFound;
	}
	public void setNumberofConsecutiveHomozygotsFound(int numberofConsecutiveHomozygotsFound) {
		this.numberofConsecutiveHomozygotsFound = numberofConsecutiveHomozygotsFound;
	}
	public BufferedWriter getBufferedWriter() {
		return bufferedWriter;
	}
	public void setBufferedWriter(BufferedWriter bufferedWriter) {
		this.bufferedWriter = bufferedWriter;
	}
	public List<IntervalTreeNode> getIntervalTreeNodeList() {
		return intervalTreeNodeList;
	}
	public void setIntervalTreeNodeList(List<IntervalTreeNode> intervalTreeNodeList) {
		this.intervalTreeNodeList = intervalTreeNodeList;
	}
	public Map<ChromosomeName, IntervalTree> getChromosomeName2IntervalTreeMap() {
		return chromosomeName2IntervalTreeMap;
	}
	public void setChromosomeName2IntervalTreeMap(Map<ChromosomeName, IntervalTree> chromosomeName2IntervalTreeMap) {
		this.chromosomeName2IntervalTreeMap = chromosomeName2IntervalTreeMap;
	}
	
	

}
