/**
 * 
 */
package auxiliary;

import enumtypes.ChromosomeName;
import gnu.trove.list.TIntList;

/**
 * @author Burçak Otlu
 * @date Dec 1, 2016
 * @project Collaborations 
 *
 */
public class HG38RefSeqGeneInformation {
	
	String refSeqGeneName; 
	ChromosomeName chromName;
	char strand;
	int txStart;
	int txEnd;
	int exonCounts;
	String geneSymbol; 
	TIntList exonStartList;
	TIntList exonEndList;


	public String getRefSeqGeneName() {
		return refSeqGeneName;
	}





	public void setRefSeqGeneName(String refSeqGeneName) {
		this.refSeqGeneName = refSeqGeneName;
	}





	public ChromosomeName getChromName() {
		return chromName;
	}





	public void setChromName(ChromosomeName chromName) {
		this.chromName = chromName;
	}





	public char getStrand() {
		return strand;
	}





	public void setStrand(char strand) {
		this.strand = strand;
	}





	public int getTxStart() {
		return txStart;
	}





	public void setTxStart(int txStart) {
		this.txStart = txStart;
	}





	public int getTxEnd() {
		return txEnd;
	}





	public void setTxEnd(int txEnd) {
		this.txEnd = txEnd;
	}





	public int getExonCounts() {
		return exonCounts;
	}





	public void setExonCounts(int exonCounts) {
		this.exonCounts = exonCounts;
	}



	public String getGeneSymbol() {
		return geneSymbol;
	}





	public void setGeneSymbol(String geneSymbol) {
		this.geneSymbol = geneSymbol;
	}





	public TIntList getExonStartList() {
		return exonStartList;
	}





	public void setExonStartList(TIntList exonStartList) {
		this.exonStartList = exonStartList;
	}





	public TIntList getExonEndList() {
		return exonEndList;
	}





	public void setExonEndList(TIntList exonEndList) {
		this.exonEndList = exonEndList;
	}




	public HG38RefSeqGeneInformation(
			String refSeqGeneName, 
			ChromosomeName chromName, 
			char strand, 
			int txStart, 
			int txEnd, 
			int exonCounts, 
			String geneSymbol, 
			TIntList exonStartList, 
			TIntList exonEndList) {
		
		super();
		this.refSeqGeneName = refSeqGeneName;
		this.chromName = chromName;
		this.strand = strand;
		this.txStart = txStart;
		this.txEnd = txEnd;
		this.exonCounts = exonCounts;
		this.geneSymbol = geneSymbol;
		this.exonStartList = exonStartList;
		this.exonEndList = exonEndList;
	}





	public HG38RefSeqGeneInformation() {
		super();
	}


	

}
