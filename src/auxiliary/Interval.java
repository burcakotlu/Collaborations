/**
 * 
 */
package auxiliary;

import enumtypes.ChromosomeName;

/**
 * @author Burçak Otlu
 * @date Jan 6, 2017
 * @project Collaborations 
 *
 */
public class Interval {
	
	ChromosomeName chromosomeName;
	int low;
	int high;
	
	
	public ChromosomeName getChromosomeName() {
		return chromosomeName;
	}
	public void setChromosomeName(ChromosomeName chromosomeName) {
		this.chromosomeName = chromosomeName;
	}
	public int getLow() {
		return low;
	}
	public void setLow(int low) {
		this.low = low;
	}
	public int getHigh() {
		return high;
	}
	public void setHigh(int high) {
		this.high = high;
	}
	
	
	public Interval(ChromosomeName chromosomeName, int low, int high) {
		super();
		this.chromosomeName = chromosomeName;
		this.low = low;
		this.high = high;
	}
	
	

}
