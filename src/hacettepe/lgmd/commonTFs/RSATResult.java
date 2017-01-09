/**
 * 
 */
package hacettepe.lgmd.commonTFs;

import java.util.Comparator;

import auxiliary.FunctionalElementMinimal;

/**
 * @author Burçak Otlu
 * @date Jan 8, 2017
 * @project Collaborations 
 *
 */
public class RSATResult {
	
	String intervalName;
	String tfName;
	
	String matrixName;
	int matrixNumber;
	
	char direction;
	int start;
	int end;
	String sequence;
	Double pValue;
	
	
	
	
	public String getIntervalName() {
		return intervalName;
	}
	public void setIntervalName(String intervalName) {
		this.intervalName = intervalName;
	}
	public String getTfName() {
		return tfName;
	}
	public void setTfName(String tfName) {
		this.tfName = tfName;
	}
	public String getMatrixName() {
		return matrixName;
	}
	public void setMatrixName(String matrixName) {
		this.matrixName = matrixName;
	}
	public int getMatrixNumber() {
		return matrixNumber;
	}
	public void setMatrixNumber(int matrixNumber) {
		this.matrixNumber = matrixNumber;
	}
	public char getDirection() {
		return direction;
	}
	public void setDirection(char direction) {
		this.direction = direction;
	}
	public int getStart() {
		return start;
	}
	public void setStart(int start) {
		this.start = start;
	}
	public int getEnd() {
		return end;
	}
	public void setEnd(int end) {
		this.end = end;
	}
	public String getSequence() {
		return sequence;
	}
	public void setSequence(String sequence) {
		this.sequence = sequence;
	}
	
	
	
	
	public Double getpValue() {
		return pValue;
	}
	public void setpValue(Double pValue) {
		this.pValue = pValue;
	}
	public RSATResult(String intervalName, String tfName, String matrixName, int matrixNumber, char direction, int start, int end, String sequence, double pValue) {
		super();
		this.intervalName = intervalName;
		this.tfName = tfName;
		this.matrixName = matrixName;
		this.matrixNumber = matrixNumber;
		this.direction = direction;
		this.start = start;
		this.end = end;
		this.sequence = sequence;
		this.pValue = pValue;
	}
	
	public RSATResult() {
		super();
		// TODO Auto-generated constructor stub
	}
	
	
	
	public static Comparator<RSATResult> P_VALUE = new Comparator<RSATResult>() {

		public int compare( RSATResult element1, RSATResult element2) {

			return element1.getpValue().compareTo( element2.getpValue());

		}
	};
	

}
