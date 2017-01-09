/**
 * 
 */
package hacettepe.lgmd.commonTFs;

import java.util.Comparator;

/**
 * @author Burçak Otlu
 * @date Jan 9, 2017
 * @project Collaborations 
 *
 */
public class TFBasedResult {
	
	String tfName;
	Double avgPValue;
	Double accumulatedLogRatio;
	
	public String getTfName() {
		return tfName;
	}


	public void setTfName(String tfName) {
		this.tfName = tfName;
	}


	public Double getAvgPValue() {
		return avgPValue;
	}


	public void setAvgPValue(Double avgPValue) {
		this.avgPValue = avgPValue;
	}


	public Double getAccumulatedLogRatio() {
		return accumulatedLogRatio;
	}


	public void setAccumulatedLogRatio(Double accumulatedLogRatio) {
		this.accumulatedLogRatio = accumulatedLogRatio;
	}


	public TFBasedResult(String tfName, Double avgPValue, Double accumulatedLogRatio) {
		super();
		this.tfName = tfName;
		this.avgPValue = avgPValue;
		this.accumulatedLogRatio = accumulatedLogRatio;
	}
	
	
	public static Comparator<TFBasedResult> AVG_P_VALUE = new Comparator<TFBasedResult>() {

		public int compare( TFBasedResult element1, TFBasedResult element2) {

			return element1.getAvgPValue().compareTo( element2.getAvgPValue());

		}
	};
	

}
