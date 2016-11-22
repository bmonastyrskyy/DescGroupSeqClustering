package edu.ucdavis.gc.bm.descGroupSeqMetric;

import java.util.Set;

public interface GroupSeqSignal {
	/**
	 * 
	 * @return score of the entire group
	 */
	Double getSeqSignal();
	/**
	 * 
	 * @param cluster - set of indexes of descriptors in group
	 * @return score of the cluster
	 */
	Double getSeqSignal(Set<Integer> cluster);
	
}
