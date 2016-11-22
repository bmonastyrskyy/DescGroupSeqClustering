package edu.ucdavis.gc.bm.descGroupSeqMetric;

import java.util.Set;

import edu.ucdavis.gc.bm.descriptorGroup.Descriptor;
import edu.ucdavis.gc.bm.descriptorGroup.Group;
import edu.ucdavis.gc.pdb.PdbSubstitutionMatrix;

public class GroupSeqSignalNaiv implements GroupSeqSignal {

	private double[][] matrDistance;

	private Group group;

	private double seqSignal = 0;

	// private List<Double> seqSignalPerSegms = new ArrayList<Double>();

	public GroupSeqSignalNaiv() {

	}

	public GroupSeqSignalNaiv(Group group) {
		this.group = group;
		this.calcSeqSignal();
		this.calcMatrDist();
	}

	private void calcMatrDist() {
		this.matrDistance = new double[group.getNumberMembers()][group
				.getNumberMembers()];
		for (int i = 1; i < this.matrDistance.length; i++) {
			for (int j = 0; j < i; j++) {
				// '-' is chosed in order to return to "true" sequence score based on Blosum62
				// in class DescriptorSeqMetric the distanceTo() method used '-' in order to simulate the 
				// behavior of the distance measure (less - better)
				this.matrDistance[i][j] = -((DescriptorSeqMetric) group
						.getDescriptors().get(i))
						.distanceTo((DescriptorSeqMetric) group
								.getDescriptors().get(j));
				this.matrDistance[j][i] = this.matrDistance[i][j];
			}
		}
	}

	/**
	 * calculates signal sequence in very naive way: just apply Blosum62
	 * substitution matrix to pairs of amino acids at corresponding positions in
	 * descriptors in the group
	 */
	private void calcSeqSignal() {
		// initialize List of Segments scores

		int count = 0;
		for (int k = 0; k < group.getNumberSegments(); k++) { // loop over
																// segments
			int count_segm = 0;
			double seqSignalperSegm = 0;
			for (int i = 1; i < group.getNumberMembers(); i++) { // loop over
																	// first
																	// descriptor
				Descriptor di = group.getDescriptors().get(i);
				for (int j = 0; j < i; j++) { // loop over the second descriptor
					Descriptor dj = group.getDescriptors().get(j);
					String segmi = di.getSeqs().get(k).toUpperCase();
					String segmj = dj.getSeqs().get(k).toUpperCase();

					for (int m = 0; m < segmi.length(); m++) {
						count++;
						count_segm++;
						char ai = segmi.charAt(m);
						if (ai == '.' || ai == '-') {
							ai = '*';
						}
						char aj = segmj.charAt(m);
						if (aj == '.' || aj == '-') {
							aj = '*';
						}
						try {
							seqSignal += PdbSubstitutionMatrix.getScore(ai, aj);
							seqSignalperSegm += PdbSubstitutionMatrix.getScore(
									ai, aj);
						} catch (NullPointerException e) {
							System.out.println("ai " + ai + "   aj " + aj);
							throw e;
						}
					}
				}
			}
			// seqSignalPerSegms.add(seqSignalperSegm/count_segm);
		} // end loop over segments
		seqSignal /= count;
	}

	/**
	 * calculates signal sequence in very naive way: just apply Blosum62
	 * substitution matrix to pairs of amino acids at corresponding positions in
	 * descriptors in the group
	 */
	private Double calcSeqSignal(Set<Integer> cluster) {
		// initialize List of Segments scores
		Double result = 0.0;
		int count = 0;
		for (int i = 1; i < group.getNumberMembers(); i++) { // loop over first
																// descriptor
			if (!cluster.contains(i)) {
				continue;
			}
			for (int j = 0; j < i; j++) { // loop over the second descriptor
				if (!cluster.contains(j)) {
					continue;
				}
				result += this.matrDistance[i][j];
				count++;
			}
		}
		// seqSignalPerSegms.add(seqSignalperSegm/count_segm);
		// end loop over segments
		result /= count;
		return result;
	}

	@Override
	public Double getSeqSignal() {
		// TODO Auto-generated method stub
		return this.seqSignal;
	}

	@Override
	public Double getSeqSignal(Set<Integer> cluster) {
		// TODO Auto-generated method stub
		return this.calcSeqSignal(cluster);
	}

	/*
	 * TODO implement the code which calculates with seqSignal per segment
	 * public List<Double> getSeqSignalPerSegm(){ return this.seqSignalPerSegms;
	 * }
	 * 
	 * public String toStringSeqSignalPerSegm(){ String result = ""; for(double
	 * d: this.seqSignalPerSegms){ result += d + "\t"; } return result; }
	 */
}
