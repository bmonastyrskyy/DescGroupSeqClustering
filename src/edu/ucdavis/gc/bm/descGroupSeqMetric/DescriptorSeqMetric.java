package edu.ucdavis.gc.bm.descGroupSeqMetric;

import java.util.List;

import edu.ucdavis.gc.bm.clustering.Metricable;
import edu.ucdavis.gc.bm.descriptorGroup.Descriptor;
import edu.ucdavis.gc.pdb.PdbSubstitutionMatrix;

public class DescriptorSeqMetric extends Descriptor implements Metricable{
	
	public DescriptorSeqMetric(){
		
	}
	
	public DescriptorSeqMetric(Descriptor descriptor){
		super();
		super.setName(descriptor.getName());
		super.setSegments(descriptor.getSegments());
		super.setFoldAstral(descriptor.getFoldAstral());
		
	}

	@Override
	public Double distanceTo(Metricable o) {
		// TODO Auto-generated method stub
		Descriptor desc = (Descriptor) o;
		int numSegs = this.getNumberSegments();
		if(numSegs != desc.getNumberSegments()){
			return null;
		}
		List<String> segs1 = this.getSeqs();
		List<String> segs2 = desc.getSeqs();
		String seq1 = "";
		String seq2 = "";
		for(int i = 0; i < numSegs; i++){
			seq1 += segs1.get(i);
			seq2 += segs2.get(i);
		}
		if(seq1.length() != seq2.length()){
			return null;
		}
		Double result = 0.0;
		for(int i = 0; i < seq1.length(); i++){
			char a1 = seq1.charAt(i);
			char a2 = seq2.charAt(i);
			if (a1 == '.' || a1 == '-') {
				a1 = '*';
			}
			if (a2 == '.' || a2 == '-') {
				a2 = '*';
			}
			try {
				result += PdbSubstitutionMatrix.getScore(a1, a2);		
			} catch (NullPointerException e) {
				System.out.println("ai " + a1 + "   aj " + a2);
				throw e;
			}
		}
		return (double) -result/seq1.length(); // '-' should be selected in order to make this measure similar to distance:
		// less - better
	}

}
