package edu.ucdavis.gc.bm.starthere;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import edu.ucdavis.gc.bm.clustering.AddRedundancyClustering;
import edu.ucdavis.gc.bm.clustering.HierarchicalClustering;
import edu.ucdavis.gc.bm.clustering.Metricable;
import edu.ucdavis.gc.bm.descGroupSeqMetric.DescriptorSeqMetric;
import edu.ucdavis.gc.bm.descGroupSeqMetric.GroupSeqSignalNaiv;
import edu.ucdavis.gc.bm.descriptorGroup.Descriptor;
import edu.ucdavis.gc.bm.descriptorGroup.Group;
import edu.ucdavis.gc.bm.descriptorGroup.ParseD3KGroups;
import edu.ucdavis.gc.bm.properties.Props;

public class Main {

	//private static String descGrDir = "/home/bohdan/DescriptorsGroups/seqClustering/";
	private static String descGrDir = null; // "/home/bohdan/DescriptorsGroups/DIIIK/";

	//private static String descGrFile = "tmp_groups_ASTRAL175v3_827astral.txt";
	private static String descGrFile = null; // "DescriptorGroups.827astral.2.5.fUp.f7.s7.d7.txt.tmp1";
	
	private static String outDir = "";
	
	private static String outFile = null;

	public static TreeSet<String> problemDomains;
	/**
	 * HashMap of mapping residues addresses to fastaNo's for all astral domains 
	 */
	public static HashMap<String, HashMap<String, Integer>> hashResMapStr2Int;
	
	/**
	 * HashMap of mapping fastaNo's to residues' addresses for all astral domains
	 */
	public static HashMap<String,HashMap<Integer,String>> hashResMapInt2Str;
	
	public static HashMap<String, String> fastaSeqs;
	public static HashMap<String, String> SS_Seqs;
	
	private static Double cutOff_1 = 0.2;
	
	private  static Double cutOff_2 = 0.4;
	
	private static int minNumDescs = 1;
	
	public static boolean checkBounds = false;
	/**
	 * fOutput - parameter to handle the format of output
	 * 1 - use "Clusters" and indexes of in descriptors (should be used for debugging )
	 * 0 - format of Groups txt file (should be used for generating groups )
	 */
	private static int fOutput = 1;

	/**
	 * @param args
	 * @throws FileNotFoundException 
	 */
	public static void main(String[] args) throws FileNotFoundException {

		// read args
		for (int i = 0; i < args.length; i++) {
			if (args[i].equalsIgnoreCase("-o")
					|| args[i].equalsIgnoreCase("-outDir")
					|| args[i].equalsIgnoreCase("--outDir")) {
				outDir = args[i + 1] + "/";
			}
			if (args[i].equalsIgnoreCase("-d")
					|| args[i].equalsIgnoreCase("-descGrDir")
					|| args[i].equalsIgnoreCase("--descGrDir")) {
				descGrDir = args[i + 1] + "/";
			}
			if (args[i].equalsIgnoreCase("-i")
					|| args[i].equalsIgnoreCase("-descGrFile")
					|| args[i].equalsIgnoreCase("--descGrFile")) {
				descGrFile = args[i + 1] ;
			}
			if (args[i].equalsIgnoreCase("-c1")
					|| args[i].equalsIgnoreCase("-cutOff1")
					|| args[i].equalsIgnoreCase("--cutOff1")) {
				cutOff_1 = Double.valueOf(args[i + 1]);
			}			
			if (args[i].equalsIgnoreCase("-c2")
					|| args[i].equalsIgnoreCase("-cutOff2")
					|| args[i].equalsIgnoreCase("--cutOff2")) {
				cutOff_2 = Double.valueOf(args[i + 1]);
			}		
			if (args[i].equalsIgnoreCase("-m")
					|| args[i].equalsIgnoreCase("-minNumDescs")
					|| args[i].equalsIgnoreCase("--minNumDescs")) {
				minNumDescs = Integer.valueOf(args[i + 1]);
			}		
			if (args[i].equalsIgnoreCase("-fOutput")
					|| args[i].equalsIgnoreCase("-fOutput")
					|| args[i].equalsIgnoreCase("--fOutput")) {
				fOutput = Integer.valueOf(args[i + 1]);
			}			
			
		}
			
		// problem domains has "XXXXXXXXXX" in fasta sequence the program can
		// not handle with them
		problemDomains = new TreeSet<String>();
		problemDomains.add("1gw5a_");
		problemDomains.add("1gw5b_");
		problemDomains.add("1uc8a2");
		problemDomains.add("1pya.1");
		problemDomains.add("2fpwa1"); // this protein sequence has just one "X"
										// but the program couldn't handle the
										// problem

		hashResMapStr2Int = null;
		FileInputStream fis = null;
		ObjectInputStream in = null;
		String astralResMap_ser_filepath = Props.get("hashResMapStr2Int"); 
		//"/home/bohdan/workspace/DescGroupUtils/astral1.75.ResMap.Str2Int.java.ser"; // props.get("resMap_ser_path");
		try {
			fis = new FileInputStream(astralResMap_ser_filepath);
			in = new ObjectInputStream(fis);
			hashResMapStr2Int = (HashMap<String, HashMap<String, Integer>>) in
					.readObject();
			in.close();
		} catch (IOException ex) {
			ex.printStackTrace();
		} catch (ClassNotFoundException ex) {
			ex.printStackTrace();
		}

		hashResMapInt2Str = null;
		fis = null;
		in = null;
		astralResMap_ser_filepath = Props.get("hashResMapInt2Str"); 
		//"/home/bohdan/workspace/DescGroupUtils/astral1.75.ResMap.Int2Str.java.ser"; // props.get("resMap_ser_path");
		try {
			fis = new FileInputStream(astralResMap_ser_filepath);
			in = new ObjectInputStream(fis);
			hashResMapInt2Str = (HashMap<String, HashMap<Integer, String>>) in
					.readObject();
			in.close();
		} catch (IOException ex) {
			ex.printStackTrace();
		} catch (ClassNotFoundException ex) {
			ex.printStackTrace();
		}
		
		
		fis = null;
		in = null;
		String astralFastaSeq_ser_filepath = Props.get("fastaSeqs");
		//"/home/bohdan/workspace/DescGroupUtils/astral1.75.fastaSeqs.java.ser";// props.get("fastaSeq_ser_path");
		try {
			fis = new FileInputStream(astralFastaSeq_ser_filepath);
			in = new ObjectInputStream(fis);
			fastaSeqs = (HashMap<String, String>) in.readObject();
			in.close();
		} catch (IOException ex) {
			ex.printStackTrace();
		} catch (ClassNotFoundException ex) {
			ex.printStackTrace();
		}

		fis = null;
		in = null;
		String astralSsSeq_ser_filepath = Props.get("SsSeqs");
		//"/home/bohdan/workspace/DescGroupUtils/astral1.75.fastaSeqs.java.ser";// props.get("fastaSeq_ser_path");
		try {
			fis = new FileInputStream(astralSsSeq_ser_filepath);
			in = new ObjectInputStream(fis);
			SS_Seqs = (HashMap<String, String>) in.readObject();
			in.close();
		} catch (IOException ex) {
			ex.printStackTrace();
		} catch (ClassNotFoundException ex) {
			ex.printStackTrace();
		}
		cutOff_1 = Props.getDouble("cutOff_1");
		
		cutOff_2 = Props.getDouble("cutOff_2");
		
		minNumDescs = Props.getInt("minNumDescs");
		
		descGrDir = Props.get("descGrDir");
		descGrFile = Props.get("descGrFile");
		outFile = Props.get("outFile");
		fOutput = Props.getInt("fOutput");
		
		ParseD3KGroups parserGroups = new ParseD3KGroups( descGrDir, descGrFile, hashResMapStr2Int, hashResMapInt2Str, fastaSeqs, SS_Seqs, problemDomains);
		
		List<Group> groups = null;
		try {
			groups = parserGroups.parse();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		int all = groups.size();


		
/*		for (Group group : groups) {
			if(true)break;

			System.out.println("\n"
					+ "################################################\n\n"
					+ "The clustering for group " + group.getName() + "("
					+ group.getNumberMembers() + ")" + " is started.\n");

			GroupSeqSignalNaiv grSeqSign = new GroupSeqSignalNaiv(group);

			List<Metricable> descs = new ArrayList<Metricable>();
			for (Descriptor desc : group.getDescriptors()) {
				descs.add((DescriptorSeqMetric) desc);
			}

			for (int nC = 5; nC <= 20; nC += 5) {
				HierarchicalClustering clustering = new HierarchicalClustering(
						descs);
				
				clustering.process(nC);

				List<Set<Integer>> clusters = clustering.getClusters();

				for (int i = 0; i < clusters.size(); i++) { // Set<Integer>
															// cluster
															// : clusters){
					System.out.println("Cluster # " + (i + 1) + "(" + nC + ")"
							+ ":\t" + grSeqSign.getSeqSignal(clusters.get(i))
							+ " of " + grSeqSign.getSeqSignal());
					for (int index : clusters.get(i)) {
						System.out.println(descs.get(index).toString());
					}
					System.out.println();
				}
				System.out
						.println("===================================================");
			}

			try {
				InputStreamReader isr = new InputStreamReader(System.in);
				BufferedReader br = new BufferedReader(isr);
				System.out.println("The clusterring for group "
						+ group.getName() + " is terminated.\n");
				String s = br.readLine();

			} catch (IOException e) {
				e.printStackTrace();
			}
		}
*/		
		PrintWriter pw = new PrintWriter(new File(outFile));
		
		for (Group group : groups) {
			 
			List<Descriptor> metrDescs = new ArrayList<Descriptor>();
			for(Descriptor desc : group.getDescriptors()){
				metrDescs.add(new DescriptorSeqMetric(desc));
			}
			group.setDescriptors(metrDescs);
			if(fOutput == 1){
			System.out.println("\n"
					+ "################################################\n\n"
					+ "The clustering for group " + group.getName() + "("
					+ group.getNumberMembers() + ")" + " is started.\n");
			}

			GroupSeqSignalNaiv grSeqSign = new GroupSeqSignalNaiv(group);

			List<Metricable> descs = new ArrayList<Metricable>();
			for (Descriptor desc : group.getDescriptors()) {
				descs.add((DescriptorSeqMetric) desc);
			}

			for (double cutOff1 = cutOff_1; cutOff1 <= cutOff_1; cutOff1 += 0.2) {
				HierarchicalClustering clustering = new HierarchicalClustering(
						descs);
				
				clustering.process(-cutOff1); // in clustering process it's necessary to use "-"
				List<Set<Integer>> clusters = clustering.getClusters();
/*
				System.out.println("Phase 1");
				for (int i = 0; i < clusters.size(); i++) { // Set<Integer>
															// cluster
															// : clusters){
					System.out.printf("Cluster # %d (%2.1f):\t%7.5f of %7.5f\n", (i + 1),  cutOff1 ,
							 grSeqSign.getSeqSignal(clusters.get(i)), grSeqSign.getSeqSignal());
					for (int index : clusters.get(i)) {
						System.out.println(index + "\t" + descs.get(index).toString() );
					}
					System.out.println();
				}
	*/			
				AddRedundancyClustering adred = new AddRedundancyClustering(clusters,descs);
				//for(double cutOff2 = 0.2; cutOff2 <= cutOff1; cutOff2 += 0.2){
				double cutOff2 = cutOff_2;
					//System.out.println("Phase2");
					adred.process(-cutOff2);
					List<Set<Integer>> clustersRed = adred.getClusters();

					for (int i = 0; i < clustersRed.size(); i++) { // Set<Integer>
																// cluster
																// : clusters){
						if(clustersRed.get(i).size() < minNumDescs){
							continue;
						}
						if(fOutput == 1){
						System.out.printf("Cluster # %d (%2.1f : %2.1f):\t%7.5f of %7.5f\n", (i + 1),  cutOff1 , cutOff2,
								 grSeqSign.getSeqSignal(clustersRed.get(i)), grSeqSign.getSeqSignal());
						}
						if(fOutput == 0){
							pw.printf("GROUP: %s_%d: %d\n", group.getName(), (i+1), clustersRed.get(i).size());
						}
						
						for (int index : clustersRed.get(i)) {
							if(fOutput == 1){
								if(clusters.get(i).contains(index)){
									System.out.println(index + "*\t" + descs.get(index).toString() );
								}else{
									System.out.println(index + "\t" + descs.get(index).toString() );
								}							
							}
							if(fOutput == 0){
								pw.println(descs.get(index).toString() );
							}
						}
						System.out.println();
						pw.println();
					}
				//}
				if(fOutput == 1){
					System.out
					.println("===================================================");					
				}
			}
			//if(true)break;
			if(fOutput == 1){
				System.err.println("The clustering for group "
						+ group.getName() + " is terminated.\n\n");
				try {
					InputStreamReader isr = new InputStreamReader(System.in);
					BufferedReader br = new BufferedReader(isr);
					String s = br.readLine();
					
				} catch (IOException e) {
					e.printStackTrace();
				}			
			}
			
		}
		pw.close();
	}

}
