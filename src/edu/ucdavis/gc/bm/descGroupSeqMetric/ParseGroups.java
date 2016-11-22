package edu.ucdavis.gc.bm.descGroupSeqMetric;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.ucdavis.gc.bm.descriptorGroup.Descriptor;
import edu.ucdavis.gc.bm.descriptorGroup.Group;
import edu.ucdavis.gc.bm.descriptorGroup.Segment;
import edu.ucdavis.gc.bm.starthere.Main;

public class ParseGroups {

	private File file;

	public ParseGroups(File file) {
		this.file = file;
	}

	public ParseGroups(String filePath) {
		this.file = new File(filePath);
	}

	public ParseGroups(String dirName, String fileName) {
		File dir = new File(dirName);
		this.file = new File(dir, fileName);
	}

	/**
	 * parse file with groups of descriptors
	 * 
	 * @return list of groups
	 * @throws IOException
	 */
	public List<Group> parse() throws IOException {
		List<Group> result = new ArrayList<Group>();
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String line;
		List<String> groupLines = new ArrayList<String>(); // contains all lines
															// corresponding to
															// one group
		// boolean fStartGroup = false;
		boolean firstGroup = true;
		while ((line = reader.readLine()) != null) {
			if (line.trim().equals("")) { // if empty line - continue
				continue;
			}
			if (line.startsWith("GROUP")) { // line starts with "GROUP" -
											// indicator to start groupLines
				// fStartGroup = true;
				if (firstGroup == false) { // if it is not first group
					result.add(this.parseGroup(groupLines)); // add group to
																// result
				}
				firstGroup = false;
				groupLines.clear(); // clear groupLines
				groupLines.add(line); // add "GROUP" line
				continue;
			}
			groupLines.add(line);
		}
		result.add(this.parseGroup(groupLines)); // add last group to result
		reader.close();
		return result;
	}

	private Group parseGroup(List<String> groupLines) {
		Group result = new Group();
		List<Descriptor> descriptors = new ArrayList<Descriptor>();
		for (String line : groupLines) {
			if (line.startsWith("GROUP")) {
				result.setName(this.parseGroupLine(line));
			} else {
				DescriptorSeqMetric curDesc = this.parseDescriptorLine(line);
				if (curDesc != null) {
					descriptors.add(curDesc);
				}
			}
		}
		result.setDescriptors(descriptors);
		return result;
	}

	private String parseGroupLine(String line) { // "GROUP: 1ub3a_#203: 5"
		String[] tokens = line.trim().split("\\s+");
		return tokens[1].substring(0, tokens[1].indexOf(":"));
	}

	private DescriptorSeqMetric parseDescriptorLine(String line) { // "1c1yb_#109  97-103   AVFRLLH  107-111  GKKAR  125-131  ELQVDFL...  d.15.1.5 Human"
		String domain = line.substring(0, line.indexOf("#"));
		String chain = domain.substring(4, 5).toUpperCase();
		if (Main.problemDomains.contains(domain)) {
			return null;
		}
		DescriptorSeqMetric result = new DescriptorSeqMetric();
		String[] tokens = line.split("\\s+");
		result.setName(tokens[0]);
		List<Segment> segments = new ArrayList<Segment>();
		int i = 1;
		int previousEndFastaNo = -1000;
		while (!tokens[i].matches("[a-g]\\..*")) { // loop until meet "d.15.1.5"
													// - indicator of ASTRAL
													// fold
			/*
			 * Segment segment = new Segment();
			 * 
			 * if (!tokens[i + 1].matches("\\.*[A-Za-z]+\\.*")) { i++; }
			 * segment.setSeq(tokens[i + 1]); int indexSeparator =
			 * tokens[i].indexOf('-', 1); String start = tokens[i].substring(0,
			 * indexSeparator); String end = tokens[i].substring(indexSeparator
			 * + 1); // set fastaNoEnd and fastaNoStart - No's of start and end
			 * of // the // segment according to the enumeration of fasta
			 * sequence try { String startAddr = start; String endAddr = end; if
			 * (!chain.equals(".")) { startAddr = chain + startAddr; endAddr =
			 * chain + endAddr; } if (!startAddr.matches(".*[a-zA-Z]$")) {
			 * startAddr = startAddr + "_"; } if
			 * (!endAddr.matches(".*[a-zA-Z]$")) { endAddr = endAddr + "_"; }
			 * 
			 * if (Main.checkBounds) { int startFastaNo =
			 * Main.hashResMapStr2Int.get(domain).get( startAddr); int
			 * endFastaNo = Main.hashResMapStr2Int.get(domain).get( endAddr); //
			 * check for equality of lengths if segment without dots if
			 * ((endFastaNo - startFastaNo + 1) != tokens[i + 1] .length() &&
			 * tokens[i + 1].matches("[A-Z]+")) { System.err .println(line +
			 * "\n\t difference of length of segment sequence and calculated length based on segment bounds"
			 * ); return null; } if (endFastaNo < startFastaNo) {
			 * System.err.println(line +
			 * "\n\t right bound is less then left bound"); return null; } if
			 * (previousEndFastaNo >= startFastaNo) { System.err.println(line +
			 * "\n\t overlapping segments"); return null; } previousEndFastaNo =
			 * endFastaNo; } } catch (NullPointerException e) {
			 * System.err.println(line + "\n\t caught NullPointerException");
			 * return null; } segment.setStart(start); segment.setEnd(end);
			 */
			if (tokens[i].matches("[A-Z]*[0-9]+-[0-9]+[A-Z]*")
					&& tokens[i + 1].matches("\\.*[A-Za-z]+\\.*")) {
				// System.out.println(tokens[i] + "\t" + tokens[i+1]);
				Segment segment = parseSegment(domain, tokens[i], tokens[i + 1]);
				segments.add(segment);
				i += 2;
			} else if (tokens[i].matches("[A-Z]*[0-9]+-[0-9]+[A-Z]*")
					&& tokens[i + 1].matches("[A-Z]*[0-9]+-[0-9]+[A-Z]*")
					&& tokens[i + 2].matches("\\.*[A-Za-z]+\\.*")) {
				i += 3;
				return null;
			}

			// i += 2;
		}
		result.setSegments(segments);
		result.setFoldAstral(tokens[i]);
		return result;
	}

	private Segment parseSegment(String domain, String bounds, String sequence) {
		String chain = domain.substring(4, 5).toUpperCase();
		Segment segment = new Segment();
		int indexSeparator = bounds.indexOf('-', 1);
		String start = bounds.substring(0, indexSeparator);
		String end = bounds.substring(indexSeparator + 1);
		String startAddr = start;
		String endAddr = end;
		if (!chain.equals(".")) {
			startAddr = chain + startAddr;
			endAddr = chain + endAddr;
		}
		if (!startAddr.matches(".*[a-zA-Z]$")) {
			startAddr = startAddr + "_";
		}
		if (!endAddr.matches(".*[a-zA-Z]$")) {
			endAddr = endAddr + "_";
		}

		// adjust the bounds of segments
		Pattern patternBefore = Pattern.compile("[a-z]+[A-Z]");
		Matcher matcherBefore = patternBefore.matcher(sequence);
		Pattern patternAfter = Pattern.compile("[A-Z][a-z]+");
		Matcher matcherAfter = patternAfter.matcher(sequence);
		// find the first occurrence of the patternBefore
		int beforeLength = 0;
		if (matcherBefore.find()) {
			beforeLength = matcherBefore.end() - matcherBefore.start() - 1;
		}
		int afterLength = 0;
		while (matcherAfter.find()) {
			afterLength = matcherAfter.end() - matcherAfter.start() - 1;
		}

		int startFastaNo = Main.hashResMapStr2Int.get(domain).get(startAddr);
		int endFastaNo = Main.hashResMapStr2Int.get(domain).get(endAddr);

		// adjust indexes
		String newStart = Main.hashResMapInt2Str.get(domain).get(
				startFastaNo - beforeLength);
		String newEnd = Main.hashResMapInt2Str.get(domain).get(
				endFastaNo + afterLength);
		// if one of the residue addresses is actually corresponds to the gap
		// residue address of a gap is noted as _0_
		if (newStart.equalsIgnoreCase("_0_") || newEnd.equalsIgnoreCase("_0_")) {
			System.err.println(sequence);
			sequence = sequence.replaceAll("[a-z]", ".");
			newStart = Main.hashResMapInt2Str.get(domain).get(startFastaNo);
			newEnd = Main.hashResMapInt2Str.get(domain).get(endFastaNo);
		} else {
			startFastaNo -= beforeLength;
			endFastaNo += afterLength;
		}

		if (!chain.equals(".")) {
			newStart = newStart.substring(1);
			newEnd = newEnd.substring(1);
		}
		newStart = newStart.replaceAll("_", "");
		newEnd = newEnd.replaceAll("_", "");
		segment.setStart(newStart);
		segment.setEnd(newEnd);
		segment.setSeq(sequence);
		return segment;
	}

	/*
	 * TODO
	 */
	private Segment parseSegment(String domain, String bounds1, String bounds2,
			String sequence) {
		String chain = domain.substring(4, 5).toUpperCase();
		Segment segment = new Segment();
		int indexSeparator = bounds1.indexOf('-', 1);
		String start = bounds1.substring(0, indexSeparator);
		String end = bounds1.substring(indexSeparator + 1);
		String startAddr = start;
		String endAddr = end;
		if (!chain.equals(".")) {
			startAddr = chain + startAddr;
			endAddr = chain + endAddr;
		}
		if (!startAddr.matches(".*[a-zA-Z]$")) {
			startAddr = startAddr + "_";
		}
		if (!endAddr.matches(".*[a-zA-Z]$")) {
			endAddr = endAddr + "_";
		}

		int startFastaNo = Main.hashResMapStr2Int.get(domain).get(startAddr);
		int endFastaNo = Main.hashResMapStr2Int.get(domain).get(endAddr);
		segment.setStartFastaNo(startFastaNo);
		segment.setEndFastaNo(endFastaNo);
		return segment;
	}

	public static void main(String[] args) {

		String text = "...sdghjASDFADTDTSDFSGSCshdfjsk..";
		Pattern patternBefore = Pattern.compile("[a-z]+[A-Z]");
		Matcher matcherBefore = patternBefore.matcher(text);
		Pattern patternAfter = Pattern.compile("[A-Z][a-z]+");
		Matcher matcherAfter = patternAfter.matcher(text);
		// find the first occurrence of the patternBefore
		if (matcherBefore.find()) {
			int beforeLength = matcherBefore.end() - matcherBefore.start() - 1;
			System.out.println("StartBefore: " + matcherBefore.start());
			System.out.println("EndBefore: " + matcherBefore.end());
			System.out.println(beforeLength);

		}
		while (matcherAfter.find()) {
			System.out.println("StartAfter: " + matcherAfter.start());
			System.out.println("EndAfter: " + matcherAfter.end());
			int afterLength = matcherAfter.end() - matcherAfter.start() - 1;
			System.out.println(afterLength);
		}

	}
}
