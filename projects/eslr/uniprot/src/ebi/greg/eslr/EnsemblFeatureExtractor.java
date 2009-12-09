package ebi.greg.eslr;

import java.util.ArrayList;
import java.util.List;

import ebi.greg.eslr.FeatureExtractionUtils.RegionTag;
import ebi.greg.eslr.JAlignerUtils.AlignmentWithMaps;

public class EnsemblFeatureExtractor {
	String species;
	String genomeDb;

	int memberId;
	int nodeId;

	AlignmentResult uniEns = null;
	AlignmentWithMaps uniEnsMap = null;

	String ensAcc;
	String uniAcc;
	String pdbAcc;

	String ensSeq;
	String uniSeq;
	String pdbSeq;

	boolean hasEns;
	boolean hasUni;
	boolean hasPdb;

//	UniProtEntry uniEntry;

	public EnsemblFeatureExtractor(String ensemblAcc) {
		this.ensAcc = ensemblAcc;
		try {
			getFeatures();
		} catch (BreakException e) {
			e.printStackTrace();
			return;
		} catch (Exception ex) {
			ex.printStackTrace();
			return;
		}
	}

	public void getFeatures() throws Exception {
		String ensSeq = null;
		try {
			ensSeq = EnsemblUtils.getProteinSequence(ensAcc);
			hasEns = true;
		} catch (Exception e) {
			throw e;
		}

		species = EnsemblUtils.getSpecies(ensAcc);
		genomeDb = EnsemblUtils.getGenomeDb(species);
		if (genomeDb.length() < 1) {
			throw new RuntimeException(
					"No Ensembl genome DB found for species: " + species + "!");
		}
		System.err.println(ensAcc);
		memberId = EnsemblUtils.getMemberId(ensAcc);

		/*
		 * Find the best corresponding UniProt entry for this Ensembl peptide.
		 */
		try {
			uniEns = EnsemblUtils.getBestUniProtMatch(ensAcc, ensSeq, genomeDb);
			uniAcc = uniEns.getName1();
			uniEnsMap = new AlignmentWithMaps(uniEns);
//			uniEntry = UniProtUtils.getEntry(uniAcc);
			uniSeq = UniProtUtils.getSequence(uniAcc);
//			uniSeq = uniEntry.getSequence().getValue();
			if (uniSeq != null) {
				hasUni = true;
			}
		} catch (BreakException e) {
			// Do nothing special.
		}
		
		/*
		 * Find and align the best PDB match.
		 */
		ArrayList<String> ensemblPdbXrefs = EnsemblUtils.getPdbMatches(ensAcc,
				genomeDb);

		String pdbAcc = null;
		AlignmentResult pdbUni = null;
		AlignmentWithMaps pdbUniMap = null;
		String pdbSeq = null;
		String pdbAccess = null;
		String pdbDssp = null;
		if (hasUni) {
			try {
				pdbUni = PdbUtils.getBestPdbMatch(ensAcc, uniAcc,
						ensemblPdbXrefs);
				pdbAcc = pdbUni.getName1();
				pdbUniMap = new AlignmentWithMaps(pdbUni);
				pdbSeq = pdbUni.getOriginalSequence1();
				pdbAccess = PdbUtils.getPdbFinderLine(pdbAcc, "access");
				pdbDssp = PdbUtils.getPdbFinderLine(pdbAcc, "dssp");

				if (pdbSeq != null) {
					hasPdb = true;
				}
			} catch (BreakException be) {
				// Don't do anything special.
			} catch (Exception e) {
				throw e;
			}
		}

		StringBuffer ensBuff = new StringBuffer();
		StringBuffer uniBuff = new StringBuffer();
		StringBuffer pdbBuff = new StringBuffer();

		/*
		 * Grab the features for the given UniProt entry, if one exists.
		 */
		if (hasUni) {
			List<RegionTag> regions = FeatureExtractionUtils
					.getFeaturesFromUniProt(uniAcc);

			ArrayList<PosTagValue> tagvals = new ArrayList<PosTagValue>();
			for (int i = 0; i < ensSeq.length(); i++) {
				int ensPos = i;
				int uniPos = uniEnsMap.bToA[ensPos];
				if (uniPos >= 0) {
					uniBuff.append(uniSeq.charAt(uniPos));
					for (RegionTag region : regions) {
						if (region.seq_end >= uniPos
								&& region.seq_start <= uniPos) {
							int featurePos = uniPos - region.seq_start;
							PosTagValue ptv = new PosTagValue(ensPos, uniAcc,
									"UNIPROT", region.tag, region.value, uniSeq
											.charAt(uniPos - 1), featurePos);
							tagvals.add(ptv);
						}
					}
				} else {
					uniBuff.append("-");
				}
			}
			if (tagvals.size() > 0) {
				if (Main.DEBUG) {
					System.err.println("  --> Found " + tagvals.size()
							+ " UniProt tags.");
					System.out.println(PosTagValue.toTSV(tagvals));
				}
			}
		}
		
		if (hasPdb) {
			// This call does most of the work in terms of getting the domain
			// annotations.
			List<RegionTag> regions = FeatureExtractionUtils
					.getFeaturesFromPdb(pdbAcc, pdbDssp, pdbAccess);

			ArrayList<PosTagValue> tagvals = new ArrayList<PosTagValue>();
			for (int i = 0; i < ensSeq.length(); i++) {
				int ensPos = i;
				int uniPos = uniEnsMap.bToA[ensPos];
				if (uniPos >= 0) {
					int pdbPos = pdbUniMap.bToA[uniPos];
					if (pdbPos >= 0) {
						char pdbChar = pdbSeq.charAt(pdbPos);
						pdbBuff.append(pdbChar);
						for (RegionTag region : regions) {
							if (region.seq_end >= pdbPos
									&& region.seq_start <= pdbPos) // If the
																	// regio
																	// noverlaps
							{
								int featurePos = pdbPos - region.seq_start;
								PosTagValue ptv = new PosTagValue(ensPos + 1,
										pdbAcc, "PDB", region.tag,
										region.value, pdbSeq.charAt(pdbPos),
										featurePos);
								tagvals.add(ptv);
							}
						}
					} else {
						pdbBuff.append("-");
					}

				} else {
					pdbBuff.append("-");
				}
			}
			if (tagvals.size() > 0) {

				if (Main.DEBUG) {
					System.err.println("  --> Found " + tagvals.size()
							+ " PDB tags.");
					System.out.println(PosTagValue.toTSV(tagvals));
				}

			}
		}

		// Print out a portion of each sequence, just for visual effect.
		int len = ensSeq.length();
		int lo = 12;
		int hi = ensSeq.length() - 12;
		System.err.println("  " + ensSeq.substring(0, lo) + "..."
				+ ensSeq.substring(hi, len) + "  " + ensAcc);
		if (hasUni)
			System.err.println("  " + uniBuff.substring(0, lo) + "..."
					+ uniBuff.substring(hi, len) + "  " + uniAcc);
		if (hasPdb)
			System.err.println("  " + pdbBuff.substring(0, lo) + "..."
					+ pdbBuff.substring(hi, len) + "  " + pdbAcc);
	}

}
