package ebi.greg.eslr;

import java.sql.ResultSet;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

public class PdbUtils
{

	public static AlignmentResult getBestPdbMatch(String ensAcc, String uniAcc, ArrayList<String> ensemblPdbs)
			throws Exception
	{
		HashSet<String> pdbAccs = new HashSet<String>();
		pdbAccs.addAll(ensemblPdbs);

		List<String> xrefs = UniProtUtils.getPDBCrossReferences(uniAcc);
		for (String ref : xrefs)
		{
			pdbAccs.add(ref.toLowerCase());
		}

		if (pdbAccs.size() == 0)
			throw new BreakException("No PDB entries found for UniProt entry " + uniAcc);

		AlignmentResult bestAln = null;
		String bestPdbSeq = null;
		double bestScore = -10000;
		for (String pdbAcc : pdbAccs)
		{
			System.err.println("  "+pdbAcc);
			try
			{
				String pdbSeq = getPdbFinderLine(pdbAcc, "seq");
				/*
				 * Align the pdb sequence with the uniprot sequence.
				 * If this PDB's score is the best seen so far, store the relevant info.
				 */
				String uniSeq = UniProtUtils.getSequence(uniAcc);
				final AlignmentResult aln = JAlignerUtils.alignProteins(pdbAcc, uniAcc, pdbSeq, uniSeq);
				if (aln.getScore() >= bestScore)
				{
					bestAln = aln;
					bestScore = aln.getScore();
					bestPdbSeq = pdbSeq;
				}
			} catch (Exception e)
			{
				e.printStackTrace();
				continue;
			}
		}
		return bestAln;
	}

	public static String getPdbFinderLine(String pdbAcc, String lineKey) throws Exception
	{
		String returnMe = "";
		try
		{
			Statement st = EnsemblUtils.comparaConnection().createStatement();
			st.execute("use "+Config.inputDb+";");
			ResultSet rs = st.executeQuery("select " + lineKey + " from pdbfinder WHERE id=\"" + pdbAcc + "\";");
			rs.next();
			returnMe = rs.getString(1);
		} catch (Exception e)
		{
			throw new RuntimeException("Couldn't find Pdbfinder line " + lineKey);
		}
		return returnMe;
	}
}
