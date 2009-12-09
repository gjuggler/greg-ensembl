package ebi.greg.eslr;

import java.sql.ResultSet;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import uk.ac.ebi.kraken.interfaces.uniprot.DatabaseCrossReference;
import uk.ac.ebi.kraken.interfaces.uniprot.DatabaseType;
import uk.ac.ebi.kraken.interfaces.uniprot.UniProtEntry;
import uk.ac.ebi.kraken.interfaces.uniprot.dbx.hssp.Hssp;
import uk.ac.ebi.kraken.interfaces.uniprot.dbx.pdb.Pdb;

public class PdbUtils
{

	public static AlignmentResult getBestPdbMatch(String ensAcc, UniProtEntry entry, ArrayList<String> ensemblPdbs)
			throws Exception
	{
		HashSet<String> pdbAccs = new HashSet<String>();
		pdbAccs.addAll(ensemblPdbs);

		List<DatabaseCrossReference> xrefs = entry.getDatabaseCrossReferences(DatabaseType.HSSP);
		xrefs.addAll(entry.getDatabaseCrossReferences(DatabaseType.PDB));
		for (DatabaseCrossReference ref : xrefs)
		{
			if (ref instanceof Hssp)
			{
				if (Main.DEBUG)
					System.err.println("  hssp");
				Hssp hssp = (Hssp) ref;
				String pdbAcc = hssp.getHsspDescription().getValue();
				pdbAccs.add(pdbAcc.toLowerCase());
			} else if (ref instanceof Pdb)
			{
				if (Main.DEBUG)
					System.err.println("  pdb");
				Pdb pdb = (Pdb) ref;
				String pdbAcc = pdb.getPdbAccessionNumber().getValue().toLowerCase();
				pdbAccs.add(pdbAcc.toLowerCase());
			}
		}

		if (pdbAccs.size() == 0)
			throw new BreakException("No PDB entries found for UniProt entry " + entry);

		AlignmentResult bestAln = null;
		String bestPdbSeq = null;
		double bestScore = -10000;
		for (String pdbAcc : pdbAccs)
		{
			System.err.println("  "+pdbAcc);
			try
			{
				String uniAcc = UniProtUtils.getAccession(entry);
				String pdbSeq = getPdbFinderLine(pdbAcc, "seq");
				/*
				 * Align the pdb sequence with the uniprot sequence.
				 * If this PDB's score is the best seen so far, store the relevant info.
				 */
				String uniSeq = entry.getSequence().getValue();
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
