package ebi.greg.eslr;

import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import uk.ac.ebi.kraken.interfaces.uniprot.UniProtEntry;

/**
 * A collection of utility methods for Ensembl stuff.
 * 
 * @author Greg
 * 
 */
public class EnsemblUtils {
	/**
	 * Connection properties.
	 */
	static String anonDb = "ensembl_compara_50";

	/**
	 * Connection objects.
	 */
	static Connection anonConnection;
	static Connection comparaConnection;
	static ArrayList<String> databases = new ArrayList<String>();

	static PreparedStatement preparedStatement;

	/**
	 * Connect to the database.
	 */
	static {
		try {
			if (Main.DEBUG) {
				System.err.println("Starting Ensembl...");
			}

			// comparaConnection = SQLUtils
			// .getConnectionFromConfigFile(Config.comparaConfig);
			// Cache a list of all the Ensembl database names FROM THE ANONYMOUS
			// DATABASE

		} catch (Exception e) {
			e.printStackTrace();
			System.err
					.println("Could not connect to the Ensembl database(s). Quitting...");
			System.exit(0);
		}
		if (Main.DEBUG)
			System.err.println("Done!");
	}

	public static Connection anonConnection() {
		if (anonConnection == null) {
			try {
				anonConnection = SQLUtils
						.getConnectionFromString("mysql://anonymous:@ensembldb.ensembl.org:5306/");
				Statement s = anonConnection.createStatement();
				ResultSet rs = s.executeQuery("show databases;");
				while (rs.next()) {
					String db = rs.getString(1);
					// System.err.println(db);
					databases.add(db);
				}
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		return anonConnection;
	}

	public static Connection comparaConnection() {
		if (comparaConnection == null) {
			try {
				System.err.println("Compara connection...");
				comparaConnection = SQLUtils
						.getConnectionFromString(Config.comparaUrl);
			} catch (Exception e) {

			}
		}
		return comparaConnection;

	}

	public static String comparaDb() {
		return Config.inputDb;
	}

	/**
	 * Finds the best UniProt match, given an ENSP accession and a list of
	 * potential UniProt accessions.
	 * 
	 * Returns an array of strings: array[0] is the best UniProt accession
	 * array[1] is the UniProt alignment line array[2] is the Ensembl alignment
	 * line
	 * 
	 * @param ensemblPep
	 * @param uAccs
	 */
	static synchronized AlignmentResult getBestMatchImpl(String eAcc,
			String ensemblPep, String[] uAccs) throws Exception {
		AlignmentResult bestAln = null;
		double bestScore = -10000;
		for (String uniProt : uAccs) {
			// System.out.println("  "+uniProt);
			UniProtEntry entry = UniProtUtils.getEntry(uniProt);
			if (entry == null) {
				System.err.println("Null entry:" + uniProt);
				continue;
			} else {
				String uniPep = entry.getSequence().getValue();
				String thisAcc = entry.getPrimaryUniProtAccession().getValue();
				AlignmentResult aln = JAlignerUtils.alignProteins(uniProt,
						eAcc, uniPep, ensemblPep);
				if (aln.score >= bestScore) {
					bestAln = aln;
					bestScore = aln.score;
				}
				aln = null;
				uniPep = null;
			}
		}
		if (bestAln == null)
			throw new BreakException("No good UniProt alignment found!");
		return bestAln;
	}

	public static String getGenomeDb(String orgName) {
		String dbName = orgName.toLowerCase().replaceAll(" ", "_");
		/*
		 * Find the latest database in Ensembl for the given taxon.
		 */
		String latestDb = "";
		for (String db : databases) {
		    //System.out.println(db);
			db = db.toLowerCase();
			if (db.contains(dbName) && db.contains("_core_")
					&& db.matches(dbName + "_core_\\d+.*")) {
				if (db.compareTo(latestDb) > 0 || latestDb.length() == 0) {
					latestDb = db;
				}
			}
		}
		return latestDb;
	}

	static Pattern p = Pattern.compile("(\\d*?)([MD])",
			Pattern.CASE_INSENSITIVE);

	public static ArrayList<Integer> memberToAlignmentMap(int memberId,
			String ensSeq) throws Exception {
		String cigarLine = getCigarLine(memberId);
		if (cigarLine == null) {
			System.err.println("Null cigar!");
		}
		// Parse the cigar line to be able to map Ensembl coords to sitewise_aln
		// coords.
		Matcher m = p.matcher(cigarLine);
		// System.out.println(cigarLine);
		// StringBuffer sb = new StringBuffer();
		ArrayList<Integer> memberToAlignment = new ArrayList<Integer>();
		int count = 0;
		while (m.find()) {
			// System.out.println(m.group());
			String numS = m.group(1);
			int num = 1;
			if (numS.length() > 0)
				num = Integer.parseInt(numS);
			String type = m.group(2);
			for (int i = 0; i < num; i++) {
				if (type.equals("M")) {
					memberToAlignment.add(count);
					// sb.append('#');
				} else {
					// sb.append('-');
				}
				count++;
			}
		}
		return memberToAlignment;
	}

	public static AlignmentResult getBestUniProtMatch(String tx_stable_id,
			String tx_seq, String latestDb) throws Exception {
		/*
		 * Fetch the list of potential UniProt sequences.
		 */
		Statement s = anonConnection().createStatement();
		s.execute("use " + latestDb + ";");
		String query = "select t.stable_id,xdb.db_name,x.dbprimary_acc,x.display_label from translation_stable_id t, object_xref ox, xref x, external_db xdb"
				+ " where t.stable_id="
				+ qw(tx_stable_id)
				+ " and ox.ensembl_id=t.translation_id"
				+ " and ox.ensembl_object_type="
				+ qw("Translation")
				+ " and ox.xref_id=x.xref_id"
				+ " and xdb.external_db_id=x.external_db_id"
				+ " and xdb.db_name like " + qw("%Uniprot%") + " ;";
		// System.out.println("use "+latestDb+";"+query);
		ResultSet rs = s.executeQuery(query);
		int count = SQLUtils.rowCount(rs);
		String[] uAccs = new String[count];
		int i = 0;
		while (rs.next()) {
			uAccs[i] = rs.getString("x.dbprimary_acc");
			i++;
		}
		if (uAccs.length == 0)
			throw new BreakException("No UniProt matches found for "
					+ tx_stable_id);

		/*
		 * Get the best UniProt match for the given Ensembl peptide.
		 */
		if (tx_seq == null)
			tx_seq = EnsemblUtils.getProteinSequence(tx_stable_id);
		AlignmentResult result = getBestMatchImpl(tx_stable_id, tx_seq, uAccs);
		return result;
	}

	public static ArrayList<String> getPdbMatches(String stable_id,
			String latestDb) throws Exception {
		Statement s = anonConnection().createStatement();
		s.execute("use " + latestDb + ";");
		String query = "select t.stable_id,xdb.db_name,x.dbprimary_acc,x.display_label from translation_stable_id t, object_xref ox, xref x, external_db xdb"
				+ " where t.stable_id="
				+ qw(stable_id)
				+ " and ox.ensembl_id=t.translation_id"
				+ " and ox.ensembl_object_type="
				+ qw("Translation")
				+ " and ox.xref_id=x.xref_id"
				+ " and xdb.external_db_id=x.external_db_id"
				+ " and xdb.db_name like " + qw("%PDB%") + " ;";
		ResultSet rs = s.executeQuery(query);
		ArrayList<String> pdbAccs = new ArrayList<String>();
		int i = 0;
		while (rs.next()) {
			pdbAccs.add(rs.getString("x.dbprimary_acc").toLowerCase());
		}
		return pdbAccs;
	}

	/**
	 * Returns the peptide sequence for a given stable_id, using the compara
	 * tables.
	 * 
	 * @param stable_id
	 * @return
	 */
	public static String getSequence(String stable_id) throws Exception {
		Statement s = comparaConnection().createStatement();
		s.execute("use " + comparaDb());
		ResultSet rs = s
				.executeQuery("select s.sequence from member as m, sequence as s"
						+ " where m.sequence_id=s.sequence_id and m.stable_id="
						+ qw(stable_id) + ";");
		rs.next();
		return rs.getString(1);
	}

	/**
	 * Returns the protein sequence of a given ensembl peptide (i.e. ENSP...)
	 * 
	 * @param pepAcc
	 * @return
	 */
	static String getProteinSequence(String pepAcc) throws Exception {
		Statement s = comparaConnection().createStatement();
		s.execute("use " + comparaDb());
		ResultSet rs;
		if (pepAcc.contains("G00")) {
			rs = s
					.executeQuery("select s.sequence from member as gm, member as pm, sequence as s"
							+ " where pm.sequence_id=s.sequence_id and gm.stable_id="
							+ qw(pepAcc)
							+ " and gm.member_id=pm.gene_member_id");
		} else {
			rs = s
					.executeQuery("select s.sequence from member as m, sequence as s "
							+ "where m.sequence_id=s.sequence_id and m.stable_id="
							+ qw(pepAcc));
		}
		if (!rs.last())
			throw new BreakException("Cannot find peptide " + pepAcc);
		return rs.getString(1);
	}

	static String getProteinAcc(String geneAcc) throws Exception {
		Statement s = comparaConnection().createStatement();
		s.execute("use " + comparaDb());
		ResultSet rs;
		rs = s
				.executeQuery("select pm.member_id from member as gm, member as pm where gm.stable_id="
						+ qw(geneAcc) + " and gm.member_id=pm.gene_member_id");
		return rs.getString(1);
	}

	static String getSpecies(String pepAcc) throws Exception {
		Statement s = comparaConnection().createStatement();
		s.execute("use " + comparaDb());
		ResultSet rs = s
				.executeQuery("select db.name from member as m, genome_db as db where db.taxon_id=m.taxon_id and m.stable_id="
						+ qw(pepAcc));
		rs.next();
		try {
			return rs.getString(1);
		} catch (Exception e) {
			throw new BreakException("No species db for " + pepAcc);
		}
	}

	static String qw(String s) {
		return "\"" + s + "\"";
	}

	public static int getMemberId(String eAcc) throws Exception {
		Statement s = EnsemblUtils.comparaConnection().createStatement();
		s.execute("use " + comparaDb());
		ResultSet rs = s
				.executeQuery("select m.member_id from member as m where m.stable_id="
						+ qw(eAcc) + ";");
		rs.next();
		return rs.getInt(1);
	}

	public static String getCigarLine(int memberId) throws Exception {
		Statement s = EnsemblUtils.comparaConnection().createStatement();
		s.execute("use " + comparaDb());
		ResultSet rs = s.executeQuery("select m.cigar_line from "
				+ "protein_tree_member" + " as m where m.member_id=" + memberId
				+ " and m.cigar_line is not null;");
		rs.next();
		return rs.getString(1);
	}

	public static int getTreeNodeId(int memberId) throws Exception {
		// Return the treeNodeId of the tree that contains this member protein.
		Statement s = EnsemblUtils.comparaConnection().createStatement();
		s.execute("use " + comparaDb());
		String cmd = "select parent.node_id from protein_tree_node parent, protein_tree_node child, protein_tree_member ptm "
				+ "where child.left_index between parent.left_index and parent.right_index "
				+ "and ptm.node_id=child.node_id "
				+ "and parent.parent_id=parent.root_id "
				+ "and ptm.member_id="
				+ memberId + " and ptm.cigar_line is not null;";

		ResultSet rs = s.executeQuery(cmd);
		rs.next();
		return rs.getInt(1);
	}

	public static List<Integer> getRootNodesForMember(int memberId)
			throws Exception {
		Statement s = EnsemblUtils.comparaConnection().createStatement();
		s.execute("use " + comparaDb());

		String cmd = "select node_id from protein_tree_member ptm where node_id="
				+ memberId + ";";
		ResultSet rs = s.executeQuery(cmd);
		int nodeId = rs.getInt(1);

		// rs.last(); // Should have ~7,000 human proteins involved in SLR
		// calcs.
		return new ArrayList<Integer>();
	}

	// GJ 2009-01-19
	public static Iterator<String> getStableIdIteratorForProteinTree(int rootId)
			throws Exception {
		Statement s = EnsemblUtils.comparaConnection().createStatement();
		s.execute("use " + comparaDb());

		// Collect the list of member_ids from this protein tree.
		String cmd = "select m.stable_id from protein_tree_node n, protein_tree_node n2, protein_tree_member ptm, member m "
				+ " where n2.left_index between n.left_index and n.right_index and"
				+ " n.node_id="
				+ rootId
				+ " and ptm.node_id=n2.node_id and m.member_id=ptm.member_id;";
		ResultSet rs = s.executeQuery(cmd);

		Iterator<String> intIterator = SQLUtils.getIterator(rs, 1);
		return intIterator;
	}

	public static ResultSet getHumanComparaPeptides() throws Exception {
		Statement s = EnsemblUtils.comparaConnection().createStatement();
		s.execute("use " + comparaDb());

		// Human taxon = 9606
		String cmd = "select distinct m.stable_id from member m, sitewise_member sm where "
				+ " sm.member_id=m.member_id and m.taxon_id=9606;";
		ResultSet rs = s.executeQuery(cmd);

		rs.last(); // Should have ~7,000 human proteins involved in SLR calcs.
		System.out.println("Number of members found: " + rs.getRow());
		rs.first();
		return rs;
	}
}
