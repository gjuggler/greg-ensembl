package ebi.greg.eslr;

import java.io.File;
import java.io.FileInputStream;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.Statement;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class GregDBUtils
{
//	static String WORKING_DB = "";
//	static String GREG_DB = "";

	static Connection connection;
	static
	{
		if (Main.DEBUG)
			System.err.println("Starting GregDB...");
		try
		{
			connection = EnsemblUtils.comparaConnection();
		} catch (Exception e)
		{
			e.printStackTrace();
			System.err.println("Could not connect to the database. Quitting...");
			System.exit(0);
		}
		if (Main.DEBUG)
			System.err.println("Done!");
	}

	/**
	 * Creates the residue-feature table.
	 */
	public static void createResidueFeatureTable() throws Exception
	{
		System.out.println("Creating tables if necessary...");
		Statement s = connection.createStatement();
		/*
		 * Execute the SQL commands stored in a file. Easier than dealing with Java strings.
		 */
		s.executeUpdate("use " + outputDb() + ";");
		String sql = JavaUtils.streamToString(new FileInputStream(new File("sitewise_aln_tag.sql")));
		String[] statements = sql.split(";");
		for (String cmd : statements)
		{
			if (cmd.trim().length() == 0)
				continue;
			String newCmd = cmd.replaceAll("aln_tag", uniprotOutputTableName());
			s.executeUpdate(newCmd);
			newCmd = cmd.replaceAll("aln_tag", domainOutputTableName());
			s.executeUpdate(newCmd);
			newCmd = cmd.replaceAll("aln_tag", pdbOutputTableName());
			s.executeUpdate(newCmd);
		}
		
		
	}

	public static String workingDb()
	{
		return Config.inputDb;
	}
	
	public static String outputDb()
	{
		return "";
	}
	
	public static String uniprotOutputTableName()
	{
		return "";
	}
	
	public static String domainOutputTableName()
	{
		return "";
	}

	public static String pdbOutputTableName()
	{
		return "";
	}
		
	private static String getTagValueRow(int memberId, int treeNodeId, String uniProtId, int ensemblPos,
			int sitewisePos, int sitewiseId, String columnTag, String columnValue) throws Exception
	{
		// Returns the SQL command that will store the given tag value.
		//		int sitewiseId = EnsemblUtils.getSitewiseId(treeNodeId, sitewisePos);
		//		System.out.println(treeNodeId + "  " + sitewiseId);
		String eq = "=";
		String c = ", ";

		String cmd = memberId + c;
		cmd += qw(uniProtId) + c;
		cmd += ensemblPos + c;
		cmd += sitewiseId + c;
		cmd += qw(columnTag) + c;
		cmd += qw(columnValue);

		cmd = "(" + cmd + ")";
		return cmd;
	}

	static String qw(String s)
	{
		return "\"" + s + "\"";
	}

	static String qqw(String s)
	{
		return "`" + s + "`";
	}
	
	public static class PosTagValue
	{
		public int aln_position;
		public int feature_position;
		public int member_id;
		public String stable_id;
		public int node_id;
		public String source;
		public String tag;
		public String value;

		public PosTagValue(int node_id, int aln_position, int member_id, String stable_id, String source, String tag, String value)
		{
			this.node_id = node_id;
			this.member_id = member_id;
			this.stable_id = stable_id;
			this.source = source;
			this.tag = tag;
			this.value = value;
			this.aln_position = aln_position;
			this.feature_position = 0;
		}

		public PosTagValue(int node_id, int aln_position, int member_id, String stable_id, String source, String tag, String value, int feature_position)
		{
			this(node_id,aln_position,member_id,stable_id,source,tag,value);
			
			this.feature_position = feature_position;
		}
		
		public static String getInsertHeader()
		{
			return " (node_id,aln_position,member_id,source,tag,value,feature_position) ";
		}
		
		public String toString()
		{
			return getInsertRow();
		}
		
		public String getInsertRow()
		{
			String c = ",";
			String[] cols = new String[]{node_id+"",aln_position+"",member_id+"",qw(source),qw(tag),qw(value),feature_position+""};
			return "(" + JavaUtils.join(",",cols) + ")";
		}
		
		public static String getInsertRows(List<PosTagValue> tagvals)
		{
			return JavaUtils.join(tagvals, ",");
		}
	}

	
	
}
