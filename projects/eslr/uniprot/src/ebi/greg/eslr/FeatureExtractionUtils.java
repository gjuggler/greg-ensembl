package ebi.greg.eslr;

import java.sql.ResultSet;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.List;

import uk.ac.ebi.kraken.interfaces.uniprot.UniProtEntry;
import uk.ac.ebi.kraken.interfaces.uniprot.features.Feature;
import uk.ac.ebi.kraken.interfaces.uniprot.features.FeatureType;

public class FeatureExtractionUtils
{

	public static ArrayList<RegionTag> getDomainsFromEnsembl(String stable_id, String latestDb) throws Exception
	{
		Statement s = EnsemblUtils.anonConnection.createStatement();
		s.execute("use " + latestDb + ";");
		
		String query =
				"select t.stable_id,pf.hit_name,pf.seq_start,pf.seq_end,pf.hit_start,pf.hit_end from translation_stable_id t, protein_feature pf "
						+ " where t.stable_id=" + EnsemblUtils.qw(stable_id)
						+ " and pf.translation_id=t.translation_id" + " and pf.hit_name LIKE " + EnsemblUtils.qw("PF%")
						+ " ;";
		ResultSet rs = s.executeQuery(query);
		ArrayList<RegionTag> tags = new ArrayList<RegionTag>();
		while (rs.next())
		{
			tags.add(new RegionTag("PFAM", rs.getString("pf.hit_name"), rs.getInt("pf.seq_start"), rs
					.getInt("pf.seq_end"), rs.getInt("pf.hit_start"), rs.getInt("pf.hit_end")));
		}

		// We won't use Smart for now.
//		query =
//				"select t.stable_id,pf.hit_id,pf.seq_start,pf.seq_end from translation_stable_id t, protein_feature pf "
//						+ " where t.stable_id=" + EnsemblUtils.qw(stable_id)
//						+ " and pf.translation_id=t.translation_id" + " and pf.hit_id LIKE " + EnsemblUtils.qw("SM%")
//						+ " ;";
//		rs = s.executeQuery(query);
//		while (rs.next())
//		{
//			tags.add(new RegionTag("DM_SMART", rs.getString("pf.hit_id"), rs.getInt("pf.seq_start"), rs
//					.getInt("pf.seq_end")));
//		}
		return tags;
	}

	

	public static class RegionTag
	{
		public String tag = null;
		public String value = null;
		public int seq_start = 0;
		public int seq_end = 0;
		public int hit_start = 0;
		public int hit_end = 0;
		
		public RegionTag(String tag, String value, int start, int end, int hit_start, int hit_end)
		{
			this(tag,value,start,end);
			
			this.hit_start = hit_start;
			this.hit_end = hit_end;
		}
		
		public RegionTag(String tag, String value, int start, int end)
		{
			this.tag = tag;
			this.value = value;
			this.seq_start = start;
			this.seq_end = end;	
		}
	}

	public static List<RegionTag> getFeaturesFromUniProt(UniProtEntry uniEntry)
	{
		ArrayList<Feature> features = new ArrayList<Feature>();
		for (FeatureType ft : UniProtUtils.featureTypes)
		{
			features.addAll(uniEntry.getFeatures(ft));
		}
		ArrayList<RegionTag> regionTags = new ArrayList<RegionTag>();
		for (Feature f : features)
		{
			String[] tagval = UniProtUtils.featureToTagValue(f);
			RegionTag rt = new RegionTag(tagval[0],tagval[1],f.getFeatureLocation().getStart(),f.getFeatureLocation().getEnd());
			regionTags.add(rt);
		}
		return regionTags;
	}

	public static List<RegionTag> getFeaturesFromPdb(String pdbAcc,String pdbDssp, String pdbAccess)
	{
		ArrayList<RegionTag> regions = new ArrayList<RegionTag>();
		
		for (int i=0; i < pdbDssp.length(); i++)
		{
			RegionTag rt = new RegionTag("DSSP",""+pdbDssp.charAt(i),i,i);
			regions.add(rt);
		}
		for (int i=0; i < pdbAccess.length(); i++)
		{
			RegionTag rt = new RegionTag("ACC",""+pdbAccess.charAt(i),i,i);
			regions.add(rt);
		}
		return regions;
	}

}
