package ebi.greg.eslr;

import java.net.URL;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import uk.ac.ebi.kraken.interfaces.uniprot.features.FeatureType;
import ebi.greg.eslr.FeatureExtractionUtils.RegionTag;

/**
 * Extracts protein features from the UniProt database, given a list of Ensembl
 * IDs.
 * 
 * Some information on the feature types within UniProt:
 * http://www.ebi.ac.uk/uniprot
 * /remotingAPI/javadoc/domain/api/uk/ac/ebi/kraken/interfaces
 * /uniprot/features/FeatureType.html
 * 
 * UniProt example file: http://www.genome.jp/dbget-bin/www_bget?uniprot+P78536
 * UniProt user manual: http://expasy.org/sprot/userman.html (especially note
 * the description of turning dssp sec-struct results into simple 3-character
 * factors)
 * 
 * @author Greg
 * 
 */
public class UniProtUtils {
	static String[] featureTypes = new String[] {
		"ACT_SITE",
		"DNA_BIND",
		"METAL",
		"TRANSMEM",
		"CARBOHYD",
		"MOD_RES",
		"DISULFID",
		"BINDING",
		"CROSSLINK",
		"SIGNAL"
	};
	
	static Map<String,String> entryHash;
	public static Map<String,String> getEntryHash() {
		if (entryHash == null)
			entryHash = new HashMap<String,String>();
		return entryHash;
	}
	
	public static String getEntryString(String acc) {
		if (!getEntryHash().containsKey(acc)) {
			// Fetch the entry string online.
			try {
				URL url = new URL("http://www.uniprot.org/uniprot/" + acc + ".txt");
				String s = JavaUtils.streamToString(url.openStream());
				getEntryHash().put(acc, s);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		return getEntryHash().get(acc);
	}
	
	public static String getSequence(String acc) {
		String entry = getEntryString(acc);
		String seq = null;
		
		Pattern p = Pattern.compile("SQ.*?\\n(.*)\\/\\/", Pattern.DOTALL);
		Matcher m = p.matcher(entry);
		if (m.find()) {
			seq = m.group(1);
		}
		String[] lines = seq.split("\n");
		seq = JavaUtils.join("", lines).replaceAll("\\s", "");
//		System.out.println(seq);
		return seq;
	}
	
	public static List<String> getPDBCrossReferences(String acc) {
		String entry = getEntryString(acc);
		List<String> pdbRefs = new ArrayList<String>();
		
		Pattern p = Pattern.compile("DR\\s+PDB;\\s*(\\S+);");
		Matcher m = p.matcher(entry);
		while (m.find()) {
			String pdbEntry = m.group(1);
			pdbRefs.add(pdbEntry);
		}
		return pdbRefs;
	}
	
	public static List<RegionTag> getUniProtFeatures(String featureType, String acc) {
		ArrayList<RegionTag> tags = new ArrayList<RegionTag>();
		
		String entry = getEntryString(acc);
		Pattern p = Pattern.compile("FT\\s{3}("+featureType+".*?)\\n()\\S\\S\\s{3}\\S+",Pattern.DOTALL);
		Matcher m = p.matcher(entry);
		int index = 0;
		while (m.find(index)) {
			String s = m.group(1);
			index = m.start(2);
			s = s.replaceAll("FT\\s", "");
			s = s.replaceAll("\\n\\s*"," ");
			String[] toks = s.split("\\s+");
			String type = toks[0];
			String start = toks[1];
			String end = toks[2];
			String desc = "";
			for (int i=3; i < toks.length; i++) {
				desc += " "+toks[i];
			}
			RegionTag rt = new RegionTag(type,desc,Integer.parseInt(start),Integer.parseInt(end));
			tags.add(rt);
		}
		return tags;
	}


	public static String getResidue(String acc, int index) {
		return getSequence(acc).charAt(index+1)+"";
	}


}
