package ebi.greg.eslr;

import java.io.File;
import java.lang.reflect.Method;
import java.util.ArrayList;

import uk.ac.ebi.kraken.interfaces.common.Sequence;
import uk.ac.ebi.kraken.interfaces.uniprot.UniProtEntry;
import uk.ac.ebi.kraken.interfaces.uniprot.features.Feature;
import uk.ac.ebi.kraken.interfaces.uniprot.features.FeatureDescription;
import uk.ac.ebi.kraken.interfaces.uniprot.features.FeatureType;
import uk.ac.ebi.kraken.interfaces.uniprot.features.HelixFeature;
import uk.ac.ebi.kraken.uuw.services.remoting.EntryIterator;
import uk.ac.ebi.kraken.uuw.services.remoting.EntryRetrievalService;
import uk.ac.ebi.kraken.uuw.services.remoting.Query;
import uk.ac.ebi.kraken.uuw.services.remoting.UniProtQueryBuilder;
import uk.ac.ebi.kraken.uuw.services.remoting.UniProtQueryService;
import uk.ac.ebi.kraken.uuw.services.remoting.UniProtRemoteServiceFactory;

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
	static FeatureType[] featureTypes = new FeatureType[] {
			FeatureType.ACT_SITE, 
			FeatureType.DNA_BIND, 
			FeatureType.METAL, // Metal-binding
			FeatureType.TRANSMEM, // Transmembrane domain
			FeatureType.CARBOHYD, // Glycosylated
			FeatureType.MOD_RES, // Other modifications
			FeatureType.DISULFID,
			FeatureType.BINDING,
			FeatureType.CROSSLNK,
			FeatureType.SIGNAL
	};

	static UniProtRemoteServiceFactory factory;
	static {
		// Try to set up a proxy using uniprot.properties config. See:
		// http://www.ebi.ac.uk/uniprot/remotingAPI/faq.html#p
		try {
			new File("uniprot.properties").createNewFile();
			factory = new UniProtRemoteServiceFactory(new File(
					"uniprot.properties"));
		} catch (Exception e) {
			e.printStackTrace();
			factory = new UniProtRemoteServiceFactory();
		}
	}
	static EntryRetrievalService ers = factory.getEntryRetrievalService();
	static UniProtQueryService uqs = factory.getUniProtQueryService();

	public static final String[] featureToTagValue(Feature f) {
		String tag = f.getType().getName();
		String value = "";
		String s;
		switch (f.getType()) {
		default:
			try {
				Method m = f.getClass().getMethod("getFeatureDescription", new Class[]{});
				FeatureDescription fd = (FeatureDescription) m.invoke(f, new Object[]{});
				value = fd.getValue().toLowerCase();
			} catch (Exception e) {
//				e.printStackTrace();
				value = "NA";
			}
			
			break;
		}
		return new String[] { tag, value };
	}

	/**
	 * Returns the entries for a given arraylist of UniProt accessions.
	 * 
	 * @param accs
	 * @return
	 */
	public static EntryIterator<UniProtEntry> getEntries(ArrayList<String> accs) {
		Query query = UniProtQueryBuilder.buildIDListQuery(accs);

		EntryIterator<UniProtEntry> entries = uqs.getEntryIterator(query);
		return entries;
	}

    public static Map<String,Object> getEntryHash(String acc) {
	

    }

	public static UniProtEntry getEntry(String acc) {
	    
		UniProtEntry entry = ers.getUniProtEntry(acc);
		if (entry == null) {
			Query query = UniProtQueryBuilder.buildQuery(acc);
			EntryIterator<UniProtEntry> entries = uqs.getEntryIterator(query);
			if (entries.hasNext())
				entry = entries.next();
		}
		return entry;
	}

	public static String getResidue(UniProtEntry entry, int index) {
		Sequence sub = entry.getSequence().subSequence(index, index + 1);
		return sub.getValue();
	}

	public static String getAccession(UniProtEntry entry) {
		return entry.getPrimaryUniProtAccession().getValue();
	}

}
