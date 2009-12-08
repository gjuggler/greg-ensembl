package ebi.greg.eslr;

import java.io.File;
import java.io.FileInputStream;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;

import uk.ac.ebi.kraken.interfaces.uniprot.UniProtEntry;
import ebi.greg.eslr.JAlignerUtils.AlignmentWithMaps;

public class Main {
	public static final boolean DEBUG = true;

	public static void main(String[] args) throws Exception {
		String ensAcc = "ENSP00000369497";
		String ensUrl = "mysql://ensadmin:ensembl@127.0.0.1:5425/gj1_57";
		if (args.length >= 1)
			ensAcc = args[0];
		if (args.length >= 2)
			ensUrl = args[1];

		Config.comparaUrl = ensUrl;
		Config.inputDb = SQLUtils.getDBFromString(ensUrl);

		EnsemblUtils.anonConnection();
		EnsemblFeatureExtractor extractor = new EnsemblFeatureExtractor(ensAcc);
	}

	public static float outputAlign(String base, String aligned,
			AlignmentWithMaps alignToBase) {
		System.out.print("uniprot: ");
		int coverCount = 0;
		for (int i = 0; i < base.length(); i++) {
			int ba = alignToBase.bToA[i];
			if (ba >= 0) {
				char uniChar = aligned.charAt(ba);
				System.out.print(uniChar);
				coverCount++;
			} else
				System.out.print("-");
		}
		System.out.println();
		return ((float) coverCount / base.length());
	}
}
