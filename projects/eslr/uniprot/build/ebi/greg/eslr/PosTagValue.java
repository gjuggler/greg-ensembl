package ebi.greg.eslr;

import java.util.List;

public class PosTagValue {
	public int aln_position;
	public int feature_position;
	public String stable_id;
	public String source;
	public String tag;
	public String value;
	public char residue;

	public PosTagValue(int aln_position, String stable_id, String source,
			String tag, String value, char residue) {
		this.stable_id = stable_id;
		this.source = source;
		this.tag = tag;
		this.value = value;
		this.aln_position = aln_position;
		this.feature_position = 0;
		this.residue = residue;
	}

	public PosTagValue(int aln_position, String stable_id, String source,
			String tag, String value, char residue, int feature_position) {
		this(aln_position, stable_id, source, tag, value, residue);

		this.feature_position = feature_position;
	}

	public String toString() {
		return getRow();
	}

	public String getRow() {
		String[] cols = new String[] { stable_id, aln_position + "", source,
				tag, value, residue + "" };
		return JavaUtils.join("\t", cols);
	}

	static String qqw(String s) {
		return "`" + s + "`";
	}

	static String qw(String s) {
		return "\"" + s + "\"";
	}

	public static String tsvHeader() {
		return JavaUtils.join("\t", new String[] { "stable_id", "aln_position",
				"source", "tag", "value","residue" });
	}

	public static String toTSV(List<PosTagValue> tagvals) {
		return JavaUtils.join(tagvals, "\n");
	}

	public static String getInsertRows(List<PosTagValue> tagvals) {
		return JavaUtils.join(tagvals, ",");
	}
}
