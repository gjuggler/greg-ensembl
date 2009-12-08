package ebi.greg.eslr;

public class AlignmentResult
{
	public String labelA;
	public String labelB;
	
	public String seqA;
	public String seqB;
	
	double score;
	
	public AlignmentResult(String labelA, String labelB, String seqA, String seqB, double score)
	{
		this.labelA = labelA;
		this.labelB = labelB;
		this.seqA = seqA;
		this.seqB = seqB;
		this.score = score;
	}
	
	public double getScore()
	{
		return score;
	}
	
	public String getName1()
	{
		return labelA;
	}
	
	public String getName2()
	{
		return labelB;
	}
	
	public char[] getSequence1()
	{
		return seqA.toCharArray();
	}
	
	public char[] getSequence2()
	{
		return seqB.toCharArray();
	}
	
	public String getOriginalSequence1()
	{
		return seqA.replaceAll("-", "");
	}
	
	public String getOriginalSequence2()
	{
		return seqB.replaceAll("-", "");
	}
}
