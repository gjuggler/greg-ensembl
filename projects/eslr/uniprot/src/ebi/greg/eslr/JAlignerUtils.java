package ebi.greg.eslr;

import jaligner.Alignment;
import jaligner.NeedlemanWunschGotoh;
import jaligner.Sequence;
import jaligner.SmithWatermanGotoh;
import jaligner.matrix.Matrix;
import jaligner.matrix.MatrixLoader;

import java.io.FileInputStream;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Date;

public class JAlignerUtils
{

	public static AlignmentResult alignProteins(String l1, String l2, String s1, String s2)
	{
		try
		{
			if (s1.length() > 8000 || s2.length() > 8000)
				throw new RuntimeException("Sequence longer than 8000!");
			Matrix matrix = MatrixLoader.load("BLOSUM90");
			Sequence seq1 = new Sequence(s1, l1, l1, Sequence.PROTEIN);
			Sequence seq2 = new Sequence(s2, l2, l2, Sequence.PROTEIN);
//			Alignment aln = SmithWatermanGotoh.align(seq1,seq2,matrix,10,1); // Local alignment.
			Alignment aln = NeedlemanWunschGotoh.align(seq1, seq2, matrix, 10, 1); // Global alignment.
			return new AlignmentResult(l1,l2,new String(aln.getSequence1()),new String(aln.getSequence2()),aln.getScore());
		} catch (Exception e)
		{
			e.printStackTrace();
			return null;
		}
	}

	/**
	 * Maps sequence A to sequence B, returning an int[] showing the
	 * corresponding index (0-based) of position A[i] in sequence B. Some
	 * indices may contain a value of -1, indicating that this position
	 * corresponds to a gap in B.
	 * 
	 * @param a
	 * @return
	 */
	static int[] mapSequences(AlignmentResult aln, boolean anchorA)
	{
		char[] a, b;
		int mapLength;
		if (anchorA)
		{
			a = aln.getSequence1();
			b = aln.getSequence2();
			mapLength = aln.getOriginalSequence1().length();
		} else
		{
			a = aln.getSequence2();
			b = aln.getSequence1();
			mapLength = aln.getOriginalSequence2().length();
		}

		// Go through the alignment columns, spitting out the appropriate output as we go.
		// NB: we're going to make this mapping 0-based, so the first position in both sequences is 0.
//		int[] mappedIndices = new int[aln.getOriginalSequence1().length()];
		int[] mappedIndices = new int[mapLength];
		int aInd = 0;
		int bInd = 0;
		for (int i = 0; i < a.length; i++)
		{
			char aC = a[i];
			char bC = b[i];
			if (aC == Alignment.GAP)
			{
				bInd++;
			} else if (bC == Alignment.GAP)
			{
				mappedIndices[aInd] = -1;
				aInd++;
			} else
			{
				mappedIndices[aInd] = bInd;
				aInd++;
				bInd++;
			}
		}
		return mappedIndices;
	}

	static public class AlignmentWithMaps
	{
		public int[] aToB;
		public int[] bToA;
		public AlignmentResult aln;

		public AlignmentWithMaps(AlignmentResult aln)
		{
			this.aln = aln;
			aToB = mapSequences(aln, true);
			bToA = mapSequences(aln, false);
		}
	}

}
