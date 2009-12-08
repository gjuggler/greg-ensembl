package ebi.greg.eslr;

import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.Writer;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;

public class LogUtils
{
	private static File logFile;
	private static Writer logWriter;
	static final String SEP = File.separator;

	private static void init()
	{
		if (logFile == null)
		{
			File logDir = new File("logs");
			logDir.mkdirs();

			logFile = new File(logDir.getAbsolutePath() + SEP + dateStamp() + ".txt");
			try
			{
				logWriter = new FileWriter(logFile);
			} catch (Exception e)
			{

			}
		}
	}

	private static String dateStamp()
	{
		Date now = new Date();
		DateFormat df = new SimpleDateFormat("MM-dd HH'h'mm'm'ss's'");
		return df.format(now);
	}
	
	private static String timeStamp()
	{
		Date now = new Date();
		DateFormat df = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss.S");
		return df.format(now);
	}
	
	public static void log(String s)
	{
		if (true)
			return;
		init();
		try
		{
			StackTraceElement[] st = Thread.currentThread().getStackTrace();
			StackTraceElement ste = st[st.length-1];
			String classMethod = ste.getClassName()+"."+ste.getMethodName();
			logWriter.write("###### "+timeStamp()+" "+classMethod+"\n");
			logWriter.write(s + "\n");
			logWriter.flush();
		} catch (Exception e)
		{
			e.printStackTrace();
		}
	}
}
