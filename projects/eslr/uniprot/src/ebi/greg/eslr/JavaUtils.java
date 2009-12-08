package ebi.greg.eslr;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.lang.reflect.Field;
import java.net.URL;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Properties;

public class JavaUtils
{

	public static void loadPropertiesIntoClass(String filename, Class c) throws Exception
	{
		Iterator<String> it = getFileLineIterator(filename);
		while (it.hasNext())
		{
			String s = it.next();
			String[] propval = s.split("=", 2);
			if (propval.length < 2)
				continue;

			String prop = propval[0].trim();
			String val = propval[1].trim();

			try
			{
				Field f = c.getField(prop);
				f.set(null, val);
			} catch (Exception e)
			{
				continue;
			}
		}
	}

	public static String urlToString(String s) throws Exception
	{
		URL url = new URL(s);
		InputStream in = url.openStream();
		return streamToString(in);
	}

	public static Iterator<String> getFileLineIterator(String filename) throws Exception
	{
		File f = new File(filename);
		final BufferedReader in = new BufferedReader(new FileReader(f));
		return new Iterator<String>()
		{
			String nextLine = null;

			public boolean hasNext()
			{
				try
				{
					if (nextLine == null)
						nextLine = in.readLine();
					if (nextLine != null)
						return true;
				} catch (Exception e)
				{
					return false;
				}
				return false;
			}

			public String next()
			{
				try
				{
					if (nextLine == null)
						nextLine = in.readLine();
					String ln = nextLine;
					nextLine = null;
					return ln;
				} catch (Exception e)
				{
					throw new RuntimeException(e.getMessage());
				}
			}

			public void remove()
			{
			}

		};
	}

	// UNIX
	public Properties getEnvironment() throws java.io.IOException
	{
		Properties env = new Properties();
		env.load(Runtime.getRuntime().exec("env").getInputStream());
		return env;
	}

	/**
	 * Here's a nifty one -- Prints out a list of field names and values for a
	 * given class.
	 * 
	 * Handy for collecting and reporting error counts while maintaining Java
	 * refactoring abilities!
	 * 
	 * @param c
	 * @param o
	 */
	public static String getClassFields(Class c, Object o)
	{
		// NOTE: the Object may be null if we're dealing with a static class.
		Field[] fields = c.getFields();
		StringBuffer sb = new StringBuffer();
		sb.append(c.getName() + ":" + "\n");
		for (Field f : fields)
		{
			try
			{
				String name = f.getName();
				Object val = f.get(o);
				if (val instanceof ArrayList)
					val = getArrayListString((ArrayList) val);
				sb.append(name + "\t" + val + "\n");
			} catch (Exception e)
			{
				continue;
			}
		}
		return sb.toString();
	}

	private static String getArrayListString(ArrayList a)
	{
		StringBuffer sb = new StringBuffer();
		for (Object o : a)
		{
			sb.append(o + "\t");
		}
		return sb.toString();
	}

	public static String streamToString(InputStream in) throws Exception
	{
		StringBuffer b = new StringBuffer();
		String s;
		BufferedReader r = new BufferedReader(new InputStreamReader(in));
		while ((s = r.readLine()) != null)
		{
			b.append(s);
			b.append("\n");
		}
		return b.toString();
	}

	/**
	 * Retrieves an accession list from a file.
	 * 
	 * @param f
	 * @return
	 * @throws Exception
	 */
	public static ArrayList<String> streamToArray(InputStream in)
	{
		ArrayList<String> strings = new ArrayList<String>();
		try
		{
			BufferedReader r = new BufferedReader(new InputStreamReader(in));
			String line;
			while ((line = r.readLine()) != null)
			{
				strings.add(line);
			}
		} catch (Exception e)
		{
			e.printStackTrace();
		}
		return strings;
	}

	public static String join(Collection s, String delimiter)
	{
		StringBuffer buffer = new StringBuffer();
		Iterator iter = s.iterator();
		while (iter.hasNext())
		{
			buffer.append(iter.next().toString());
			if (iter.hasNext())
			{
				buffer.append(delimiter);
			}
		}
		return buffer.toString();
	}

	public static String join(String delimiter, String... array)
	{
		StringBuffer buffer = new StringBuffer();
		int len = array.length;
		for (int i = 0; i < len; i++)
		{
			buffer.append(array[i]);
			if (i < len - 1)
				buffer.append(delimiter);
		}
		return buffer.toString();
	}

	public static int[] toInts(List<Integer> ints)
	{
		int[] intArray = new int[ints.size()];
		for (int i = 0; i < ints.size(); i++)
		{
			intArray[i] = ints.get(i).intValue();
		}
		return intArray;
	}

	public static float[] tofloats(List<Float> floats)
	{
		float[] floatArray = new float[floats.size()];
		for (int i = 0; i < floats.size(); i++)
		{
			floatArray[i] = floats.get(i).floatValue();
		}
		return floatArray;
	}
}
