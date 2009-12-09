package ebi.greg.eslr;

import java.io.File;
import java.io.FileInputStream;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class SQLUtils {

	public static Map<String, String> getConnectionMapFromString(String url) {
		Pattern p = Pattern.compile("mysql://(.*?):(.*?)@(.*?)/(.*)");
		Matcher m = p.matcher(url);
		HashMap<String, String> map = new HashMap<String, String>();
		if (m.find()) {
			String user = m.group(1);
			String password = m.group(2);
			String host = m.group(3);
			String port = "";
			if (host.contains(":")) {
				String[] toks = host.split(":");
				host = toks[0];
				port = toks[1];
			}
			String db = m.group(4);
			map.put("user", user);
			map.put("password", password);
			map.put("host", host);
			map.put("port", port);
			map.put("db", db);
		}
		return map;
	}

	public static String getDBFromString(String url) {
		Map<String, String> map = getConnectionMapFromString(url);
		return map.get("db");
	}

	public static Connection getConnectionFromString(String url)
			throws Exception {
		Map<String, String> map = getConnectionMapFromString(url);
		String user = map.get("user");
		String password = map.get("password");
		String host = map.get("host");
		String port = map.get("port");
		String db = map.get("db");
		Connection c = null;
		Class.forName("com.mysql.jdbc.Driver").newInstance();
		if (port.equals("")) {
			c = DriverManager.getConnection("jdbc:mysql://" + host + "/" + db,
					user, password);
		} else {
			c = DriverManager.getConnection("jdbc:mysql://" + host + ":" + port
					+ "/" + db, user, password);
		}
		return c;
	}

	public static Connection getConnectionFromConfigFile(String filename)
			throws Exception {

		System.setProperty("http.proxyHost", "wwwcache.sanger.ac.uk");
		System.setProperty("http.proxyPort", "3128");
		File f = new File(filename);
		String s = JavaUtils.streamToString(new FileInputStream(f));
		String[] lines = s.split("\n");
		HashMap<String, String> props = new HashMap<String, String>();
		String host = "";
		String user = "";
		String password = "";
		String port = "";
		String db = "";
		for (String line : lines) {
			if (line.startsWith("#"))
				continue;
			String[] keyval = line.split("=");
			if (keyval.length < 2)
				continue;

			String k = keyval[0].trim().toLowerCase();
			String v = keyval[1].trim();

			if (k.equals("host")) {
				host = v;
			} else if (k.equals("user")) {
				user = v;
			} else if (k.equals("password")) {
				password = v;
			} else if (k.equals("port")) {
				port = v;
			} else if (k.equals("db")) {
				db = v;
			}
		}
		Connection c = null;
		Class.forName("com.mysql.jdbc.Driver").newInstance();

		System.out.println("User: " + user + "\n" + "Host: " + host + "\nDB: "
				+ db + "\nPW:" + password + "\n");

		if (port.equals("")) {
			c = DriverManager.getConnection("jdbc:mysql://" + host + "/" + db,
					user, password);
		} else {
			c = DriverManager.getConnection("jdbc:mysql://" + host + ":" + port
					+ "/" + db, user, password);
		}
		return c;
	}

	public static int rowCount(ResultSet rs) throws Exception {
		int oldRow = rs.getRow();
		rs.last();
		int count = rs.getRow();
		if (oldRow == 0)
			rs.beforeFirst();
		else
			rs.absolute(oldRow);
		return count;
	}

	static String rowString(ResultSet rs) {
		String row = "";
		int i = 1;
		while (true) {
			try {
				String s = rs.getString(i);
				row += s + "\t";
				i++;
			} catch (Exception e) {
				break;
			}
		}
		return row;
	}

	/**
	 * Returns an iterator wrapped around a ResultSet.
	 * 
	 * @param rs
	 * @param colIndex
	 * @return
	 */
	public static Iterator<String> getIterator(final ResultSet rs, int colIndex) {
		if (colIndex == 0)
			colIndex = 1;

		final int ind = colIndex;

		return new Iterator<String>() {
			private String next;

			public boolean hasNext() {
				if (next == null) {
					try {
						if (!rs.next()) {
							return false;
						}
						next = rs.getString(ind);
					} catch (Exception e) {
						e.printStackTrace();
						return false;
					}
				}
				return true;
			}

			public String next() {
				if (!hasNext()) {
					throw new NoSuchElementException();
				}
				String retval = next;
				next = null;
				return retval;
			}

			public void remove() {
				throw new UnsupportedOperationException("no remove allowed");
			}
		};
	}

	public static ArrayList<String> getColumnStrings(ResultSet rs, int col)
			throws Exception {
		ArrayList<String> strings = (ArrayList<String>) getColumn(rs, col);
		return strings;
	}

	public static ArrayList<Float> getColumnFloats(ResultSet rs, int col)
			throws Exception {
		ArrayList<Float> floats = (ArrayList<Float>) getColumn(rs, col);
		return floats;
	}

	public static ArrayList<? extends Object> getColumn(ResultSet rs, int col)
			throws Exception {
		ArrayList<Object> list = new ArrayList<Object>();
		while (rs.next()) {
			Object o = rs.getObject(col);
			list.add(o);
		}
		return list;
	}
}
