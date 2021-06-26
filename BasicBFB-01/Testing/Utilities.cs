using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace BasicBFB.Testing
{
	public static class Utilities
	{
		// Save a string to a new file, or overwrite an existing one
		public static void saveTextAs(string txt, string path)
		{
			File.WriteAllText(path, txt);
		}



	}
}
