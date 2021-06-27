using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.Serialization;

namespace BasicBFB
{    

    [Serializable()]
    public class NotConvergedException : Exception
    {
        private double value1;
        private double value2;

        protected NotConvergedException()
           : base()
        { }

        public NotConvergedException(double value1, double value2) :
           base(String.Format("Not converged after MAXITER.  " +
               "                Last bracket: {0:g3}, {1:g3}", value1, value2))
        {
            this.value1 = value1;
            this.value2 = value2;
        }

        public NotConvergedException(double value1, double value2, string message)
           : base(message)
        {
            this.value1 = value1;
            this.value2 = value2;
        }

        public NotConvergedException(double value, string message, Exception innerException) :
           base(message, innerException)
        {
            value1 = value;
        }   

        protected NotConvergedException(SerializationInfo info,
                                    StreamingContext context)
           : base(info, context)
        { }

        public double bracket1
        { 
            get { return this.value1; } 
        }

        public double bracket2
		{
            get { return this.value2; }
		}

        public (double, double) lastBracket
		{
            get { return (value1, value2); }
		}
    }


    // ==================================================================================



    [Serializable()]
    public class UnexpectedValueException : Exception
    {
        private double value1;
        private string objName;

        protected UnexpectedValueException()
           : base()
        { }

        public UnexpectedValueException(string name, double value1) :
           base(String.Format("Unexpected value: {0} = {1g3} falls outside model expectations", name, value1))
        {
            this.value1 = value1;
            this.objName = name;
        }

        public UnexpectedValueException(string name, double value1, string message)
           : base(message)
        {
            this.value1 = value1;
            this.objName = name;
        }

        public UnexpectedValueException(int value, string message, Exception innerException) :
           base(message, innerException)
        {
            value1 = value;
        }

        protected UnexpectedValueException(SerializationInfo info,
                                    StreamingContext context)
           : base(info, context)
        { }
    }


}
