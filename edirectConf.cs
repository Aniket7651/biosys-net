using IronPython.Hosting;
using Microsoft.Scripting.Hosting;
using System;
using System.Diagnostics;

namespace BioSySNet
{
    public class edirectConf
    {
        void Configuration(){}
    }
}


namespace PyConfiguration
{
    /// <summary>
    /// Related to all python configurations and setup programs
    /// </summary>
    public class Config
    {
        Process process = new Process();
        string conf_pyC3 = AppDomain.CurrentDomain.BaseDirectory + @"pyScript\\pyC3.conf";

        public void ConfPython(string pythonPath)
        {
            // string configure_pyC3 = AppDomain.CurrentDomain.BaseDirectory + @"pyScript\\pyC3.conf";
            StreamWriter writer = new(conf_pyC3, false);
            writer.WriteLine($"ENVPATH=\"{pythonPath}\"");
            writer.Close();
        }

        private string ReadENV()
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Install pip requirements for executing graphs using python interpereter
        /// </summary>
        public void RequiredModules()  // no module name 'os' ERROR ........
        {
            string path = AppDomain.CurrentDomain.BaseDirectory + "pyscript\\requirements.txt";
            process.StartInfo = new ProcessStartInfo($"python -m pip install -r {path}")
            {
                UseShellExecute = false,
                CreateNoWindow = true,
                RedirectStandardOutput = true

            }; process.Start();
            string output = process.StandardOutput.ReadToEnd();
            process.WaitForExit(); process.Close(); Console.WriteLine(output);
        }

        // local function for making connection between pyscript folder located in bin\Debug\net6.0
        // with parameters of the python script and second one, the function build in perticular script
        public dynamic Make_con_pyScript(string pyscript, string funcName)
        {
            string path = AppDomain.CurrentDomain.BaseDirectory + "pyscript\\" + pyscript;
            ScriptEngine engine = Python.CreateEngine();
            ScriptScope scope = engine.CreateScope();
            engine.ExecuteFile(path, scope);
            dynamic funcSelector = scope.GetVariable(funcName);
            return funcSelector;
        }
    }
}