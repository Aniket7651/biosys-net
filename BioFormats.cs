/*
~BioFormats.cs 
file contains programs related to the biological file formats such as,
FASTA, FASTQ, PDB etc.. also the programs that's contains in this file can perform
on creating alignment file as *.sysalign, csv file from FASTA and FASTQ file for making dataset

@ Author: ANIKET YADAV (aniketyadav8687@gmail.com)
*/

using System;
using System.IO.Compression;
using System.Net;
using System.Diagnostics;
using IronPython.Compiler.Ast;
using IronPython.Runtime.Operations;
using HtmlAgilityPack;
using System.Reflection.Metadata;


namespace BioSySNet
{
    #pragma warning disable

    /// <summary>
    /// class <c>BioFormat</c> contain programs related to organization and manipulation of biological
    /// file formats, ie. <i>PDB</i>, <i>FASTA</i>, <i>FASTQ</i>, etc..
    /// </summary>
    public class BioFormats
    {
        /// <summary>
        /// Function for automatic selection of format
        /// </summary>
        /// <param name="path">Full path of any file</param>
        /// <returns>Name of the .format as string</returns>
        public string AutoFormatSelector(string path) {
            string[] form = path.Split('.');
            return form.Last();
        }

        /// <summary> structure for multiple value return type from PDB format</summary>
        public struct PDBvars
        {
            public string? title, header;
            public List<List<string>>? atoms;
        }

        /// <summary>
        /// Function can read PDB (Protein Data Bank) file format
        /// </summary>
        /// <param name="path_pdbFile">Full path of the <b>PDB</b> file</param>
        /// <returns>
        /// Multiple values like; <c>.title</c> of the PDB, <c>.atoms</c> of the Amino acid as List Of List datafram,
        /// <c>.header</c> of the file
        /// </returns>
        public PDBvars PDBReader(string path_pdbFile)
        {
            List<List<string>> atom = new();
            string title = ""; string head = "";

            static string TagsToString(string line, string T)
            {
                string removedline = line.Replace("\n", " ");
                string content = removedline.Replace(T, "");
                return content;
            }

            string[] lines = File.ReadAllLines(path_pdbFile);
            foreach (string line in lines)
            {
                string tag = line.Split(' ')[0];
                if (tag.Contains("ATOM"))
                {
                    List<string> lineLis = new();
                    foreach (string item in line.Split(" "))
                    {
                        if (item != "") {
                            lineLis.Add(item);
                            Console.Write(item + '\t');
                        }
                    } atom.Add(lineLis);
                    Console.Write('\n');
                }
                if (tag.Contains("HEADER"))
                {
                    head += TagsToString(line, "HEADER");
                }
                if (tag.Contains("TITLE"))
                {
                    title += TagsToString(line, "TITLE");
                }
            }
            var pdbvars = new PDBvars() { atoms = atom, header = head, title = title };
            return pdbvars;
        }

        /// <summary> structure store multiple return type value </summary>
        public class FASTA {
            public string? Header, Seq;
            public string[]? MultiHeader, MultiSeq;
            public int Len, readLen;
        }

        /// <summary>
        /// A Temporary FASTA format reader for <b>only single sequence in a file</b>
        /// </summary>
        /// <param name="path">full path of the <b>FASTA/TXT</b> file</param>
        /// <returns>multi variable return like; <c>.Header</c>, <c>.Seq</c> and <c>.Len</c></returns>
        public FASTA Readfasta(string path) {
            string[] reads = File.ReadAllLines(path);
            string seq = "";
            for (int i = 1; i <= reads.Length - 1; i++) {
                string replacedNewLine = reads[i].Replace("\n", "");
                seq += replacedNewLine;
            }
            var fastaInfo = new FASTA {
                Header = reads[0].Replace("\n", ""),
                Seq = seq,
                Len = seq.Length };
            return fastaInfo;
        }

        /// <summary>
        /// Read FASTA file, which will contain Multiple FASTA sequence
        /// </summary>
        /// <param name="path">full path of the <b>FASTA/TXT</b> file</param>
        /// <returns>
        /// multi variable return like; <c>.MultiHeader</c>, <c>.MultiSeq</c> and <c>.readLen</c> No. of seq in file
        /// </returns>
        public FASTA multiFASTAread(string path) {
            List<string> header = new List<string>();
            string seqs = "";
            string[] reads = File.ReadAllLines(path);
            foreach (string line in reads) {
                if (line.Contains('>')) {
                    header.Add(line); seqs += " ";
                }
                else { seqs += line; }
            }
            string[] seqarr = seqs.Split(' ');
            var fastainfo = new FASTA {
                readLen = seqs.Split(' ').Length - 1,
                MultiHeader = header.ToArray(),
                MultiSeq = seqarr.Skip(1).ToArray()
            };
            return fastainfo;
        }

        /// <summary>
        /// Download PDB or CIF file format from RCSB.org (make sure, internet should be connected)
        /// </summary>
        /// <param name="pdbId">Alphanumeric PDB id from RCSB</param>
        /// <param name="format">file format you want to download means, format=<b>"pdb" or "cif"</b></param>
        public void PDBFileFormat(string pdbId, string format = "pdb")
        {
            // <a href="//files.rcsb.org/download/4U5X.cif">PDBx/mmCIF Format</a> pdb file download link
            if (format == "cif") {
                Process.Start(new ProcessStartInfo() {
                    FileName = $"https://files.rcsb.org/download/{pdbId}.cif",
                    UseShellExecute = true
                });
            } else {
                Process.Start(new ProcessStartInfo() {
                    FileName = $"https://files.rcsb.org/download/{pdbId}.pdb",
                    UseShellExecute = true
                });
            }
        }

        /// <summary>
        /// Read FASTQ file as sequence and ASCII of pairs
        /// </summary>
        /// <param name="path">full path of the <b>FASTQ</b> file</param>
        /// <returns>Dictionary of strings where key as sequence and value as ASCII code</returns>
        public Dictionary<string, string> ReadFASTQ(string path) {
            List<string> seqs = new(); List<string> ascii = new();
            var pairsDNA_ASCII = new Dictionary<string, string>();
            string[] reads = File.ReadAllLines(path);
            bool cont = true;
            for (int line = 1; line < reads.Length + 1; line++) {
                if (line % 2 != 0) {
                    if (cont == true) {
                        seqs.Add(reads[line]);
                        cont = false;
                    } else { ascii.Add(reads[line]); cont = true; }
                } }
            for (int i = 0; i < seqs.Count; i++) {
                try { pairsDNA_ASCII.Add(seqs[i], ascii[i]); }
                catch (System.ArgumentException) { continue; }   // Argument exception occurs in the case of duplicate key found
            }
            return pairsDNA_ASCII;
        }

        /// <summary>
        /// Applicable for Multiple FASTA sequence file to convert numeric dataset and save as CSV file.
        /// columns of the CSV contain; id of the sequence, GC%, AT%, Length of the sequence number of A,T,G and C
        /// </summary>
        /// <param name="fastaPath"><b>FASTA/TXT</b> file path</param>
        /// <param name="path_csv">output location where you want to save CSV</param>
        public void FASTAToDataset(string fastaPath, string path_csv) {
            BioTools BioInstanse = new();
            var aId = new List<string>(); List<int> nT = new();
            var gcp = new List<float>(); List<int> nG = new();
            var atp = new List<float>(); List<int> nC = new();
            var lens = new List<int>(); List<int> nA = new();

            var fasta = multiFASTAread(fastaPath);
            foreach (string header in fasta.MultiHeader) {
                string[] assId = header.Split(" ");
                aId.Add(assId[0].Replace(">", ""));
            }
            foreach (string seq in fasta.MultiSeq) {
                var content = BioInstanse.ATGC_Content(seq);
                float gc_percent = BioInstanse.GCP(seq); float at_percent = BioInstanse.ATP(seq);
                gcp.Add(gc_percent); lens.Add(seq.Length); atp.Add(at_percent);
                nA.Add(content.Item1); nT.Add(content.Item2); nG.Add(content.Item3); nC.Add(content.Item4);
            }
            using StreamWriter writer = new(path_csv);
            writer.WriteLine(" ,Aid,GC_percent,AT_percent,IndividualLength,no_A,no_T,no_G,no_C");
            for (int i = 0; i < aId.Count; i++)
            {
                writer.WriteLine($"{i},{aId[i]},{gcp[i]},{atp[i]},{lens[i]},{nA[i]},{nT[i]},{nG[i]},{nC[i]}");
                if (i < 5) { Console.WriteLine($"{i},{aId[i]},{gcp[i]},{atp[i]},{lens[i]},{nA[i]},{nT[i]},{nG[i]},{nC[i]}"); }
            }
        }

        /// <summary>
        /// You can save alignment output as <i>sysali</i> or <i>ali</i> or <i>txt</i> file
        /// </summary>
        /// <param name="fileLoc">file name where you want to save alignment</param>
        /// <param name="singleFASTALoc"><b>FASTA</b> file location for single target sequence</param>
        /// <param name="multifastaRefFileLoc"><b>FASTA</b> file location for multiple template sequence</param>
        public void AlignmentFile(string fileLoc, string singleFASTALoc, string multifastaRefFileLoc) {
            BioTools tools = new();
            var target = Readfasta(singleFASTALoc); var template = multiFASTAread(multifastaRefFileLoc);
            var alignment = tools.MultipleAlignment(target.Seq, template.MultiSeq);
            using StreamWriter file = new(fileLoc, false);
            for (int i = 0; i < alignment.similarity.Length; i++) {

                file.WriteLine('>' + template.MultiHeader[i] + " target_len " + target.Seq.Length + " temp_len " + template.MultiSeq[i].Length);
                file.WriteLine("~ Alignment_length: " + alignment.target[i].Length + $" temp count {i}");
                file.WriteLine("~ target_gap " + alignment.target[i].Count(gap => gap == '-') +
                    " temp_gap " + alignment.temp[i].Count(gap => gap == '-') + " identical " + alignment.similarity[i].Count(s => s == '|'));

                file.WriteLine("target    :" + alignment.target[i]);
                file.WriteLine("identical :" + alignment.similarity[i]);
                file.WriteLine($"temp      :" + alignment.temp[i]);
                file.WriteLine("~ encs" + i + " encode_len " + alignment.encodes[i].Length);
                file.WriteLine("encodes   :" + alignment.encodes[i]);
                file.WriteLine(">>>" + '\n');
            }
        }

        /// <summary>
        /// To download GSE (Gene Expression Series) file or dataset
        /// </summary>
        /// <param name="gseAccession">GSE id</param>
        /// <param name="outputLoc"> location of the output will be store</param>
        /// <param name="soft">true if SOFT format needed or otherwise false</param>
        public void GetGEO(string gseAccession, string outputLoc, bool soft = false)
        {                                               // AdditionalFiles={soft, miniml, matrix}
            string dynGSEcode = GSEXXnnn(gseAccession);
            var GEOFile = GetGEOFileName(gseAccession);

            foreach (string filename in GEOFile)
            {
                string url = $"https://ftp.ncbi.nlm.nih.gov/geo/series/{dynGSEcode}/{gseAccession}/suppl/{filename}";
                using (var client = new WebClient())
                {
                    client.DownloadFile(url, outputLoc + '\\' + filename);
                }
            }

            if (soft)
            {
                string gunZipFile = outputLoc + '\\' + gseAccession + "_family.soft.gz";
                using (var client = new WebClient())
                {
                    client.DownloadFile($"https://ftp.ncbi.nlm.nih.gov/geo/series/{dynGSEcode}/{gseAccession}/soft/{gseAccession}_family.soft.gz",
                        gunZipFile);
                }
            }
        }

        private string GSEXXnnn(string gseAccession)
        {
            string dynamicGSEnnn = "";
            for (int i = gseAccession.Length - 4; i > -1; i--) { dynamicGSEnnn += gseAccession[i]; }
            char[] dynToChar = dynamicGSEnnn.ToCharArray();
            Array.Reverse(dynToChar);
            string dynGSEcode = new string(dynToChar) + "nnn";
            return dynGSEcode;
        }

        private List<string> GetGEOFileName(string gseID)
        {
            descriptiveAnalysis duplicateRemoval = new descriptiveAnalysis();
            var doc = new HtmlWeb();
            var con = doc.Load($"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gseID}");
            var node = con.DocumentNode.SelectNodes("html/body/table/tr/td/table");
            List<string> fileName = new List<string>();
            string[] CacheArr = node[5].InnerText.Trim().Split('\n');
            for (int i = 0; i < CacheArr.Length; i++) {
                if (CacheArr[i] == "File type/resource")
                {
                    fileName.Add(CacheArr[i + 2]);
                }
                if (CacheArr[i].Contains(".gz") || CacheArr[i].Contains(".tar") ||
                    CacheArr[i].Contains(".xlsx")) { fileName.Add(CacheArr[i]); }
            }
            List<string> removedDuplicateFiles = duplicateRemoval.RemoveDuplicate(fileName);
            return removedDuplicateFiles;
        }

        /// <summary>
        /// Download GSE (Gene Expression Series) file by reading SOFT file, it will download SOFT file by default
        /// </summary>
        /// <param name="gseAccession">GSE id</param>
        /// <param name="outputLoc">output where you want to store</param>
        public void GetGEOBySOFT(string gseAccession, string outputLoc)
        {
            string dynGSEcode = GSEXXnnn(gseAccession);
            string gunZipFile = outputLoc + '\\' + gseAccession + "_family.soft.gz";
            using (var client = new WebClient())
            {
                client.DownloadFile($"https://ftp.ncbi.nlm.nih.gov/geo/series/{dynGSEcode}/{gseAccession}/soft/{gseAccession}_family.soft.gz",
                    gunZipFile);
            }
            string remoteSoftFileLoc = outputLoc + '\\' + gseAccession + "_family.soft";
            Ungzip(gunZipFile, remoteSoftFileLoc);
            var softReader = ReadMetaDataBySOFT(remoteSoftFileLoc);
            string[] filename = softReader.SeriesSupplementaryFile.Split('|');

            foreach (string file in filename)
            {
                Console.WriteLine(file.Trim());
                using (var client = new WebClient())
                {
                    client.DownloadFile(file.Trim(), outputLoc + '\\' + file.Split('/').Last());
                }
            }
        }

        /// <summary>
        /// Get single Microarray expression
        /// </summary>
        /// <param name="GSMID">GEO sample id</param>
        /// <returns>List Of List, of prob ids and expression values</returns>
        public List<List<object>> IGetExpression(string GSMID)
        {
            var doc = new HtmlWeb();
            var con = doc.Load($"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&acc={GSMID}");
            var node = con.DocumentNode.SelectNodes("/html/body/font/pre/text()[2]");

            List<List<object>> resultExpression = new();
            List<object> head = new(){"ID", GSMID};
            resultExpression.Add(head);

            List<string> affyValues = node[0].InnerText.Replace("\n", "\t").Split('\t').ToList();

            List<object> affyDoubleValues = new();
            foreach (string v in affyValues) {
                if (double.TryParse(v, out _)) {
                    affyDoubleValues.Add(v); // double.Parse(v)
                } else { affyDoubleValues.Add(v); }
            }

            for (int s = 1; s < affyDoubleValues.Count-1; s+=2) 
            {
                var affyEmpty = new List<object>()
                { affyDoubleValues[s], affyDoubleValues[s + 1] }; 
                resultExpression.Add(affyEmpty);
            }
            return resultExpression;
        }

        /// <summary>
        /// Array
        /// </summary>
        /// <param name="SoftFilePath"></param>
        /// <param name="outputPath"></param>
        public void ArrayExpressions(string SoftFilePath, string outputPath)
        {
            // object[ , ] ExpressionValue = new object[ , ];
            Console.WriteLine("1. Work loaded under process..");

            var METADATA = ReadMetaDataBySOFT(SoftFilePath);
            string[] GSMs = METADATA.SeriesSampleID.Split(',').SkipLast(1).ToArray();
            Console.WriteLine($"2. GSM Readed For {METADATA.SeriesTitle.Trim()}..");

            BioSySNet.DataFrame expr = new(IGetExpression(GSMs[0]));
            Console.WriteLine("3. Expression Extracting..");

            Console.Write(GSMs[0]+" => ");
            for (int i = 1; i < GSMs.Length; i++)
            {
                List<List<object>> Iexpr = IGetExpression(GSMs[i]);
                expr.Merge(Iexpr, ByColName: "ID");
                Console.Write(GSMs[i]+" => ");
            }

            Console.WriteLine($"\n4. Created: {outputPath+METADATA.SeriesAccession.Trim()}_Expr.csv");
            using(var csv = File.CreateText(outputPath+METADATA.SeriesAccession.Trim()+"_Expr.csv"))
            {
                foreach (object head in expr.HEADER) { csv.Write(head+","); }
                csv.WriteLine();
                Console.WriteLine("Shape: "+expr.Data.Count+" x "+expr.Data[1].Count);
                for (int i = 0; i < expr.Data.Count; i++)
                {
                    for (int j = 0; j < expr.Data[i].Count; j++){ csv.Write(expr.Data[i][j]+","); }
                    csv.WriteLine();
                }
                csv.Close();
            }
            Console.WriteLine("-- Finished --");
        }

        public void GetGEOPlatform(string platformID)
        {
            
        }

        public struct SoftFileData
        {
            public string? SeriesTitle, SeriesAccession, SeriesSummary;
            public string? SeriesSampleID, SeriesSupplementaryFile;

            public string? PlatformID, PlatformTitle, PlatformTechnology, PlatformOrganism;

            public Dictionary<string, string> SAMPLES;
        }

        public SoftFileData ReadMetaDataBySOFT(string filename)
        {
            string[] lines = File.ReadAllLines(filename);
            SoftFileData softFileData = new SoftFileData();
            for (int i = 0; i < lines.Length; i++)
            {
                if (lines[i].Contains("=")) { 
                    string ValueLine = lines[i].Split('=')[1].Replace('\n', ' ');
                    if (lines[i].Contains("!Series_title")) {
                        softFileData.SeriesTitle = ValueLine; }
                    else if (lines[i].Contains("!Series_geo_accession")) {
                        softFileData.SeriesAccession = ValueLine; }
                    else if (lines[i].Contains("!Series_summary")) {
                        softFileData.SeriesSummary = ValueLine; }
                    else if (lines[i].Contains("!Series_sample_id")) {
                        softFileData.SeriesSampleID += ValueLine.Trim()+','; }
                    else if (lines[i].Contains("!Series_supplementary_file")) {
                        softFileData.SeriesSupplementaryFile += ValueLine + '|'; }
                    else if (lines[i].Contains("!Platform_geo_accession")) {
                        softFileData.PlatformID = ValueLine; }
                    else if (lines[i].Contains("!Platform_title")) {
                        softFileData.PlatformTitle = ValueLine; }
                    else if (lines[i].Contains("!Platform_technology")) {
                        softFileData.PlatformTechnology = ValueLine; }
                    else if (lines[i].Contains("!Platform_organism")) {
                        softFileData.PlatformOrganism = ValueLine; }
                }
                else { continue; }
            }
            return softFileData;
        }

        private void Ungzip(string gzipLocation, string Newloc)
        {
            FileInfo file = new FileInfo(gzipLocation);
            using (FileStream originalStream = file.OpenRead()) {
                using (FileStream outputStream = File.Create(Newloc)) {
                    using (GZipStream gZipStream = new GZipStream(originalStream,
                                                                 CompressionMode.Decompress))
                    {
                        gZipStream.CopyTo(outputStream);
                    }
                }
            }
        } 
        // matrix file, SAM file
    }
}