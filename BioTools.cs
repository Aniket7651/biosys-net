/*
~BioTools.cs

@ Author: ANIKET YADAV (aniketyadav8687@gmail.com)
*/

using IronPython.Hosting;
using Microsoft.Scripting.Hosting;
using IronPython.Runtime;
using PyConfiguration;

namespace BioSySNet
{
    public class BioTools
    {
        //public static dynamic Make_con_pyScript(string pyscript, string funcName){
        //    string path = AppDomain.CurrentDomain.BaseDirectory + "pyscript\\"+pyscript;
        //    ScriptEngine engine = Python.CreateEngine();
        //    ScriptScope scope = engine.CreateScope();
        //    engine.ExecuteFile(path, scope);
        //    dynamic funcSelector = scope.GetVariable(funcName);
        //    return funcSelector;
        //}

        public Tuple<int, int, int, int> ATGC_Content(string ntSequence){
            int a = ntSequence.Count(a => a=='A');
            int t = ntSequence.Count(t => t=='T');
            int g = ntSequence.Count(g => g=='G');
            int c = ntSequence.Count(c => c=='C');
            return Tuple.Create(a, t, g, c);
        }

        public Dictionary<char, int> ProteinContent(string ptSequnce){
            List<char> seqLis = new();
            foreach(char seq in ptSequnce) { seqLis.Add(seq); }

            Dictionary<char, int> Content = new Dictionary<char, int>();
            descriptiveAnalysis analysis = new descriptiveAnalysis();
            char[] uniqueAA = analysis.UniqueValue(seqLis.ToArray());
            
            foreach (char aa in uniqueAA) {
                int count = 0;
                for (int i = 0; i < ptSequnce.Length; i++) {
                    if (ptSequnce[i] == aa) { 
                        count ++;
                    }
                } try { Content.Add(aa, count); } catch (ArgumentException) { continue; }
            } return Content;
        }
        
        public float GCP(string DNAsequence){
            string cleanCase = DNAsequence.ToUpper();
            float gc = cleanCase.Count(g => g == 'G')+cleanCase.Count(c => c == 'C');
            float gcp = (gc/cleanCase.Length)*100;
            return gcp;
        }

        public float ATP(string DNAsequence){
            string cleanCase = DNAsequence.ToUpper();
            float at = cleanCase.Count(a => a == 'A')+cleanCase.Count(t => t == 'T');
            float atp = (at/cleanCase.Length)*100;
            return atp;
        }

        public string Complementary(string ntSequence, bool reverse=false){
            string complement = "";
            var compKey = new Dictionary<char, char>(){{'A', 'T'},{'T','A'},{'G','C'},{'C','G'}};
            foreach(char nt in ntSequence){
                complement += compKey[nt];
            }
            if (reverse == true){
                string rev = "";
                for(int i=complement.Length-1; i>-1; i--){
                    rev += complement[i];
                } complement = rev;
            }
            return complement;
        }

        public List<string> K_merSequence(string ntSequence, int K=6){
            List<string> kmer = new();
            for(int i=0; i<ntSequence.Length; i++){
                try{ kmer.Add(ntSequence.Substring(i, K)); } catch(ArgumentOutOfRangeException){ continue; }
            }
            return kmer;
        }

        public struct Cds{
            public string? CDSSeq;
            public int startPosition, stopPosition;
        }
        public Cds CDS(string ntSequence, string stopSelection="TGA"){
            string substringCDS = "";
            int startCodon = ntSequence.IndexOf("ATG");
            int stopCodon = ntSequence.LastIndexOf(stopSelection);

            if (startCodon<stopCodon){
                substringCDS = ntSequence.Substring(startCodon, stopCodon);
            }
            var cdsInfo = new Cds{CDSSeq = substringCDS,
                                startPosition = startCodon, stopPosition = stopCodon};
            return cdsInfo;
        }

        public struct Alignment{
            public string target, template, pairs; public string encode; 
        }

        public double[ , ] SmithMatrixAlign(string seq, string tempseq, int match=1, int mismatch=-1, int gap=-2){
            descriptiveAnalysis array = new();
            int Match(char s1, char s2){
                int x = s1==s2 ? match : mismatch;
                return x;
            }
            int m = seq.Length+1; int n = tempseq.Length+1;
            double[ , ] ScoringMatrix = array.Matrix(m, n);
            for(int i = 1; i < m; i++){
                for(int j = 1; j < n; j++){
                    double Digon = ScoringMatrix[i-1, j-1]+Match(seq[i-1], tempseq[j-1]);
                    double Vertic = ScoringMatrix[i-1, j]+gap;
                    double Horizon = ScoringMatrix[i, j-1]+gap;
                    double[] argsHVD0 = new double[4]{Digon, Vertic, Horizon, 0};
                    ScoringMatrix[i, j] = array.Argmax(argsHVD0);
                }
            }
            return ScoringMatrix;
        }

        public static Tuple<double[], int[]> Traceback(double[ , ] scoreMatrix){
            descriptiveAnalysis arraysFunc = new();
            List<int> MaxIndexList = new();
            var maxIndex = arraysFunc.GetMaxInMatrix(scoreMatrix);
            int i = scoreMatrix.GetLength(0)-1; int j = scoreMatrix.GetLength(1)-1;
            List<double> maxList = new();
            for(;i>-1;){
                try{
                    double H = scoreMatrix[i,j-1]; double V = scoreMatrix[i-1,j]; double D = scoreMatrix[i-1,j-1];
                    //              0=horizon(H), 1=verical(V), 2=digonal(D)
                    double[] HVDMAX = new double[3]{H, V, D}; 
                    double max = arraysFunc.Argmax(HVDMAX); int MaxIndex = Array.IndexOf(HVDMAX, max);
                    maxList.Add(max);
                    MaxIndexList.Add(MaxIndex);
                    // Console.WriteLine(max+"  "+MaxIndex);

                    if(MaxIndex == 2){ i--;j--; }
                    if(MaxIndex ==1){ i--; }
                    if(MaxIndex == 0){ j--; }
                }
                catch(IndexOutOfRangeException){ break; }
            }
            return Tuple.Create(maxList.SkipLast(1).ToArray(), MaxIndexList.SkipLast(1).ToArray());
        }
        static string EncodeChanger(int[] encodeArray){
            string encode = "";
            for(int i = encodeArray.Length-1; i>-1; i--){
                encode += encodeArray[i];
            } return encode;
        }
        
        public Alignment SmithWatermanAlignment(string target, string temp, int match=1, int mismatch=-1, int gap=-2, bool matrix=false){
            
            double[ , ] mat = SmithMatrixAlign(target, temp, match, mismatch, gap);
            int[] indexScore = Traceback(mat).Item2;
            if(matrix){
                for (int i = 0; i < mat.GetLength(0); i++){
                    for (int j = 0; j < mat.GetLength(1); j++){
                        Console.Write(mat[i,j]+"   ");
                    }Console.WriteLine('\n');
                }
            }
            static string reverse(string str){
                string rev = "";
                for(int a=str.Length-1;a>-1;a--){
                    rev += str[a];
                } return rev;
            }
            static string pairsAlign(string s1, string s2){
                string pairs = "";
                for(int i=0;i<s1.Length; i++){
                    try{pairs += (s1[i]==s2[i]&& s1[i] !='-') ? "|": " ";} catch{break;} 
                }
                return pairs;
            }
            string targetreverse = reverse(target); string tempreverse = reverse(temp);
            for (int i = 1; i < indexScore.Length + 1; i++)
            {
                if(indexScore[i-1]==0){
                    targetreverse = targetreverse.Insert(i, "-");
                }else if(indexScore[i-1]==1){
                    tempreverse = tempreverse.Insert(i, "-");
                }
            }
            string targetAligned = reverse(targetreverse); string tempAligned = reverse(tempreverse);
            string pairs = pairsAlign(targetAligned, tempAligned);
            var alignment = new Alignment{ 
                target = reverse(targetreverse), pairs = pairs, template = reverse(tempreverse), encode = EncodeChanger(indexScore)
            };
            return alignment;
        }

        public List<string> SequenceDistributor(string ntSequence, int Alter=3, int len=3){
            List<string> tripCode = new();
            for(int i=0;i<ntSequence.Length;i+=Alter){
                try{
                    string triplet = ntSequence.Substring(i, len);
                    tripCode.Add(triplet);
                }
                catch(ArgumentOutOfRangeException){ continue; }  
            }
            return tripCode;
        }

        public string Translation(string seq, int startIndex=0, int stopIndex=0){
            var code = new Dictionary<string, string>(){
                {"TTT","F"},{"TTC","F"},{"TTA","L"},{"TTG","L"},{"TCT","S"},{"TCC","S"},{"TCA","S"},{"TCG","S"},
                {"TAT","Y"},{"TAC","Y"},{"TAA","_"},{"TAG","_"},{"TGT","C"},{"TGC","C"},{"TGA","_"},{"TGG","W"},

                {"CTT","L"},{"CTC","L"},{"CTA","L"},{"CTG","L"},{"CCT","P"},{"CCC","P"},{"CCA","P"},{"CCG","P"},
                {"CAT","H"},{"CAC","H"},{"CAA","Q"},{"CAG","Q"},{"CGT","R"},{"CGC","R"},{"CGA","R"},{"CGG","R"},

                {"ATT","I"},{"ATC","I"},{"ATA","I"},{"ATG","M"},{"ACT","T"},{"ACC","T"},{"ACA","T"},{"ACG","T"},
                {"AAT","N"},{"AAC","N"},{"AAA","K"},{"AAG","K"},{"AGT","S"},{"AGC","S"},{"AGA","R"},{"AGG","R"},

                {"GTT","V"},{"GTC","V"},{"GTA","V"},{"GTG","V"},{"GCT","A"},{"GCC","A"},{"GCA","A"},{"GCG","A"},
                {"GAT","D"},{"GAC","D"},{"GAA","E"},{"GAG","E"},{"GGT","G"},{"GGC","G"},{"GGA","G"},{"GGG","G"}
            }; // int len = stopIndex == 0 ? seq.Length: stopIndex;
            string translate = "";
            // string subs = "";    // expected subs- " CGT AGC TGC TAT CGT ACG ATC GAT TTA TAG CAT AGA A"
            BioTools instance = new();
            string[] seqCode = instance.SequenceDistributor(ntSequence: seq).ToArray();
        
            foreach(string C in seqCode){
                if (code.ContainsKey(C)){ translate += code[C]; }
            }
            return translate;
        }

        public struct MultiAlign{
            public string[] target, similarity, temp; public string[] encodes;
        }

        public MultiAlign MultipleAlignment(string target, string[] reference, int match=1, int Mis=-1, int gap=-2){
            string[] targetalign = new string[reference.Length]; string[] simalign = new string[reference.Length];
            string[] tempalign = new string[reference.Length]; string[] encodealign = new string[reference.Length];

            for(int i = 0; i < reference.Length; i++){
                var alignment = SmithWatermanAlignment(target, reference[i], match, Mis, gap);
                targetalign[i] = alignment.target;
                simalign[i] = alignment.pairs;
                tempalign[i] = alignment.template;
                encodealign[i] = alignment.encode;

            }var Multialign = new MultiAlign{
                target = targetalign, similarity = simalign, temp = tempalign, encodes = encodealign
            }; return Multialign;
        }

        public class QualityControlling{
            descriptiveAnalysis transformFunction = new();
            private static int[] GetPhred_33(string asciiChar){
                int[] asciiValue = new int[asciiChar.Length];
                for(int c=0; c<asciiChar.Length; c++){
                    asciiValue[c] = (int)(asciiChar[c])-33;
                } return asciiValue;
            }

            private static int[] GetPhred_64(string asciiChar){
                int[] asciiValue = new int[asciiChar.Length];
                for(int c=0; c<asciiChar.Length; c++){
                    asciiValue[c] = (int)(asciiChar[c])-64;
                } return asciiValue;
            }

            public int[ , ] Quality(Dictionary<string,string> dic){
                static bool IsBased64(Dictionary<string,string> dic){
                    foreach(string k in dic.Keys){
                        string ASCII = dic[k];
                        foreach(char c in ASCII){
                            if(char.IsLower(c)){
                                return true; }
                        }
                    } return false;
                }
                List<int[]> ScoringList = new List<int[]>();
                foreach(string k in dic.Keys){
                    if(IsBased64(dic)){
                        ScoringList.Add(GetPhred_64(dic[k]));
                    }
                    else{ ScoringList.Add(GetPhred_33(dic[k])); }
                } int[][] score = ScoringList.ToArray(); int[,] To2D = transformFunction.JaggedTo2DArray(score);
                int[] shape = new int[2]{To2D.GetLength(0), To2D.GetLength(1)};
                int[ , ] scoreMat = transformFunction.Transpose(To2D);
                
                return scoreMat;
            }

            public float[] MeanScore(int[ , ] scoreMat){
                List<float> mean = new();
                for(int i = 0;i<scoreMat.GetLength(0); i++){
                    float sum = 0;
                    for(int j = 0;j<scoreMat.GetLength(1); j++){
                        sum += scoreMat[i,j];
                    } mean.Add(sum/scoreMat.GetLength(1));
                } return mean.ToArray();
            }

            public void SaveScores(int[ , ] scoreMatrix, float[] meanScore, string outputSaveLoc){
                using StreamWriter writeCSV = new StreamWriter(outputSaveLoc, false);
                for(int i = 0; i<scoreMatrix.GetLength(0); i++){
                    for(int j = 0; j<scoreMatrix.GetLength(1); j++){
                        writeCSV.Write(scoreMatrix[i,j]+",");
                    } writeCSV.WriteLine(meanScore[i]);
                } Console.WriteLine("Created score as CSV path "+outputSaveLoc);
            }
            
        }

        Config config = new Config();

        public float SsNtMolecularWeight(string nt_seq){
            
            dynamic ntWeight = config.Make_con_pyScript(@"script1.py", "ssDNA_MoleculerWeight");
            var res = ntWeight(nt_seq);
            return (float)res;
        }

        public float MeltingPoint(string nt_seq){
            dynamic ntWeight = config.Make_con_pyScript(@"script1.py", "Melting_Temperature");
            var res = ntWeight(nt_seq);
            return (float)res;
        }

        public float OrganicMolecular_Weight(string emprical_formula){
            dynamic ntWeight = config.Make_con_pyScript(@"script1.py", "Chem_MolecularWeight");
            var res = ntWeight(emprical_formula);
            return (float)res;
        }

        public float Peptide_Mweight(string prot_seq){
            dynamic peptideWeight = config.Make_con_pyScript(@"script1.py", "peptide_Molecular_weight");
            var res = peptideWeight(prot_seq);
            return (float)res;
        }

        public string DNA_to_RNA(string ntSequence){
            return ntSequence.Replace('T', 'U');
        }

        public string LabelEncoding(string nt_seq){
            string labels = "";
            foreach(char c in nt_seq){
                if(c == 'A'){
                    labels += "1"; }
                if(c == 'T'){
                    labels += "2"; }
                if(c == 'G'){
                    labels += "3"; }
                if(c == 'C'){
                    labels += "4"; }
            }
            return labels;
        }
    }
}