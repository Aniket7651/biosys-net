
//////////// ////// //////// STATEMENTS UNDER TESTING BY CONSOLE APPLICATION //////// ////// ///////////////
//////////// ////// ////////        PAKAGING PERFORM BY CLASS LIBRARY        //////// ////// ///////////////
# pragma warning disable

string GlobalOutPut = "A:\\PROJECTS\\C#_Programs\\BioSySNet\\BioSySNet\\bin\\OutPut\\";

//BioSySNet.drugScraping scrap = new BioSySNet.drugScraping();
BioSySNet.BioFormats fileformat = new BioSySNet.BioFormats();
//BioSySNet.BioTools seqtools = new BioSySNet.BioTools();
// PyConfiguration.Config config = new PyConfiguration.Config();
BioSySNet.CSVArray csvArr = new();
// BioSySNet.BioFormats formats = new();
// var soft = formats.ReadSoftMetaData("C:\\Users\\lenovo\\Downloads\\GSE17193_family.soft");
// foreach (string v in soft.SeriesSampleID.Split(',').SkipLast(1)){ Console.WriteLine(v); }

BioSySNet.ImportData import = new();
// List<object> newl = new List<object>(){2,6,8,3,2};
// fileformat.ArrayExpressions("C:\\Users\\lenovo\\Downloads\\GSE17193_family.soft", GlobalOutPut);
// var data = import.FromCommaDelim("A:\\PROJECTS\\C#_Programs\\BioSySNet\\BioSySNet\\sample.csv");
// var data2 = import.FromCommaDelim("A:\\PROJECTS\\C#_Programs\\BioSySNet\\BioSySNet\\s2.csv");
// BioSySNet.DataFrame affyFrame = new(data, defaultTypes: true);
// Console.WriteLine(data[1][1].GetType());
// BioSySNet.DataFrame dataFrame = new();
// dataFrame.Merge(data2, "h1");
// // dataFrame.Stats(saveLogs: "A:\\PROJECTS\\C#_Programs\\BioSySNet\\BioSySNet\\");
// foreach (var row in dataFrame.Stats())
// {
//      foreach (object item in row)
//      {
//           Console.Write(item + "\t\t");
//      } Console.WriteLine();
// }
// affyFrame.SortedIndxInCol(newl);
// foreach (object item in data[0])
//      {
//           Console.Write(item + "\t\t");
//      }
//BioSySNet.BioTools.QualityControlling QC = new();
//BioSySNet.descriptiveAnalysis arraysFunc = new BioSySNet.descriptiveAnalysis();
//BioSySNet.GEOCSVTable csv = new();
//BioSySNet.Statistics sts = new BioSySNet.Statistics();
// BioSySNet.GEOAnalysis geo = new();
// config.RequiredModules();
// List<List<object>> affyData = geo.IGetExpression(GSMID: "GSM3462691"); // GSM3462714

// affyFrame.DuplicatesIndx(newl);
// affyFrame.SaveFrameAsCSV("A:\\PROJECTS\\C#_Programs\\BioSySNet\\BioSySNet\\expr.csv");
// BioSySNet.DataFrame dataFrame = new();
// Console.WriteLine(dataFrame.Data);

// affyFrame.SortFrame(1);
// object[ , ] arr = affyFrame.ListToArray();
// csvArr.visualizeCSV(arr);
// foreach (var row in affyFrame.SortFrame(1))
// {
//      foreach (object item in row)
//      {
//             Console.Write(item + "\t");
//      } Console.WriteLine();
// }
// affyFrame.Info();
// var soft = fileformat.ReadMetaDataBySOFT("B:\\Major_project\\dataFolder\\GSE129046_family\\GSE129046_family.soft");
//Console.WriteLine(soft.SeriesTitle);
//Console.WriteLine(soft.SeriesSampleID);
//Console.WriteLine(soft.PlatformOrganism);
// Console.WriteLine(seqtools.OrganicMolecular_Weight("C3H3"));
//string softfile = "B:\\Major_project\\dataFolder\\GSE129046_family\\GSE129046_family.soft";
// var soft = fileformat.ReadSoftFile(softfile);
// Console.WriteLine(soft.SeriesSupplementaryFile);
//double[,] data = new double[,] {{5,4,3}, 
//                                {2,1,4},
//                                {3,4,6},
//                                {4,2,8}};

// string InputFile = "B:\\Major_project\\dataFolder\\GSE154844_family\\GSE154844_Raw_counts.csv";
// string OutputLoc = "B:\\Major_project\\dataFolder\\GSE154844_family\\Analysis";
// BioSyS_Tool.ExpressionNormalization expression = new();
// double[,] trns = expression.FPKMNormalization(ct, new double[] {20, 10, 10, 30}, 10);
// var sample = csv.IndexSeparatedCounts("B:\\Major_project\\dataFolder\\GSE154844_family\\Analysis\\GEO_QNNormal.csv", new int[] {1,2,3,4,5,6});
// var table = csv.IndexSeparatedCounts("A:\\GSE129046_gene.csv", new int[] {6,7,8,9,10,11,12,13,14,15,16,17});
// geo.GEOParameters(InputFile, new int[] {1,2,3,4,5,6}, OutputLoc, "QN");
// double[,] corr = sts.CorrelationMatrix(arraysFunc.Transpose(sample.CountData), OutputLoc, sample.Accession);
// var hell = geo.sampleSplit(6);
//Console.WriteLine(5 - (Math.Abs(-5))); // to convert negative number to positive
// geo.PadjValue(new double[] {0.01,0.05,0.10,0.01,0.07,0.03});
// foreach(double val in table.GeneLength) { Console.WriteLine(val); }
// foreach (string col in ctt.Accession) { Console.Write(col + "\t"); }
//Console.WriteLine("\n\n");

//Console.WriteLine('\n');
//for (int row = 0; row < sample.Sample2.GetLength(0); row++)
//{
//    for (int col = 0; col < sample.Sample2.GetLength(1); col++)
//    {
//        Console.Write(sample.Sample2[row, col] + "\t");
//    }
//    Console.WriteLine();
//}

// List<double> list = new List<double> { 4, 1, 4, 2 };
// foreach(double i in arraysFunc.RemoveDuplicate(list)){ Console.Write(i+", "); }
//int[]idx = arraysFunc.ArgSort(list);
//for (int i = 0; i < idx.Length; i++) { Console.WriteLine(idx[i]); }
//fileformat.GetGEOBySOFT("GSE154844", "A:\\outp");
// string gseAccession = "GSE129046";
//double[] ar = new double[]{2.3, 3.4, 4.32, 5.32};
//double[] ar2 = new double[] { 2, 3 };
//double[] s1 = new double[] { 3841.49, 4203.07, 5182.89 };
//double[] s2 = new double[] { 3396.28, 3488.78, 4660.11 };
// Console.WriteLine(geo.Pvalue(1.49, 6));
// Console.WriteLine(geo.TwoTailedPvalue(1.49, 6));

// foreach(double i in geo.ValueBetween(ar, 3.0)) { Console.WriteLine(i); }
//double fc = geo.TStatistics(new double[] { 5,4}, new double[] { 3, 3 }).TStatistics;


//int n = 2;    //size of sample
//double hypotheziedMean = sts.mean(new double[] { 3, 3 });
//double sampleMean = sts.mean(new double[] { 5, 4 });
//double sampleSD = sts.stdDev(new double[] {5,4});    //sample standard deviation (n-1)

//double stdErr = sampleSD / Math.Sqrt(n);    //standard error of the mean
//double t = (sampleMean - hypotheziedMean) / stdErr;   //convert to a standard mean of 0 and SD of 1

//Console.WriteLine();
//Console.WriteLine(fc.log2FoldChange + ": " + fc.regulation + ": " + fc.foldChange);
//string seq = "HJAHSHCHHCSYSYSSHDBDYATSTDUTUA";
//var sys = seqtools.ProteinContent(seq);
//foreach (var a in sys)
//{
//    Console.WriteLine(a);
//}
//string MultyFASTAfile = @"A:\\PROJECTS\\C#_Programs\\Bioinfo\\BioSyS_Tool\\bin\\Debug\\net6.0\\docs\\samples\\sequence.fasta";
//string SingleFASTAfile = @"A:\\PROJECTS\\C#_Programs\\Bioinfo\\BioSyS_Tool\\bin\\Debug\\net6.0\\docs\\samples\\samples\\dengu3_strain-D00-0107.fasta";
// string pdbfile = "A:/PROJECTS/4u5x.pdb";
// fileformat.PDBReader(pdbfile);
// string fileM = "A:/PROJECTS/test.txt";
// var tarfasta = funcs.multiFASTAread(fileM).MultiSeq;
// var tempfasta = funcs.multiFASTAread(file);

// var fasta = funcs.Readfasta(file);
// float ecd = arraysFunc.EuclideanDistance("ATGTGTG", "CATGTG");
// Console.WriteLine(" " + ecd);
// var lis = seqtools.TripletFormer("ATGTGTG", Alter: 1);
// foreach(string s in lis){
//     Console.WriteLine(s);
// }
// float gc_pecent = seqtools.GCP("GACGTGCATGCATTATGCGAATTAGCTATAAGCT");
// Console.WriteLine(gc_pecent);
// var cds = seqtools.CDS(fasta.Seq);
// Console.WriteLine(cds.stopPosition);
// var dict = funcs.readFASTQ(@"A:\\PROJECTS\\out_1000_ERR101899.1.fastq");

// int[,] arr = new int[,]{{0,1,2,3}, {2,3,4,5}, {4,5,3,5}, {4,5,6,7}};
// float[] me = seqtools.MeanScore(arr);
// seqtools.SaveScores(arr, me, "outscore.csv");
// int[] shape = new int[] { 2, 4 };
// int[,] Tarr = arraysFunc.Transpose(arr, shape);
// for (int i = 0; i < 4; i++)
// {
//     for (int j = 0; j < 2; j++)
//     {
//         Console.WriteLine($"row {i}, col {j}, val: {Tarr[i, j]}");     ACGAA  AGCGA
//     }
// }"GTAAGAA"   "GTAAACGA"                                                A-CACACTA
//                                                                        | ||||| |
//                                                                        AGCACAC-A   RIGHT CASE
//string align_file = "align.sysalign";                                //  CCGT---AC-TA
// funcs.AlignmentFile(align_file, SingleFASTAfile, MultyFASTAfile);      // CAGACCTA      WRONG CASE
// var align = seqtools.SmithWatermanAlignment("ACGAA", "AGCGA");
// funcs.AlignmentFile(align_file, align.target, align.pairs, align.template);
// for(int i = 0; i < align.GetLength(0); i++){
//     for(int j = 0; j < align.GetLength(1); j++){
//         Console.Write(align[i,j]+"  ");
//     }Console.WriteLine();
// }  
// "ATGGTTGCCACGTAGGGCGGTCGAAAGTCGCCCCCCTCGTTCAA"
// var fastq = fileformat.ReadFASTQ("A:/PROJECTS/out_1000_ERR101899.1.fastq");
// int[ , ] qaulityProbapility = QC.Quality(fastq);

// for(int i = 0; i < qaulityProbapility.GetLength(0); i++){
//     for(int j = 0; j < qaulityProbapility.GetLength(1); j++){
//          Console.Write(qaulityProbapility[i,j]+"   ");
//     } Console.WriteLine();
// }
// scrap.drugId("paracetamol");
// var maxI = arraysFunc.GetMaxInMatrix(align);
// Console.WriteLine("r " + maxI.Item1+" c "+maxI.Item2+" value "+maxI.Item3);
// Console.WriteLine(align.target);
// Console.WriteLine(align.pairs);
// Console.WriteLine(align.template);
// Console.WriteLine(align.encode);
// int[,] arr = new int[,] { { 2, 3, 4, 5 }, { 2, 3, 4, 5 },{3,4,5,6} };

// funcs.PDBFileFormat("4U5X");
// scrap.UniProtSeq("P05067");

// funcs.FastaDataset(fastaPath: "A:/PROJECTS/test.txt", path_csv:"A:/PROJECTS/sample.csv");
//float Mw = seqtools.MeltingPoint("ACGTACGTACG");
//Console.WriteLine(" " + Mw);
// scrap.fastaFile("OQ622003.1");