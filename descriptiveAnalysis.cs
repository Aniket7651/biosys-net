/*
~descriptiveAnalysis.cs
class file contains all the logic behind arrays and methods like max value in array and min 
values et. like, all the operation perform by this script, is for arrays and csv
hence, this script is use for the suppoter progam script for other logics also.

@ Author: ANIKET YADAV (aniketyadav8687@gmail.com)
*/

using System;
using System.Collections.Generic;
using System.Collections.Immutable;
using System.Linq;
using MathNet.Numerics.Distributions;
using System.Net;
using System.Net.NetworkInformation;
using Microsoft.Scripting.ComInterop;

namespace BioSySNet
{
    public class descriptiveAnalysis
    {
        // minimum value in doubled array if index=true; return index number in double, where most smallest 
        // value exist.
        public double Argmin(double[] col, bool index=false){
            int arrLen = col.Length;
            double min = col[0];
            for (int i=0; i<arrLen; i++){
                if (min>col[i]){ if (index==false){min = col[i]; } else{ min=i; }}
            }
            return min;
        }

        // max value in doubled array if index=true; return index number in double, where most largest 
        // value exist.
        public double Argmax(double[] col, bool index=false){
            int arrLen = col.Length;
            double max = col[0];
            for (int i=0; i<arrLen; i++){
                if (max<col[i]){
                    if (index == false){ max = col[i]; }
                    else{ max = i; }
                }
            }
            return max;
        }

        // applied sorting by index in the place of values
        public int[] ArgSort(List<double> list, string order="asnd") // order: des
        {
            double[] arr = new double[list.Count];
            list.CopyTo(arr);
            if (order == "des") { list.OrderByDescending(x => x); }
            else { list.Sort(); }

            // list.Sort((a, b) => b.CompareTo(a));
            int[] indx = new int[list.Count];

            for (int i = 0; i < arr.Length; i++)
            {
                indx[i] = list.IndexOf(arr[i]);
            }
            return indx;
        }

        // make row to column and column to row by Flip matrix in 2D array; return floating 2D array
        public float[ , ] Flip(float[ , ] matrix){
            int row = matrix.GetLength(0); int col = matrix.GetLength(1);
            float[ , ] emptyFliped = new float[row, col];
            for(int i = row-1; i>-1; i--){
                for(int j = col-1; j>-1; j--){
                    emptyFliped[i, j] = matrix[i, j];
                    Console.Write(matrix[i, j]+"   ");
                }Console.WriteLine('\n');  
            } return emptyFliped;
        }

        public float EuclideanDistance(string seq1, string seq2){
            BioTools tool = new();
            string[] s1 = tool.SequenceDistributor(seq1, 1).ToArray(); string[] s2 = tool.SequenceDistributor(seq2, 1).ToArray();
            string[] seqUnion = s1.Union(s2).ToArray();
            int[] s1Count = new int[seqUnion.Length]; int[] s2Count = new int[seqUnion.Length];
            for(int i=0; i<seqUnion.Length; i++){
                int count_s1 = 0; int count_s2 = 0;
                for(int j=0; j<s1.Length; j++){
                    if(seqUnion[i] == s1[j]){ count_s1++; }
                } s1Count[i] = count_s1;
                for(int j=0; j<s2.Length; j++){
                    if(seqUnion[i] == s2[j]){ count_s2++; }
                } s2Count[i] = count_s2;
            } double sum = 0;
            for(int i=0; i<s1Count.Length; i++){
                int diff = s1Count[i] - s2Count[i];
                sum += diff*diff;
            } return (float)Math.Sqrt(sum);
        }

        public Type[] UniqueValue<Type>(Type[] sting_array){
            int upos = 0;  // Unique string arry's position
            // initialize an empty array
            var UniquestringArr = new Type[sting_array.Length];
            for (int i=0; i<sting_array.Length; i++){
                if (!UniquestringArr.Contains(sting_array[i])){
                    UniquestringArr[upos] = sting_array[i];
                    upos++; }
            }
            return UniquestringArr;
        }

        public float[] OneDarray(int n, float matrixOf=0){
            List<float> zeroList = new List<float>();
            for(int i=0; i<=n; i++){ zeroList.Add(matrixOf); }
            return zeroList.ToArray();
        }

        public Tuple<int, int, double> GetMaxInMatrix(double[ , ] mat){
            double max = 0; int rowIndex = 0; int colIndex = 0;
            for(int i=0; i<mat.GetLength(0); i++){
                for(int j=0; j<mat.GetLength(1); j++){
                    if(mat[i,j] > max){ max = mat[i,j]; rowIndex = i; colIndex = j; }
                }
            }
            return Tuple.Create(rowIndex, colIndex, max);
        }

        public double[,] Matrix(int row, int col, float matrixOf=0){
            double[,] zeromat = new double[row, col];
            for(int i=0; i<row; i++){
                for(int j=0; j<col; j++){
                    zeromat[i, j] = matrixOf;
                }
            }
            return zeromat;
        }

        public Type[ , ] Transpose<Type>(Type[ , ] matrix){
            Type[ , ] Transp = new Type[matrix.GetLength(1), matrix.GetLength(0)];
            for(int i=0; i<matrix.GetLength(1); i++){
                for(int j=0; j<matrix.GetLength(0); j++){
                    Transp[i,j] = matrix[j,i];
                }
            }
            return Transp;
        }

        public Type[ , ] JaggedTo2DArray<Type>(Type[][] matrix){
            int row = matrix.Length;
            int col = matrix[0].Length;
            Type[,] To2DA = new Type[row, col];
            for(int i=0; i<row; i++){
                for(int j=0; j<col; j++){
                    To2DA[i,j] = matrix[i][j];
                }
            } return To2DA;
        }

        public List<Type> RemoveDuplicate<Type>(List<Type> list)
        {
            List<Type> removed = new List<Type>();
            for (int i=0; i<list.Count; i++)
            {
                if (!removed.Contains(list[i])) { removed.Add(list[i]); }
            }
            return removed;
        }
    }


    public class Statistics
    {
        public double mean(double[] x)
        {
            double sum = 0;
            for(int i=0; i < x.Length; i++)
            {
                sum += x[i];
            }
            return sum/x.Length;
        }

        public double median(double[] x)
        {
            double median = (x.Length % 2 == 0)? (x[x.Length / 2] + x[(x.Length / 2) + 1]) / 2
                : x[(x.Length + 1) / 2];
            return median;
        }

        public double stdDev(double[] values)
        {
            int n = values.Length;
            double Xmean = mean(values);
            double sum_mean = 0.0;
            for(int i=0; i<n; i++)
            {
                sum_mean += Math.Pow(values[i] - Xmean, 2.0);
                
            }
            return Math.Sqrt(sum_mean/(n-1));
        }

        public double[] Zscore(double[] values)
        {
            double[] z = new double[values.Length];
            double m = mean(values); double sd = stdDev(values);
            for (int i =0; i<values.Length; i++)
            {
                z[i] = (values[i] - m) / sd;
            }
            return z;
        }

        private double Quartile(double[] x, int quartile)
        {
            Array.Sort(x);
            int n = x.Length;
            int q_i = (int)(n + quartile) / 4;
            double quar = (n % 2 == 0) ? x[q_i - 1] + x[q_i]/2: x[q_i - 1];
            return quar;
        }

        /// <summary>
        /// Identify Potential Outliers in data array by interquartile range method
        /// </summary>
        /// <param name="data">data column can be double only</param>
        /// <param name="threshold"> set threshold; can be critical or range between value (<i>default will be 1.5</i>)</param>
        /// <returns>Return <b>index</b> and <b>values</b> which are considered as outliers in data array</returns>
        public Tuple<int[], double[]> PotentialOutliersIQR(double[] data, double threshold=1.5)
        {
            double Q1 = Quartile(data, 1);
            double Q3 = Quartile(data, 3);
            double IQR = Q1 - Q3; List<int> PotentialOutliers = new(); List<double> OutlierValues = new();
            for (int i = 0; i < data.Length; i++) {
                if (data[i] < Q1 - (threshold * IQR) || data[i] > Q3 - (threshold * IQR)) { 
                    PotentialOutliers.Add(i); OutlierValues.Add(data[i]);
                }
            }
            return Tuple.Create(PotentialOutliers.ToArray(), OutlierValues.ToArray());
        }

        /// <summary>
        /// Identify Potential Outliers in data array by Z-score mathod
        /// </summary>
        /// <param name="data"></param>
        /// <param name="threshold"></param>
        /// <returns>Return <b>index</b> and <b>values</b> which are considered as outliers in data array</returns>
        public Tuple<int[], double[]> PotentialOutliersZscore(double[] data, double threshold=3)
        {
            List<int> indxPO = new(); List<double> valPO = new();
            double[] zScore = Zscore(data);
            for (int i = 0;i < data.Length;i++)
            {
                if (zScore[i] > threshold || zScore[i] < -threshold) {
                    indxPO.Add(i); valPO.Add(data[i]);
                }
            }
            return Tuple.Create(indxPO.ToArray(), valPO.ToArray());
        }

        public double PearsonCorrelation(double[] X, double[] Y) 
        {
            double X_ = mean(X); double Y_ = mean(Y);
            double Numrator = 0; double Denomrator1 = 0; double Denomrator2 = 0;
            for (int i = 0; i < X.Length; i++) {
                Numrator += (X[i] - X_) * (Y[i] - Y_);
                Denomrator1 += Math.Pow(X[i] - X_, 2.0);
                Denomrator2 += Math.Pow(Y[i] - Y_, 2.0);
            }
            return Numrator / Math.Sqrt(Denomrator1*Denomrator2);
        }

        public double[,] CorrelationMatrix(double[,] data, string? WriteCSV=null, string[]? Annotation=null, double threshold=0.0)
        {                                                 // WriteCSV is the path of the output we want to save
            descriptiveAnalysis arrays = new(); csvReader csvArr = new();
            double[,] corr = arrays.Matrix(data.GetLength(1), data.GetLength(1));
            for (int i = 0; i < data.GetLength(1); i++)
            {
                for (int j = 0; j <= i; j++) // Avoid duplicate calculations (symmetric matrix)
                {
                    double correlation = PearsonCorrelation(csvArr.SelectCol(data, i), csvArr.SelectCol(data, j));
                    if (threshold == 0.0) {
                        corr[i, j] = correlation;
                        corr[j, i] = correlation; // Fill the symmetric counterpart
                    }
                    else {
                        if (correlation < threshold) {
                            corr[i, j] = 1;
                            corr[j, i] = 1; 
                        }
                        else {
                            corr[i, j] = 0;
                            corr[j, i] = 0;
                        }
                        
                    }
                    
                }
            }
            if (WriteCSV != null) {
                using (var csv = File.CreateText(WriteCSV+$"\\Correlation.csv"))
                {
                    if (Annotation != null) { csv.WriteLine(","+string.Join(",", Annotation)); }
                    for(int i = 0;i < corr.GetLength(0); i++) {
                        if (Annotation != null) { csv.Write(Annotation[i]+","); }
                        for (int j = 0;j < corr.GetLength(1); j++) {
                            csv.Write(corr[i, j]+",");
                        }
                        csv.WriteLine();
                    }
                }
            }
            return corr;
        }

        // Lanczos approximation. {using bard} ..
        public double GammaFunction(double z)
        {
            const double g = 5.1;
            double[] p = new double[] { 0.999999997385850393468, 5.88235294112859429497, -6.50593479999979434044,
                -0.00616504938802987295963, 0.0000700350306153087951438, -0.00000040082302615802308144,1.20881100000000227329,
            };
            double x = 1.0 + (z - 0.5); double a = 0.0;
            for (int i = 0; i < p.Length; i++)
            {
                a += p[i] / (x + i);
            }
            return Math.Exp(-z + g) * Math.Sqrt(2 * Math.PI) * a / x;
        }

        // {using bard} ..
        public double NormalCDF(double x)
        {
            // Standard normal cumulative distribution function (CDF)
            double z = (x - 0.0) / Math.Sqrt(1.0);
            double t = 1.0 / (1.0 + 0.5 * Math.Abs(z));

            // Abramowitz and Stegun approximation
            double n = t * Math.Exp(-z * z - 1.26551223 + t * (1.00002368 +
                           t * (0.37409196 +
                               t * (0.09678418 +
                                   t * (-0.18628806 +
                                       t * (0.27886807 +
                                           t * (-1.13520398 +
                                               t * (1.48851587 +
                                                   t * (-0.82215223 +
                                                       t * (0.17087277))))))))));

            return (1.0 + Math.Sign(z)) * (0.5 - n);
        }

        // {using Math.Net}
        public double TwoTailedPvalue(double t, double df)
        {
            double cdf;
            var t_distribution = new StudentT(0, 1, df); // using Math.Net
            if (t < 0)
            {
                cdf = t_distribution.CumulativeDistribution(Math.Abs(t));
            }
            else { cdf = t_distribution.CumulativeDistribution(t); }
            
            return 2 * Math.Min(1 - cdf, cdf);
        }
    }


    public class GEOAnalysis: Statistics
    {
        public struct TScore
        {
            public double TStatistics;
            public int DegreeOfFreedom;
            public double AvgExpression;
        }

        public TScore TStatistics(double[] array1, double[] array2)
        {
            double x1 = mean(array1); int n1 = array1.Length;
            double x2 = mean(array2); int n2 = array2.Length;
            double sqrtN = Math.Sqrt((Math.Pow(stdDev(array1), 2.0) / n1) + 
                (Math.Pow(stdDev(array2), 2.0) / n2));
            double sum = array1.Sum() + array2.Sum();
            var ts = new TScore
            {
                TStatistics = (x1 - x2) / sqrtN,
                AvgExpression = sum / n1 + n2,
                DegreeOfFreedom = (n1 + n2) - 2, // two sample with equal variances
            };
            return ts;
        }

        private double pooledStdDev(double[] array1, double[] array2)
        {
            int n2 = array2.Length;
            int n1 = array1.Length;
            double s1 = (n1-1)*stdDev(array1); 
            double s2 = (n2-1)*stdDev(array2);

            double stdN = Math.Pow(s1, 2.0) + Math.Pow(s2, 2.0);
            return Math.Sqrt(stdN/n1+n2-2);
        }

        private double[] ValueBetween(double[] arr, double val)
        {
            double[] tmp = new double[2];
            
            for (int i = 0; i<arr.Length-1; i++)
            {
                if (val < arr[i]){
                    tmp[1] = arr[i];
                    tmp[0] = arr[i-1];
                    break;
                }
            }
            return tmp;
        }

        public string Significance(double padjValue, double threshold = 0.05)
        {
            string sign = (padjValue <= threshold) ? "TRUE" : "FALSE";
            return sign;
        }

        public struct foldChanges
        {
            public double foldChange;
            public double log2FoldChange;
            public string regulation;
        }

        public foldChanges FoldChange(double[] S1, double[] S2)
        {
            double m1S = mean(S1); double m2S = mean(S2);
            double fc = (m2S - m1S) / m1S;
            string reg = (fc < 0)? "DOWN": "UP";

            var fC = new foldChanges
            {
                foldChange = fc,
                log2FoldChange = Math.Log2(m2S/m1S),
                regulation = reg
            };
            return fC;
        }

        public double[] PadjValue(double[] pvalue)
        {
            double[] sortp = pvalue.OrderBy(x => x).ToArray();
            int[] index = new int[sortp.Length];
            double[] padj = new double[sortp.Length];
            for(int i = 0; i < sortp.Length; i++)
            {
                index[i] = Array.IndexOf(pvalue, sortp[i]);
            }
            for(int i = 0;i < sortp.Length; i++)
            {
                padj[index[i]] = (sortp[i]*pvalue.Length)/(i + 1);
                // Console.WriteLine(index[i] +"\t"+(sortp[i] * pvalue.Length) / (i + 1));
            }
            return padj;
        }

        public struct GEOBatch
        {
            public double[] fc, ts, log2fc; public string[] regs;
        }

        private Tuple<int[], int[]> sampleSplit(int totalSample)
        {
            int half = totalSample / 2; int[] s1 = new int[half]; int[] s2 = new int[half];
            for (int i = 0; i < half; i++)
            {
                s1[i] = i; s2[i] = i + half;
            }
            return Tuple.Create(s1, s2);
        }

#pragma warning disable

        /// <summary>
        /// Function for Creating CSV file for <b>normalization</b> and <b>pre-processing</b>, 
        /// In pre-processing contains <b>foldChange, log2FC, regulation, t statistics, p value and adjusted p value</b> as column.
        /// </summary>
        /// <param name="Countcsv">Full path of CSV Count Table</param>
        /// <param name="TotalSamplesIndex">Array of the indexes of combined samples</param>
        /// <param name="OutPutLoc">Full path of the folder where you want to save your output</param>
        /// <param name="NormType">Normalization type can be "QN", "RPKM", "FPKM", "TCN" (default)</param>
        /// <param name="GeneLen">Length of the gene in bp; mandatory field for "FPKM" and "RPKM" normalization type</param>
        /// <param name="AvgFragmentLen">Average of the fragment length</param>
        /// <param name="accessionIndex">By default 0 can be index for accessionId</param>
        /// <return>FileName which is saved as desired Location of the folder</return>
        public void GEOParameters(string Countcsv, int[] TotalSamplesIndex, string OutPutLoc, string NormType="TCN",
            double[] GeneLen= null, double AvgFragmentLen=0, int accessionIndex=0)
        {
            GEOCSVTable gEOCSV = new GEOCSVTable();
            ExpressionNormalization normalization = new ExpressionNormalization();

            var sampleReads = gEOCSV.IndexSeparatedCounts(Countcsv, TotalSamplesIndex);
            double[,] Rawcounts = sampleReads.CountData; string[] accessionNos = sampleReads.Accession;
            double[,] counts = new double[Rawcounts.GetLength(0), Rawcounts.GetLength(1)];

            if (NormType == "QN") { counts = normalization.QuantileNormalization(Rawcounts); }
            else if (NormType == "RPKM")
            {
                counts = normalization.RPKMNormalization(Rawcounts, GeneLen);
            }
            else if (NormType == "FPKM")
            {
                counts = normalization.FPKMNormalization(Rawcounts, GeneLen, AvgFragmentLen);
            }
            else { counts = normalization.TotalCountNormalization(Rawcounts); }

            var SampleSize = sampleSplit(TotalSamplesIndex.Length);
            var samples = gEOCSV.TwoSampleClassification(counts, SampleSize.Item1, SampleSize.Item2);

            double[] fc = new double[counts.GetLength(0)]; double[] ts = new double[counts.GetLength(0)];
            double[] logfc = new double[counts.GetLength(0)]; string[] regs = new string[counts.GetLength(0)];
            double[] pval = new double[counts.GetLength(0)]; double[] adjP = new double[counts.GetLength(0)];
            double[] avg = new double[counts.GetLength(0)];

            int DegreeOfFreedom = (SampleSize.Item1.Length + SampleSize.Item2.Length) - 2;

            for (int row = 0; row < counts.GetLength(0); row++)
            {
                double[] s1 = new double[SampleSize.Item1.Length]; double[] s2 = new double[SampleSize.Item2.Length];
                for (int col = 0; col < SampleSize.Item1.Length; col++)
                {
                    s1[col] = samples.Sample1[row, col];
                    s2[col] = samples.Sample2[row, col];
                }
                double tsi = TStatistics(s1, s2).TStatistics;
                double avgi = TStatistics(s1, s2).AvgExpression;
                fc[row] = FoldChange(s1, s2).foldChange; logfc[row] = FoldChange(s1, s2).log2FoldChange;
                ts[row] = tsi; regs[row] = FoldChange(s1, s2).regulation;
                pval[row] = TwoTailedPvalue(tsi, DegreeOfFreedom); avg[row] = avgi;
            }
            double[] adjp = PadjValue(pval);

            using (var csv = File.CreateText(OutPutLoc+$"\\GEO_{NormType}Param.csv"))
            {
                string header = $"Accession,FoldChange,log2FoldChange,Regulation,MeanExpression,TStatistics,Pvalue,Adj. Pvalue,Sig.";
                csv.WriteLine(header);
                for (int i = 1; i < counts.GetLength(0)+1; i++)
                {
                    string line = $"{accessionNos[i - 1]},{fc[i - 1]},{logfc[i - 1]},{regs[i - 1]},{avg[i - 1]},{ts[i - 1]},{pval[i - 1]},{adjp[i - 1]},{Significance(adjp[i - 1])}";
                    csv.WriteLine(line);
                }
                csv.Close();
            }
            using (var csv = File.CreateText(OutPutLoc + $"\\GEO_{NormType}Normal.csv"))
            {
                csv.WriteLine(File.ReadAllLines(Countcsv).First());
                for(int i = 1;i< counts.GetLength(0)+1; i++)
                {
                    csv.Write(accessionNos[i - 1]+",");
                    for (int j = 0;j < counts.GetLength(1); j++)
                    {
                        csv.Write(counts[i-1,j]+",");
                    }
                    csv.WriteLine();
                }
                csv.Close();
            }
            Console.WriteLine($"{counts.GetLength(0)} Rows written on {OutPutLoc}\\GEO_{NormType}Normal.csv " +
                $"and Pre processing on {OutPutLoc}\\GEO_{NormType}Param.csv");
        }

        public double tDistributionPDF()
        {
            throw new NotImplementedException();
        }

        private double[] TwoTpValues = new double[] { 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001 };

        private double[][] T_Distribution = new double[][]
        {                                                                                   // df
            new double[]{ 3.078, 6.314, 12.706, 31.82, 63.657, 127.321, 318.309, 636.619 }, // 1
            new double[]{ 1.886, 2.92, 4.303, 6.965, 9.925, 14.089, 22.327, 31.599 },       // 2
            new double[]{ 1.638, 2.353, 3.182, 4.541, 5.841, 7.453, 10.215, 12.924 },       // 3
            new double[]{ 1.533, 2.132, 2.776, 3.747, 4.604, 5.598, 7.173, 8.61 },          // 4
            new double[]{ 1.476, 2.015, 2.571, 3.365, 4.032, 4.773, 5.893, 6.869 },         // 5
            new double[]{ 1.44, 1.943, 2.447, 3.143, 3.707, 4.317, 5.208, 5.959 },          // 6
            new double[]{ 1.415, 1.895, 2.365, 2.998, 3.499, 4.029, 4.785, 5.408 },         // 7
            new double[]{ 1.397, 1.86, 2.306, 2.897, 3.355, 3.833, 4.501, 5.041 },          // 8
            new double[]{ 1.383, 1.833, 2.262, 2.821, 3.25, 3.69, 4.297, 4.781 },           // 9
            new double[]{ 1.372, 1.812, 2.228, 2.764, 3.169, 3.581, 4.144, 4.587 },         // 10
        };

        // Approximation mathod using T distribution table, valid for only df value up to 10
        public double AproxPvalue(double tScore, int df)
        {
            df = df - 1;
            double pValue;
            if (T_Distribution[df].Contains(tScore))
            {
                int indexOftScore = Array.IndexOf(T_Distribution[df], tScore);
                pValue = TwoTpValues[indexOftScore];
            }
            else
            {
                double[] distributionValueRange = ValueBetween(T_Distribution[df], tScore);
                int indexOftRange1 = Array.IndexOf(T_Distribution[df], distributionValueRange[0]);
                int indexOftRange2 = Array.IndexOf(T_Distribution[df], distributionValueRange[1]);
                pValue = (TwoTpValues[indexOftRange1] + TwoTpValues[indexOftRange2]) / 2;
            }
            return pValue;
        }
    }

    /// <summary>
    /// Different Normalization methods like Total count, Quantile, RPKM, FPKM, normalization;
    /// Use for performing normalization on Gene expression counts of the different samples.
    /// </summary>
    public class ExpressionNormalization: descriptiveAnalysis
    {
        /// <summary>
        /// 
        /// </summary>
        /// <param name="mat">2D matrix of counted samples</param>
        /// <returns>Normalized 2D array</returns>
        public double[,] TotalCountNormalization(double[,] mat)
        {
            double[,] Normalize = new double[mat.GetLength(0), mat.GetLength(1)];
            for(int i = 0; i<mat.GetLength(0); i++)
            {
                double sum = 0;
                for (int j = 0; j<mat.GetLength(1); j++)
                {
                    sum += mat[i, j];
                }
                for(int j = 0; j<mat.GetLength(1); j++)
                {
                    Normalize[i, j] = mat[i, j]/sum;
                }
            }
            return Normalize;
        }

        public double[,] RPKMNormalization(double[,] mat, double[] geneLen)
        {
            double[,] RPKMNormalize = new double[mat.GetLength(0),mat.GetLength(1)];
            for(int i = 0; i<mat.GetLength(0); i++)
            {
                double sum = 0;
                for (int j = 0;j<mat.GetLength(1); j++)
                {
                    sum += mat[i, j];
                }
                for (int j = 0; j<mat.GetLength(1); j++)
                {
                    RPKMNormalize[i, j] = (mat[i, j] / geneLen[i])*1000000/(sum/1000000);
                }
            }
            return RPKMNormalize;
        }

        public double[,] FPKMNormalization(double[,] mat, double[] geneLen, double avgFragmentSize)
        {
            double[,] RPKMNormalize = new double[mat.GetLength(0), mat.GetLength(1)];
            for (int i = 0; i < mat.GetLength(0); i++)
            {
                double sum = 0;
                for (int j = 0; j < mat.GetLength(1); j++)
                {
                    sum += mat[i, j];
                }
                for (int j = 0; j < mat.GetLength(1); j++)
                {
                    RPKMNormalize[i, j] = ((mat[i, j] / geneLen[i]) * 1000000 / (sum / 1000000))*avgFragmentSize;
                }
            }
            return RPKMNormalize;
        }

        public double[,] QuantileNormalization(double[,] mat)
        {
            double[,] transp = Transpose(mat);
            int[,] IndexSorted = new int[transp.GetLength(0), transp.GetLength(1)];
            
            for (int eachRow = 0; eachRow < transp.GetLength(0); eachRow++)
            {
                double[] Row = new double[transp.GetLength(1)];
                for(int eachCol = 0;eachCol < transp.GetLength(1); eachCol++)
                {
                    Row[eachCol] = transp[eachRow, eachCol];
                }
                int[] IndexSort = ArgSort(Row.ToList());
                Array.Sort(Row);
                for (int eachCol = 0; eachCol<transp.GetLength(1); eachCol++)
                {
                    IndexSorted[eachRow, eachCol] = IndexSort[eachCol];
                    transp[eachRow, eachCol] = Row[eachCol];
                }
            }

            double[,] LastTransposed = Transpose(transp);
            double[] normArray = new double[LastTransposed.GetLength(0)];

            for (int i = 0; i < LastTransposed.GetLength(0); i++)
            {
                double row = 0;
                for(int j = 0; j< LastTransposed.GetLength(1); j++)
                {
                    row += LastTransposed[i, j];
                }
                double diff = row / LastTransposed.GetLength(1);
                normArray[i] = diff;
            }

            double[,] Normalize = new double[IndexSorted.GetLength(0),IndexSorted.GetLength(1)];

            for (int i=0; i<IndexSorted.GetLength(0); i++)
            {
                for(int j=0; j< IndexSorted.GetLength(1); j++)
                {
                    Normalize[i, j] = normArray[IndexSorted[i, j]];
                }
            }
            return Transpose(Normalize);
        }

    }

    // Bag of words is the simple way to vectorize string or sequence data
    public class BagOfWords
    {
        descriptiveAnalysis Unique = new descriptiveAnalysis();
        private string[] KmerAndUniqueModifire(string[] sequnce, int K){
            BioTools tool = new();
            List<string> String = new();

            for(int i=0; i<sequnce.Length; i++){
                var kmer = tool.K_merSequence(sequnce[i], K);
                foreach(string s in kmer){ String.Add(s); }
            }
            string[] uniques = Unique.UniqueValue(String.ToArray());
            return uniques;
        }

        public void WordCounter(string[] sequence, int K=6){
            string[] uniqueWords = KmerAndUniqueModifire(sequence, K);
            
            throw new NotImplementedException();
        }
    }
}