using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace BioSySNet
{

    # pragma warning disable 8602
    
    public struct DataFrame
    {
        public float[ , ] X_Train;
        public float[] Y_Train;
        public float[ , ] X_Test;
        public float[] Y_Test;
    } 


    public class csvReader
    {
        public int[] CSVshape(string location){
            /*
                int[] (1D array) return type of function for get the shape of CSV file, rows and Columns
                These 1D array returns two value of row and column of a CSV file like:
                int[] shape = {row, col};
                takes the input which is the location of csv file as string variable.
            */
            var reader = new StreamReader(location);
            int row = 0;
            int colL = 0;
            while(!reader.EndOfStream){
                var line = reader.ReadLine();
                string[] col = line.Split(',').ToArray<string>();
                colL = col.Length;
                row++;
            }
            int[] shape = {row, colL};
            return shape;
        }
        
        public double[ , ] readCSV(string location){
            int[] shape = CSVshape(location);
            // declare 2D array                 R         C
            double[ , ] TwoDArr = new double[shape[0], shape[1]];

            var reader = new StreamReader(location);
            int i = 0;
            while(!reader.EndOfStream){
                if(i<shape[0]){
                    var line = reader.ReadLine();
                    var String_value_array = line.Split(',');

                    for(int j=0; j<shape[1]; j++){
                        TwoDArr[i, j] = double.Parse(String_value_array[j]);
                    }
                } 
                i++;
            }
            return TwoDArr;
        }

        public double[] SelectCol(double[,] matrix, int colIndex)
        {
            double[] SelectedCol = new double[matrix.GetLength(0)];
            for(int i = 0; i < matrix.GetLength(0); i++)
            {
                SelectedCol[i] = matrix[i, colIndex];
            }
            return SelectedCol;
        }

        public double[][] SelectMultiCols(double[,] matrix, int[] colIndex)
        {
            double[][] SelectedCols = new double[colIndex.Length][];
            for(int i = 0;i < colIndex.Length; i++)
            {
                double[] col = new double[matrix.GetLength(0)];
                for (int j = 0; j < matrix.GetLength(0); j++)
                {
                    col[j] = matrix[j, colIndex[i]];
                }
                SelectedCols[i] = col;
            }
            return SelectedCols;
        }
        
        public double[ , ] CSV_TopRows(string location, int top=5){
            double[ , ] csvarray = readCSV(location);
            double[ , ] seletedRowArr = new double[top, csvarray.GetLength(1)];
            for(int i=0; i<top; i++){
                for(int j=0; j<csvarray.GetLength(1); j++){
                    seletedRowArr[i,j] = csvarray[i,j];
                }
            }
            return seletedRowArr;
        }

        public DataFrame TrainTest_Split(float[,] datasetArray, double test_size=0.2){
            int row = datasetArray.GetLength(0);
            int col = datasetArray.GetLength(1); // if, total = 10
            int testRows = (int)(test_size*row); // then, test = 2
            int trainRows = (int)(row-testRows); // then, train = 8

            var DataFrames = new DataFrame{
                X_Train = new float[trainRows, col-1],
                Y_Train = new float[trainRows],
                X_Test = new float[testRows, col-1],
                Y_Test = new float[testRows]
            };
            // adding data in Train set
            for(int i=0; i<trainRows; i++){
                for(int j=0; j<col-1; j++){
                    DataFrames.X_Train[i, j] = datasetArray[i, j];  
                }
                float y = datasetArray[i, col-1];
                DataFrames.Y_Train[i] = y;
            }
            // adding data in Test set
            for(int i=row-1; i>trainRows-1; i--){
                for(int j=0; j<col-1; j++){
                    DataFrames.X_Test[(row-1)-i, j] = datasetArray[i, j]; 
                }
                float y = datasetArray[i, col-1];
                DataFrames.Y_Test[(row-1)-i] = y;
            }

            return DataFrames;
        }

        public List<object> RowSelection<Type>(Type[,] matrix, int rowIndex)
        {
            List<object> list = new List<object>();
            for (int j = 0;  j < matrix.GetLength(1); j++)
            {
                # pragma warning disable CS8604
                list.Add(matrix[rowIndex, j]);
            }
            return list;
        }
        
        public void visualizeCSV<Type>(Type[ , ] csvArray, bool dropIndexing=false){
            for (int i=0; i<csvArray.GetLength(0); i++){
                for(int j=0; j<csvArray.GetLength(1); j++){
                    Console.Write(csvArray[i,j]+" | ");
                }
                if (dropIndexing == true){
                    Console.WriteLine();
                }else{ Console.WriteLine(i); }
            }
        }
    }


    public class GEOCSVTable: csvReader
    {
        public struct GEOData
        {
            public string[]? Accession, SampleID; public int NumberOfReads;
            // public double[]? GeneLength, StartIndexes, EndIndexes;
            public double[,]? CountData;
        }

        public struct Samples { public double[,] Sample1; public double[,] Sample2; }

        public GEOData IndexSeparatedCounts(string csvLocation, int[] SampleCountsIndex, int geneAccessionIndex = 0) 
        {
            string[] reads = File.ReadAllLines(csvLocation);
            
            List<string> Accession = new List<string>();
            double[,] Counts = new double[reads.Length-1, SampleCountsIndex.Length];

            for (int r = 1; r < reads.Length; r++)
            {
                string[] row = reads[r].Split(',');
                Accession.Add(row[geneAccessionIndex]);

                for(int c = 0; c < SampleCountsIndex.Length; c++)
                {
                    Counts[r - 1, c] = double.Parse(row[SampleCountsIndex[c]]);
                }
            }
            string[] header = reads[0].Split(',').Select(ele => ele.ToLower()).ToArray();
            string[] sampl = new string[SampleCountsIndex.Length];

            for(int i = 0;  i < SampleCountsIndex.Length; i++) {
                sampl[i] = header[SampleCountsIndex[i]];
            }

            GEOData data = new GEOData
            {
                Accession = Accession.ToArray(), NumberOfReads = reads.Length,
                CountData = Counts, SampleID = sampl,
            };
            return data;
        }

        public double[] SelectGeneLen(string csvLocation)
        {
            string[] reads = File.ReadAllLines(csvLocation);
            string[] header = reads[0].Split(',').Select(ele => ele.ToLower()).ToArray();
            return SelectColumn(csvLocation, Array.IndexOf(header, "length"));
        }

        public string[] SelectColumnByName(string csvLocation, string colName)
        {
            string[] reads = File.ReadAllLines(csvLocation);
            string[] header = reads[0].Split(',');
            string[] vals = new string[reads.Length - 1];

            int lenindx = Array.IndexOf(header, colName);
            for (int r = 1; r < reads.Length; r++)
            {
                string[] row = reads[r].Split(',');
                vals[r - 1] = row[lenindx];
            }
            return vals;
        }

        public double[] SelectColumn(string csvLocation, int Index)
        {
            string[] reads = File.ReadAllLines(csvLocation);
            double[] vals = new double[reads.Length - 1];

            for (int r = 1; r < reads.Length; r++)
            {
                string[] row = reads[r].Split(',');
                vals[r - 1] = double.Parse(row[Index]);
            }
            return vals;
        }

        public string GEOHeader(string csvLocation)
        {
            return File.ReadAllLines(csvLocation).First();
        }

        public Samples TwoSampleClassification(double[,] CountTableData, int[] sampleCols1, int[] sampleCols2)
        {
            double[,] s1 = new double[CountTableData.GetLength(0),sampleCols1.Length];
            double[,] s2 = new double[CountTableData.GetLength(0), sampleCols2.Length];

            for (int i = 0; i < CountTableData.GetLength(0); i++)
            {
                for (int k = 0; k < sampleCols1.Length; k++)
                {
                    s1[i, k] = CountTableData[i, sampleCols1[k]];
                    s2[i, k] = CountTableData[i, sampleCols2[k]];
                }
                
            }
            var s12 = new Samples { Sample1 = s1, Sample2 = s2 };
            return s12;
        }


    }
}