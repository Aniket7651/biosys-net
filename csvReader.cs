using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection.Metadata;
using System.Runtime.Remoting;
using System.Threading.Tasks;

namespace BioSySNet
{

    # pragma warning disable CS8602
    
    public struct TRAIN_TEST
    {
        public double[ , ] X_Train;
        public double[] Y_Train;
        public double[ , ] X_Test;
        public double[] Y_Test;
    } 


    public class ImportData
    {
        public List<List<object>> FromCommaDelim(string filePath)
        {
            List<List<object>> data = new();
            var reader = new StreamReader(filePath);
            while(!reader.EndOfStream) 
            {
                List<object> temp = new();
                var line = reader.ReadLine();
                foreach (string str in line.Split(","))
                {
                    if (double.TryParse(str, out _)) { temp.Add(double.Parse(str)); }
                    else { temp.Add(str); }
                }
                data.Add(temp);
                
            }
            return data;
        }

        public List<List<object>> FromTabDelim(string filePath)
        {
            List<List<object>> data = new();
            var reader = new StreamReader(filePath);
            while(!reader.EndOfStream) 
            {
                List<object> temp = new();
                var line = reader.ReadLine();
                foreach (string str in line.Split("\t"))
                {
                    if (double.TryParse(str, out _)) { temp.Add(double.Parse(str)); }
                    else { temp.Add(str); }
                }
                data.Add(temp);
                
            }
            return data;
        }

        public List<List<object>> FromOtherDelim(string filePath, string sep)
        {
            List<List<object>> data = new();
            var reader = new StreamReader(filePath);
            while(!reader.EndOfStream) 
            {
                List<object> temp = new();
                var line = reader.ReadLine();
                foreach (string str in line.Split(sep))
                {
                    if (double.TryParse(str, out _)) { temp.Add(double.Parse(str)); }
                    else { temp.Add(str); }
                }
                data.Add(temp);
                
            }
            return data;
        }
    }

    public class DataFrame
    {
        public List<List<object>> Data;
        public int[] Shape;  // ROWS and COLUMNS are not updated by changes applied in dataset shape
        public int Rows, Columns;
        public List<object>? HEADER;
    
        public DataFrame(List<List<object>> ListData, bool header = true, 
            bool displayData = false, bool defaultTypes = false)
        {
            List<List<object>> temp = new List<List<object>>();
            for (int item = 1; item < ListData.Count; item++)
            {
                if (defaultTypes) { temp.Add(ListData[item]); }
                else { 
                    List<object> t1 = new();
                    foreach (object str in ListData[item])
                    {
                        try{
                        if (double.TryParse((string)str, out _)){ t1.Add(double.Parse((string)str)); }
                        else { t1.Add((string)str); }
                        }
                        catch (InvalidCastException){ Console.WriteLine(
                            "ERROR IN TYPE CONVERSION: Change the defaultType=true OR false..."
                        ); break; }
                    }
                    temp.Add(t1);
                }
            }
            Data = header? temp: ListData;

            Shape = new int[2]{Data.Count, Data[0].Count};
            Rows = Data.Count; Columns = Data[0].Count;
            HEADER = header? ListData[0]: Data[0];
            if (displayData == true)
            {

                for(int d = 0; d <= Math.Abs(Data.Count/4); d++)
                {
                    foreach (object value in Data[d])
                    {
                        Console.Write(value+"\t\t");
                    }
                    Console.WriteLine();
                }
            } 
        }

        public void Head(int upto = 5)
        {
            if (upto > Rows)
            {
                for(int d = 0; d <= Data.Count; d++)
                {
                    foreach (object value in Data[d])
                    {
                        Console.Write(value+"\t\t");
                    }
                    Console.WriteLine();
                }
            }
            else
            {
                for(int d = 0; d <= upto; d++)
                {
                    foreach (object value in Data[d])
                    {
                        Console.Write(value+"\t\t");
                    }
                    Console.WriteLine();
                }
            }
        }    // NOT TESTED

        public List<List<object>> TransposeFrame(){
        // reference:- https://stackoverflow.com/questions/39484996/rotate-transposing-a-listliststring-using-linq-c-sharp
            List<List<object>> res = Data
                .SelectMany(inner => inner.Select((item, index) => new { item, index }))
                .GroupBy(i => i.index, i => i.item)
                .Select(g => g.ToList())
                .ToList();

            return res;
        }

        public List<List<object>> SubSet(params string[] colName)
        {
            int[] listIndex = new int[colName.Length];
            for (int i = 0; i < colName.Length; i++)
            {
                listIndex[i] = HEADER.IndexOf(colName[i]);
            }
            DataFrame frame = new(Data, header: false);
            return frame.ISubSet(listIndex);
        }

        public List<List<object>> ISubSet(params int[] colIndex)
        {
            List<List<object>> Out_temp = new List<List<object>>();
            for (int i = 0; i < colIndex.Length; i++)
            {
                List<object> In_temp = new List<object>();
                for (int j = 0; j < Rows; j++)
                {
                    In_temp.Add(Data[j][colIndex[i]]);
                } Out_temp.Add(In_temp);
            }
            DataFrame tempTranspose = new(Out_temp, header: false);
            return tempTranspose.TransposeFrame();
        }

        public List<List<object>> ArrayToList(object[ , ] DataFrameInObj)
        {
            List<List<object>> list = new List<List<object>>();
            for (int i = 0; i < DataFrameInObj.GetLength(0); i++)
            {
                List<object> temp = new();
                for (int j = 0; j < DataFrameInObj.GetLength(1); j++) 
                {
                    temp.Add(DataFrameInObj[i, j]);
                }
                list.Add(temp);
            }
            return list;
        }

        public void RenameIColumn(int colIndex, string newName)
        {
            HEADER[colIndex] = newName;
        }

        public void RenameColumn(string colName, string newName)
        {
            int Ind = HEADER.IndexOf(colName);
            HEADER[Ind] = newName;
        }

        public void AddColumn(List<object> colValues, string colName, int IndexAt)
        {
            if (Columns >= IndexAt)
            {
                HEADER.Insert(IndexAt, colName);
                for (int i = 0; i < Rows; i++)
                {
                    Data[i].Insert(IndexAt, colValues[i]);
                }
            }
            else
            {
                AddColumnToLast(colValues, colName);
            }
        }

        public void AddColumnToLast(List<object> colValues, string colName)
        {
            HEADER.Insert(Columns, colName);
            for (int i = 0; i < Rows; i++)
            {
                Data[i].Insert(Columns, colValues[i]);
            }
        }

        public List<List<object>> Stats(bool display = false, string? saveLogs = null)
        {
            Statistics measures = new Statistics();
            List<List<object>> transpose = TransposeFrame();
            double[] sum = new double[transpose.Count];
            double[] max = new double[transpose.Count];
            double[] min = new double[transpose.Count];
            double[] mean = new double[transpose.Count];
            double[] Q1 = new double[transpose.Count];
            double[] Q3 = new double[transpose.Count];
            double[] median = new double[transpose.Count];
            double[] stdDev = new double[transpose.Count];
            
            for (int i = 0; i < transpose.Count; i++)
            {
                List<object> TempCols = transpose[i];
                if (TempCols[1] is not string)
                {
                    object[] TempArray = TempCols.ToArray();
                    double[] TempColsArray = Array.ConvertAll(TempArray, x => (double)x);

                    max[i] = TempColsArray.Max(); min[i] = TempColsArray.Min();
                    double sumtemp = 0;
                    foreach (object val in TempColsArray) { sumtemp += (double)val; }
                    sum[i] = sumtemp;
                    mean[i] = measures.mean(TempColsArray);
                    median[i] = measures.median(TempColsArray);
                    stdDev[i] = measures.stdDev(TempColsArray);
                    Q1[i] = measures.Quartile(TempColsArray, 1);
                    Q3[i] = measures.Quartile(TempColsArray, 3);
                }
                else
                {
                    sum[i] = max[i] = min[i] = mean[i] = median[i] = stdDev[i] = Q1[i] = Q3[i] = 00;
                }
            }
            int v = 0;
            List<List<object>> Matrix = new(); 
            List<object> head = new List<object>(){
                "Col", "Sum", "Max", "Min", "Mean", "Median", "Std", "Q1(25%)", "Q3(75%)"
            };
            Matrix.Add(head);
            // if (display) { Console.WriteLine("Col\t\tSum\t\tMax\t\tMin\t\tMean\t\tMed\t\tStd\t\tQ1_\t\tQ3_"); }
            do {
                List<object> values = new(){HEADER[v], sum[v], max[v], min[v], mean[v], median[v],
                    stdDev[v], Q1[v], Q3[v]};
                if (display) {
                    Console.WriteLine(HEADER[v]+"\t\t"+sum[v]+"\t\t"+max[v]+"\t\t"+min[v]+"\t\t"+mean[v]+
                        "\t\t"+median[v]+"\t\t"+stdDev[v]+"\t\t"+Q1[v]+"\t\t"+Q3[v]);
                }
                v++;
                Matrix.Add(values);
            } while (v < transpose.Count);
            if (saveLogs != null) { 
                using (var file = File.CreateText(saveLogs+"Info(Data).csv"))
                {
                    for (int i = 0; i < Matrix.Count; i++)
                    {
                        for (int j = 0; j < Matrix[0].Count; j++){ file.Write(Matrix[i][j]+","); }
                        file.WriteLine();
                    }
                    file.Close();
                }
            }
            return Matrix;
        }

        public List<object> IColumn(int idx)
        {
            List<object> Lis = new();
            for (int i = 0; i < Rows; i++)
            {
                Lis.Add(Data[i][idx]);
            }
            return Lis;
        }

        public object[ , ] ListToArray() 
        {
            object [ , ] Arrays = new object[Rows, Columns];
            for (int i = 0; i < Rows; i++)
            {
                for (int j = 0; j < Columns; j++) 
                {
                    Arrays[i, j] = Data[i][j];
                }
            }
            return Arrays;
        }

        private int[] SortedIndxInCol(List<object> Col) {
            descriptiveAnalysis dupl = new();
            
            List<object> TargetSorted = Col.OrderBy(o=>o).ToList();
            List<object> uniq = dupl.RemoveDuplicate(TargetSorted);
            List<int> SortIndx = new();
            for(int i = 0; i < Col.Count; i++)
            {
                int Ind = TargetSorted.IndexOf(Col[i]);
                    
                if (SortIndx.Contains(Ind)) {
                    SortIndx.Add(Ind+1);
                }
                else { SortIndx.Add(Ind); }
            }
            return SortIndx.ToArray();
            
        } // CHECKED 

        public List<List<object>> SortFrame(int ByColIndx)
        {
            object [ , ] result = new object[Rows, Columns];
            object[ , ] DataArray = ListToArray();
            int[] SortIndx = SortedIndxInCol(IColumn(ByColIndx));

            for (int i = 0; i < Rows; i++) 
            {
                for (int j = 0; j < Columns; j++) 
                {
                    result[SortIndx[i], j] = DataArray[i, j];
                }
            }
            return ArrayToList(result);
        }
        
        public void Info()
        {
            for (int j = 0; j < Columns; j++)
            {
                Console.WriteLine(HEADER[j]+"\t\t"+Data[0][j].GetType());
            }
        }

        public void Merge(List<List<object>> AnotherFrame, string ByColName)
        {
            DataFrame f1 = new(AnotherFrame, defaultTypes: true);
            int col = HEADER.IndexOf(ByColName);
            List<object> T1 = IColumn(col); List<object> T2 = f1.IColumn(col);
            
            for (int i = 0; i < AnotherFrame.Count-1; i++) {
                int Idx = T1.IndexOf(T2[i]);
                AnotherFrame[i+1].RemoveAt(col);
                Data[Idx].AddRange(AnotherFrame[i+1]);
            }
            f1.HEADER.RemoveAt(col);
            HEADER.AddRange(f1.HEADER);
        }

        public void SaveAsCSV(string fileLoc, string sep = ",")
        {
            using(var csv = File.CreateText(fileLoc))
            {
                foreach (object head in HEADER) { csv.Write(head+sep); }
                csv.WriteLine();
                for (int i = 0; i < Rows; i++)
                {
                    for (int j = 0; j < Columns; j++){ csv.Write(Data[i][j]+sep); }
                        csv.WriteLine();
                }
                csv.Close();
            }
        }
    }

    public class CSVArray
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
        
        public double[ , ] readNumericCSV(string location){
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
            double[ , ] csvarray = readNumericCSV(location);
            double[ , ] seletedRowArr = new double[top, csvarray.GetLength(1)];
            for(int i=0; i<top; i++){
                for(int j=0; j<csvarray.GetLength(1); j++){
                    seletedRowArr[i,j] = csvarray[i,j];
                }
            }
            return seletedRowArr;
        }

        public TRAIN_TEST TrainTest_Split(double[,] datasetArray, double test_size=0.2){
            int row = datasetArray.GetLength(0);
            int col = datasetArray.GetLength(1); // if, total = 10
            int testRows = (int)(test_size*row); // then, test = 2
            int trainRows = (int)(row-testRows); // then, train = 8

            var DataFrames = new TRAIN_TEST{
                X_Train = new double[trainRows, col-1],
                Y_Train = new double[trainRows],
                X_Test = new double[testRows, col-1],
                Y_Test = new double[testRows]
            };
            // adding data in Train set
            for(int i=0; i<trainRows; i++){
                for(int j=0; j<col-1; j++){
                    DataFrames.X_Train[i, j] = datasetArray[i, j];  
                }
                double y = datasetArray[i, col-1];
                DataFrames.Y_Train[i] = y;
            }
            // adding data in Test set
            for(int i=row-1; i>trainRows-1; i--){
                for(int j=0; j<col-1; j++){
                    DataFrames.X_Test[(row-1)-i, j] = datasetArray[i, j]; 
                }
                double y = datasetArray[i, col-1];
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


    public class GEOCSVTable: CSVArray
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