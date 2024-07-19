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

    /// <summary>
    /// Load delimited files as dataframe of List Of List (LoL)
    /// </summary>
    public class ImportData
    {
        /// <summary>
        /// Read comma seperated files (CSVs)
        /// </summary>
        /// <param name="filePath">Full path of the CSV file</param>
        /// <returns>List Of List as dataframe</returns>
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

        /// <summary>
        /// Read Tab seperated Files (TSVs)
        /// </summary>
        /// <param name="filePath">Full path of the TSV or text file</param>
        /// <returns>List Of List as dataframe</returns>
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

        /// <summary>
        /// Read files, seperated with any symbol as delimition
        /// </summary>
        /// <param name="filePath">Full path of the CSV file</param>
        /// <param name="sep">Symbol seperated with</param>
        /// <returns>List Of List as dataframe</returns>
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

    /// <summary>
    /// DataFrame class provides for manipulating dataset and their operation.
    /// In the form of <b>List Of List (LoL)</b> and data type as <b>object</b>
    /// </summary>
    public class DataFrame
    {
        /// <summary>Load data as List Of List</summary>
        public List<List<object>> Data;
        public int[] Shape;  // ROWS and COLUMNS are not updated by every changes applied in dataset shape
        /// <summary>Load row and column length as integer type</summary>
        public int Rows, Columns;
        /// <summary>Load Header of the dataframe</summary>
        public List<object>? HEADER;

        /// <param name="ListData">Data frame in the form of List of List, where dType would be object</param>
        /// <param name="header">A bool type parameter for confirming if header exist on first row of dataset</param>
        /// <param name="displayData">Bool typed parameter, display data after it was loaded</param>
        /// <param name="defaultTypes">Bool type parameter for leave parameter as default, likly no need to convert types.</param>
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

        /// <summary>
        /// Representative, Head of the Dataset upto desired rows which returns void
        /// </summary>
        /// <param name="upto">Top Number of rows, where default value will 5</param>
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

        /// <summary>
        /// Transpose the dataframe which makes row to column and column to row
        /// </summary>
        /// <returns>List Of List (LoL) type can store as a variable</returns>
        public List<List<object>> TransposeFrame(){
        // reference:- https://stackoverflow.com/questions/39484996/rotate-transposing-a-listliststring-using-linq-c-sharp
            // as `LINQ method`
            List<List<object>> res = Data
                .SelectMany(inner => inner.Select((item, index) => new { item, index }))
                .GroupBy(i => i.index, i => i.item)
                .Select(g => g.ToList())
                .ToList();

            return res;
        }

        /// <summary>
        /// Create a new dataframe from existing dataframe with <b>columns name</b>
        /// </summary>
        /// <param name="colName">parameter as name of the multiple columns </param>
        /// <returns>List Of List as dataframe</returns>
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

        /// <summary>
        /// Create a new dataframe from existing dataframe with <b>columns Index</b>
        /// </summary>
        /// <param name="colIndex">Index of the multiple columns</param>
        /// <returns>List Of List as dataframe</returns>
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

        /// <summary>
        /// Convert form of <c>object[ , ]</c> type dataframe to List Of List dataframe
        /// </summary>
        /// <param name="DataFrameInObj"><c>object[ , ]</c> type dataframe</param>
        /// <returns>List Of List as dataframe</returns>
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

        /// <summary>
        /// Rename a column header By <b>Index</b>
        /// </summary>
        /// <param name="colIndex">Index of the column</param>
        /// <param name="newName">New name of the column to be replaced by existing one</param>
        public void RenameIColumn(int colIndex, string newName)
        {
            HEADER[colIndex] = newName;
        }

        /// <summary>
        /// Rename a column header By <b>Name of the column</b>
        /// </summary>
        /// <param name="colIndex">name of the column</param>
        /// <param name="newName">New name of the column to be replaced by existing one</param>
        public void RenameColumn(string colName, string newName)
        {
            int Ind = HEADER.IndexOf(colName);
            HEADER[Ind] = newName;
        }

        /// <summary>
        /// Add a column at the specific index in the dataframe
        /// </summary>
        /// <param name="colValues">List of object type is use as the values of that column</param>
        /// <param name="colName">Name of the new column</param>
        /// <param name="IndexAt">Index number where you want to add this new column</param>
        public void AddColumn(List<object> colValues, string colName, int IndexAt)
        {
            if (Columns >= IndexAt)
            { // if Index is lower OR equal to length of total column
                HEADER.Insert(IndexAt, colName);
                for (int i = 0; i < Rows; i++)
                {
                    Data[i].Insert(IndexAt, colValues[i]);
                }
            }
            else
            {
                // By default add column at the last
                AddColumnToLast(colValues, colName);
            }
        }

        /// <summary>
        /// Add new column at the last of the dataframe
        /// </summary>
        /// <param name="colValues">List of object type is use as the values of that column</param>
        /// <param name="colName">Name of the new column</param>
        public void AddColumnToLast(List<object> colValues, string colName)
        {
            HEADER.Insert(Columns, colName);
            for (int i = 0; i < Rows; i++)
            {
                Data[i].Insert(Columns, colValues[i]);
            }
        }

        /// <summary>
        /// Sum of single numeric row (axis = 1)
        /// </summary>
        /// <param name="row"></param>
        /// <returns>Sum of the row as double</returns>
        public double SumOfRow(List<object> row)
        {
            double Sumation = 0;
            foreach (object v in row)
            {
                Sumation += (double)v;
            }
            return Sumation;
        }
        
        /// <summary>
        /// Statistical Information of the numeric dataframe i.e. <i>Mean</i>, <i>Median</i>, <i>Maximum</i>
        /// <i>Minimum</i>, <i>Sum</i>, <i>Q1 quartile</i>, <i>Q3 quartile</i>, <i>standared deviation</i> like information.
        /// </summary>
        /// <param name="display">Display information frame by display=true</param>
        /// <param name="saveLogs">Path of the folder, where you want to save this information as CSV</param>
        /// <returns>List Of List as dataframe</returns>
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

        /// <summary>
        /// Convert form of List Of List dataframe to <c>object[ , ]</c> type dataframe
        /// </summary>
        /// <returns>List Of List as dataframe</returns>
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

        /// <summary>Status: CHECKED</summary>
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
        }

        /// <summary>
        /// Sort dataframe through a specific column index
        /// </summary>
        /// <param name="ByColIndx">Index of the column</param>
        /// <returns>Sorted dataframe of List Of List form</returns>
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
        
        /// <summary>
        /// Data type related Information of the dataframe, print on console
        /// </summary>
        public void Info()
        {
            for (int j = 0; j < Columns; j++)
            {
                Console.WriteLine(HEADER[j]+"\t\t"+Data[0][j].GetType());
            }
        }

        

        /// <summary>
        /// Merge OR Map two dataframe into a single dataframe
        /// </summary>
        /// <param name="AnotherFrame">Another dataframe to be merge from existing dataframe</param>
        /// <param name="ByColName">Merge dataframe by a specific column</param>
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

        /// <summary>
        /// Export dataframe as file
        /// </summary>
        /// <param name="fileLoc">Full path of that file OR name of the file where you want to save dataframe</param>
        /// <param name="sep">seperation of data values to each other, by default it's " , " (Comma seperated)</param>
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

    /// <summary>
    /// Something also same as <c>ImportData</c> class, But <c>CSVArray</c> Reads only numeric CSV file with
    /// numeric return type also in 1D OR 2D (Classical functionality of <c>DataFrame</c> and <c>ImportData</c> classes)
    /// </summary>
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

    /// <summary>
    /// A class for loading and handling RNA Sequencing Counts dataset, which is inherit from <c>CSVArray</c>
    /// </summary>
    public class GEOCSVTable: CSVArray
    {
        public struct GEOData
        {
            public string[]? Accession, SampleID; public int NumberOfReads;
            // public double[]? GeneLength, StartIndexes, EndIndexes;
            public double[,]? CountData;
        }

        public struct Samples { public double[,] Sample1; public double[,] Sample2; }

        /// <summary>
        /// A method to seperate counts, accession IDs and other information from RNA Counts file
        /// </summary>
        /// <param name="csvLocation">path of the RNA sequence count file</param>
        /// <param name="SampleCountsIndex">Array of that all Indexes</param>
        /// <param name="geneAccessionIndex">Index of the Accession Id where exist in the count dataset, default as 0</param>
        /// <returns>Multiple return type values as, <c>Accession</c>, <c>SampleID</c> and <c>NumberOfReads</c></returns>
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