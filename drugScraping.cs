/*
~drugScraping.cs
class file for providing functionality from drugbank database [https://go.drugbank.com/]
by use of HtmlAgilityPack nuget package we can scrap the data from drugbank site
such informations can reterive by this script like, other ids, targets and enzymes, chem. identifires

@ Author: ANIKET YADAV (aniketyadav8687@gmail.com)
*/

using System;
using System.Collections.Generic;
using System.Linq;
using HtmlAgilityPack;
using System.Diagnostics;

namespace BioSySNet
{
    public class drugScraping
    {
        // <a class="button tertiary" href="https://rest.uniprot.org/uniprotkb/P05067.fasta">Download</a>

        public void UniProtSeq(string entryNo){
            Process.Start(new ProcessStartInfo() { 
                FileName=$"https://rest.uniprot.org/uniprotkb/{entryNo}.fasta", 
                UseShellExecute = true
            });
        }

        public string drugId(string DrugName){
            var webCon = new HtmlWeb();
            var doc = webCon.Load($"https://go.drugbank.com/unearth/q?query={DrugName}&button=&searcher=drugs");
            var node = doc.DocumentNode.SelectSingleNode("/html/body/main/div/div/div[2]/div[2]/dl[1]/dd[4]");
            // foreach(HtmlNode item in node){ Console.WriteLine(item.InnerHtml); }
            return node.InnerHtml;
        }

        // bugyy ..........................
        public void fastaFile(string id){
            var webCon = new HtmlWeb();
            var doc = webCon.Load($"https://www.ncbi.nlm.nih.gov/nuccore/OQ622003.1?report=genbank&log$=seqview");
            HtmlNodeCollection node = doc.DocumentNode.SelectNodes(@"//span[@class='ff_line']");
            // foreach(HtmlNode item in node){ Console.WriteLine(item.InnerHtml); }
            //*[@id="viewercontent1"]/div/div/pre //span[@class="ff_line"]
            
            foreach(var n in node){
                Console.WriteLine(n);
            }
            
        }

        public struct Info{
            public object IUPAC, InChIKey, SMILE, CAS;
            public object ChEMBL, ChEBI, PubChemCompound, ZINC, KEGG;
        }

        public Info drugChemIdentification(string drugbank_id){
            var webCon = new HtmlWeb();
            var doc = webCon.Load($"https://go.drugbank.com/drugs/{drugbank_id}");
            var druginfo = new Info{
                IUPAC = doc.DocumentNode.SelectSingleNode("/html/body/main/div/div/div[2]/div[2]/dl[6]/dd[5]/div").InnerHtml,
                SMILE = doc.DocumentNode.SelectSingleNode("/html/body/main/div/div/div[2]/div[2]/dl[6]/dd[6]/div").InnerHtml,
                CAS = doc.DocumentNode.SelectSingleNode("/html/body/main/div/div/div[2]/div[2]/dl[6]/dd[2]").InnerHtml,
                InChIKey = doc.DocumentNode.SelectSingleNode("/html/body/main/div/div/div[2]/div[2]/dl[6]/dd[3]").InnerHtml,
            };
            return druginfo;
        }

        // public Info otherDrugId(string drugbank_id){

        // }
    }

    public class NCBI
    {
        public void GetUniProtID(string geneName, string Organism)
        {
            var webCon = new HtmlWeb();
            var doc = webCon.Load($"");
            var node = doc.DocumentNode.SelectSingleNode("/html/body/main/div/div/div[2]/div[2]/dl[1]/dd[4]");
            Console.WriteLine(node.InnerText);
            // return node.InnerText;
        }
    }
}