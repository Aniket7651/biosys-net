# <code style="font-size: 30px">C# For Bioinformatics and Biological Data Processing: BioSyS.Net</code>

<p align="center" width="100%">
    <a href='https://www.nuget.org/packages/BioSySNet'><img src="removebg-logo.ico" alt="BioSyS_logo"></a>
</p>

![net](https://img.shields.io/badge/.NET-512BD4?style=for-the-badge&logo=dotnet&logoColor=white)
![nuget](https://img.shields.io/badge/NuGet-004880?style=for-the-badge&logo=nuget&logoColor=white)
![cs](https://img.shields.io/badge/C%23-239120?style=for-the-badge&logo=csharp&logoColor=white)
![python](https://img.shields.io/badge/Python-FFD43B?style=for-the-badge&logo=python&logoColor=blue)

`Bioinformatics` is an interdisciplinary field that analyzes biological data using computer science. I have created a new package for bioinformatics and computational analysis built in C#. This package can analyze biological data such as DNA, RNA, and protein sequences, as well as perform tasks such as Biological Database Scraping. Additionally, `'GEO Analysis'` offers a variety of analyses for RNA-seq preprocessing, parsing various biological file formats (FASTA, PDB, SOFT, FASTQ, etc.). [BioSySNet](https://www.nuget.org/packages/BioSySNet) also utilizes the Python interpreter to run certain plotting functions (i.e., it uses IronPython to execute Python functions, "so first, you need to make sure that Python is installed on your machine").

# ` Features `
- Scrape useful information from various databases like, DrugBank, UnipProt, PDB etcand PDB.
- Download and parse some of the important biological file formats such as .pdb, .fasta, .fastq, .soft, and CSV.
- Built-in statistical functions for differential expression analysis and pre-processing.
- Built-in CSV reader and manipulation.
- Graphs for fastq quality analysis, GC plots, and enrichment analysis using Python.
- Different array computation programs and matrix calculations.
- Built-in data structures for data science and statistics.


# <code>Requirements</code>
BioSyS requires Python on your machine for plotting and for providing some functionalities. First of all, download the Python 3.x version. Then, download the supported version of BioSyS and configure it to your Python interpreter.
### Softwares Requirement:
- [.Net](https://dotnet.microsoft.com/en-us/) >= 5.0.x
- [Python](https://www.python.org/) >= 3.0.x
### Packages Requirement:
- [HtmlAgilityPack](https://www.nuget.org/packages/HtmlAgilityPack) >= 1.11.x
- [IronPython](https://www.nuget.org/packages/IronPython) >= 3.4.1
- [MathNet.Numerics](https://www.nuget.org/packages/MathNet.Numerics/6.0.0-beta1) >= 4.9.0
- [Matplotlib](https://pypi.org/project/matplotlib/) >= 3.6.0
- [beautifulsoup4](https://pypi.org/project/beautifulsoup4/) == 4.11.0
- [seaborn](https://pypi.org/project/seaborn/) >= 0.12.1

# <code>Download</code>
You can easly download the package with nuget package manager:</br>
```
PM> NuGet\Install-Package BioSySNet -Version 0.0.1
```
OR</br>

Can also be installed by CLI using dotnet:</br>
```
> dotnet add package BioSySNet --version 0.0.1
```

# ` Uses `

COMING SOON....

<code style="background: red">Sometimes the features of Python functions can be defective and cause errors!</code>