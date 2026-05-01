[![Latest release](https://img.shields.io/github/release/524D/compareMS2.svg)](https://github.com/524D/compareMS2/releases/latest) [![GitHub](https://img.shields.io/badge/github-repo-000.svg?logo=github&labelColor=gray&color=blue)](https://github.com/524D/compareMS2)
[![GitHub](https://img.shields.io/github/license/524D/compareMS2)](https://github.com/524D/compareMS2/blob/main/LICENSE.txt) [![RSD](https://img.shields.io/badge/rsd-compareMS2-00a3e3.svg)](https://research-software.nl/software/comparems2) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14990300.svg)](https://doi.org/10.5281/zenodo.14990300) [![OpenSSF Best Practices](https://www.bestpractices.dev/projects/10480/badge)](https://www.bestpractices.dev/projects/10480) [![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F-green)](https://fair-software.eu) [![bio.tools](https://img.shields.io/badge/bio.tools-compareMS2-005472)](https://bio.tools/compareMS2)


# compareMS2

compareMS2 calculates the global similarity between tandem mass spectrometry datasets.

[1. Introduction](#1-Introduction)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[1.1 What is compareMS2?](#11-What-is-compareMS2)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[1.2 How does compareMS2 differ from other tools?](#12-How-does-compareMS2-differ-from-other-tools)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[1.3 What can compareMS2 be used for?](#13-What-can-compareMS2-be-used-for)  
[2. Installing compareMS2](#2-Installing-compareMS2)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[2.1 Running compareMS2 in development mode](#21-Running-compareMS2-in-development-mode)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[2.2 Building compareMS2](#22-Building-compareMS2)  
[3. Using compareMS2](#3-Using-compareMS2)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[3.1 Configuring compareMS2](#31-Configuring-compareMS2)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[3.2 Calculating distance matrices](#32-Calculating-distance-matrices)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[3.3 Running compareMS2](#33-Running-compareMS2)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[3.4 Molecular phylogenetics](#34-Molecular-phylogenetics)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[3.5 Data quality control](#35-Data-quality-control)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[3.6 Experimental features](#36-Experimental-features)  
[4. compareMS2 tutorial](#4-compareMS2-tutorial)  
[5. Acknowledgements ](#5-Acknowledgements)  
[6. Further reading](#6-Further-reading)  

## 1. Introduction

### 1.1 What is compareMS2?

compareMS2 is a tool for direct comparison of tandem mass spectrometry datasets, typically from liquid chromatography-tandem mass spectrometry (LC-MS/MS), defining similarity as a function of shared (similar) spectra and distance as the inverse of this similarity. Data with identical spectral content thus have similarity 1 and distance 0. The similarity of datasets with no similar spectra tend to 0 (distance +∞) as the size of the sets go to infinity. The extremes of none or all spectra being similar between two LC-MS/MS datasets are unlikely to occur in reality.

### 1.2 How does compareMS2 differ from other tools?

Though compareMS2 is not limited to tandem mass spectra of peptides, it has seen most application to this type of data. There are four broad categories of tools for the analysis of peptide tandem mass spectra in mass spectrometry-based proteomics based on what prior information they utilize. compareMS2 belongs to a class of tools that do not use existing sequence data or libraries of spectra assigned to a specific peptide sequence, but compare tandem mass spectra directly with other tandem mass spectra:

<table>
    <tbody>
        <tr>
            <th rowspan=2 colspan=2></th>
            <th colspan=2>(translated) genome sequences available</th>
        </tr>
        <tr>
            <td>+</td>
            <td>-</td>
        </tr>
        <tr>
            <th rowspan=2>prior/other tandem <br>mass spectra available</th>
            <td>+</td>
            <td>spectral libraries (<a href="https://doi.org/10.1021/ac060279n">BiblioSpec</a>,
                <a href="https://doi.org/10.1002/pmic.200600625">SpectraST</a>, ...)</td>
            <td>direct comparison (<a href="https://doi.org/10.1002/rcm.6162"><b>compareMS2</b></a>,
                <a href="https://doi.org/10.1186/s12859-017-1514-2">DISMS2</a>, ...)</td>
        </tr>
        <tr>
            <td>-</td>
            <td>database search (<a href="https://doi.org/10.1002/(SICI)1522-2683(19991201)20:18<3551::AID-ELPS3551>3.0.CO;2-2">Mascot</a>,
                <a href="https://doi.org/10.1038/msb4100024">Comet</a>, ...)</td>
            <td><i>de novo</i> sequencing (<a href="https://doi.org/10.1002/pro.5560010902">LUTEFISK</a>,
                <a href="https://doi.org/10.1021/ac048788h">PepNovo</a>, ...)</td>
        </tr>
    </tbody>
</table>

### 1.3 What can compareMS2 be used for?

compareMS2 (and similar tools) have extremely broad utility, but have so far seen most utility in data quality control, food/feed species identification and molecular phylogenetics. Molecular phylogenetics is the study of evolution and relatedness of organisms, genes or proteins. The field dates back to [1960](https://doi.org/10.1073/pnas.46.10.1349) using patterns of tryptic peptides separated by paper chromatography. compareMS2 is a [21st-century analogue](https://doi.org/10.1021/acs.jproteome.1c00528), comparing patterns of tryptic peptides as analyzed by tandem mass spectrometry, with the difference that it can use thousands of peptides and that the tandem mass spectra are highly peptide-specific.

However, not only the amino acid sequences of the peptides affect the distance metric in compareMS2, but also the abundance (or coverage) of the proteins. compareMS2 can also be used to quantify the similarity of proteomes from different cell lines or tissues from the same species, before and independently of any protein identification by database or spectral library search.

## 2. Installing compareMS2

The compareMS2 software can be run under Windows (64 bit AMD/Intel), Linux (64 bit AMD/Intel) and MacOS (ARM and Intel).

On Windows and Ubuntu, the easiest way to install compareMS2 is through the [installer](https://github.com/524D/compareMS2/releases/latest/) (under "assets").

Alternatively, and for other platforms, follow the instructions below.

A recent version of `NodeJs` and `yarn` are required. The versions supplied with most Linux distributions are outdated.
We recommend installing the latest versions from the NodeJs and yarn websites. NodeJs must be installed before yarn.

* For NodeJs, use the [installer](https://nodejs.org/en/download), except for Ubuntu/Debian (because no installer is provided).
  For Ubuntu Linux we used the instructions found
  [here](https://www.digitalocean.com/community/tutorials/how-to-install-node-js-on-ubuntu-22-04#option-2-installing-node-js-with-apt-using-a-nodesource-ppa)
* For yarn, we used the instructions found [yarn website](https://classic.yarnpkg.com/lang/en/docs/install/)
  under "Alternatives", installing the "Classic stable" version.

Then run the following on the command line:

```text
git clone https://github.com/524D/compareMS2
cd compareMS2
yarn
```

Ignore npm vulnerability warnings (**don't** run ```npm audit fix```). Since compareMS is not a web app, they are of limited relevance,
and they can't easily be fixed.

### 2.1 Running compareMS2 in development mode

To run compareMS2 in "development mode", simply issue:

```text
yarn start
```

For debug mode (enabling Chrome development tools):

```text
CPM_MS2_DEBUG="x" yarn start
```

### 2.2 Building compareMS2

To build a distributable package (for the platform on which this command is executed):

```text
yarn make
```

For example, the resulting Windows installer can than be found (relative to the compareMS2 main directory) in
`out\make\squirrel.windows\x64\`.

## 3. Using compareMS2

compareMS2 can be used both from the command-line interface (CLI) and through the compareMS2 GUI. Every compareMS2 analysis consists of two phases: (1) pairwise comparison of all LC-MS/MS datasets and (2) calculating a distance matrix from all pairwise comparisons. The compareMS2 GUI provides real-time feedback by continuously updating the distance matrix, and drawing a UPGMA tree at the completion of each row in the (lower triangular) distance matrix. The default distance metric D is symmetric, i.e. the distance from dataset A to dataset B is identical to the distance from dataset B to dataset A. If the distance D(A, B) has already been calculated, there is no need to calculate D(B, A). As every dataset is identical to itself, there is no point in calculating D(A, A) or D(B, B), as these are always zero.

![compareMS2 on primate datasets](./pictures/primates_circular.png)  
Figure 1. Phylogenetic tree based on sample primate [sera datasets](https://osf.io/sg796/) of 1,000 tandem mass spectra, as displayed during a compareMS2 run. This is a good test dataset for compareMS2.

See PRIDE Project [PXD034932](https://www.ebi.ac.uk/pride/archive/projects/PXD034932) for additional compareMS2 test data.

### 3.1 Configuring compareMS2  

The compareMS2 CLI has a small number of parameters, which are:

-A *first dataset filename*  
-B *second dataset filename*
-W *first scan number*, *last scan number*  
-R *first retention time*, *last retention time*  
-c *cutoff for spectral similarity*  
-o *output filename*  
-m *minimum base peak intensity*, *minimum total MS/MS intensity*  
-w *maximum scan number difference*  
-r *maximum retention time difference*  
-p *maximum difference in precursor mass*  
-e *maximum precursor mass measurement error*  
-s *intensity scaling before dot product*  
-n *noise threshold for dot product*  
-d *version of set distance metric*  
-q *version of QC metric*  
-N *include only N most intense spectra in comparison*  
-b *bin size for dot product*  
-I *minimum number of peaks for dot product*  
-L *lower m/z for dot product*  
-U *upper m/z for dot product*  
-x *experimental features*  

The compareMS2 GUI exposes some of these, and determine others automatically, e.g. the dataset filenames from a specified directory.

### 3.2 Calculating distance matrices

Distance matrices are calculated using a separate executable, compareMS2_to_distance_matrices. This can also average the distances for multiple replicates per species for more accurate molecular phylogenetic analysis. For this, a tab-delimited file with filenames and species names are required. If no such file is provided, one is created automatically, using the filenames as sample "species". The distance matrix can currently be saved in the MEGA or Nexus formats. [MEGA](https://www.megasoftware.net/) is recommended for creating trees from compareMS2 results.

### 3.3 Running compareMS2

After specifying the parameters, click on the "Start" button to run compareMS2 on all files in the specified directory. Alternatively, compareMS2 can be run on two specific files using the CLI version.

### 3.4 Molecular phylogenetics

We recommend [MEGA](https://www.megasoftware.net/) creating phylogenetic trees from compareMS2 results. However, most phylogenetic software can take distance matrices as input for UPGMA analysis. This was the original use for which compareMS2 was developed, see the [2012 paper](https://doi.org/10.1002/rcm.6162).

### 3.5 Data quality control

compareMS2 provides a very quick overview of large number of datasets to see if they cluster as expected or if there are outliers. Data of lower quality can thus be detected *before* running them through a data analysis pipeline and statistical analysis. It is not absolutely necessary to include all spectra in the analysis - major discrepancies should be detectable with ~1,000 spectra, if selected systematically. Similarly, compareMS2 can be used to determine the relative importance of factors in sample preparation and analysis, as shown in a [2016 paper](https://doi.org/10.1002/rcm.7494).

In addition, compareMS2 collects metadata on each dataset (by default the number of tandem mass spectra) and visualizes this on top of the hierarcical clustering or phylogenetic tree.

### 3.6 Experimental features

Starting in version 2.0, we have begun to include experimental features in compareMS2. These are only available on the command line, but allow extraction of additional information from the comparisons, such as the distribution of similarity between tandem mass spectra as function of precursor mass measurement error, allowing identification of isotope errors and charge state distributions *before* any database search:

![Experimental feature](./pictures/experimental_features.png)  
Figure 2. Similarity (spectral angle from 0 to 1) of tandem mass spectra plotted against precursor *m*/*z* difference, revealing isotope errors up to at least 2 (corresponding to bands at *m*/*z* difference 2/3 and 2/5) and charge states up to 6 (corresponding to the band at *m*/*z* difference 1/6).

## 4. compareMS2 tutorial

A compareMS2 workshop was held at the [EuBIC Winter School in Winterberg, Germany, January 15-19 2024](https://eubic-ms.org/events/2024-winter-school/). Some of the slides from this workshops can be found [here](https://osf.io/8qjwz), along with a tutorial [here](https://osf.io/e5j7q). All data for the tutorial are also available on [OSF](https://osf.io/sjtrm/).

## 5. Acknowledgements

The developers wish to thank Dr. Michael Dondrup at the University of Bergen for providing changes and additions to make compareMS2 work under macOS. All users and beta testers are also acknowledged for their valuable feedback that helped to improve compareMS2.

## 6. Further reading

compareMS2 and related applications have been described or used in a number of papers:

Shotgun Proteomics Protocol for Insects, Varunjikar MS, Belghit I, Oveland E, Palmblad M and Rasinger JD, *Methods Mol. Biol.* **2884**:81-98, 2025, [doi.org/10.1007/978-1-0716-4298-6_7](https://doi.org/10.1007/978-1-0716-4298-6_7)

Fish species authentication in commercial fish products using mass spectrometry and spectral library matching approach, Varunjikar MS, Pineda-Pampliega J, Belghit I, Palmblad M, Einar Grøsvik B, Meier S, Asgeir Olsvik P, Lie KK and Rasinger JD, *Food Res. Int.* **192**:114785, 2024, [doi.org/10.1016/j.foodres.2024.114785](https://doi.org/10.1016/j.foodres.2024.114785)

compareMS2 2.0: An Improved Software for Comparing Tandem Mass Spectrometry Datasets, Marissen M, Varunjikar MS, Laros JFJ, Rasinger JD, Neely BA and Palmblad M, *J. Proteome Res.*  **22(2)**:514–519, 2023, [doi.org/10.1021/acs.jproteome.2c00457](https://doi.org/10.1021/acs.jproteome.2c00457)

Shotgun proteomics approaches for authentication, biological analyses, and allergen detection in feed and food-grade insect species, Varunjikar MS, Belghit I, Gjerde J, Palmblad M, Oveland E and Rasinger JD, *Food Control* **131**, 2022, [doi.org/10.1016/j.foodcont.2022.108888](https://doi.org/10.1016/j.foodcont.2022.108888)

Comparing novel shotgun DNA sequencing and state-of-the-art proteomics approaches for authentication of fish species in mixed samples, Varunjikar MS, Moreno-Ibarguen C, Andrade-Martinez JS, Tung HS, Belghit I, Palmblad M, Olsvik PA, Reyes A, Rasinger JD and Lie KK, *Food Control* **131**:108417, 2022, [doi.org/10.1016/j.foodcont.2021.108417](https://doi.org/10.1016/j.foodcont.2021.108417)

Rewinding the molecular clock: looking at pioneering molecular phylogenetics experiments in the light of proteomics, Neely B and Palmblad M, *J. Proteome Res.* **20(10)**:4640-4645, 2021, [doi.org/10.1021/acs.jproteome.1c00528](https://doi.org/10.1021/acs.jproteome.1c00528)

Future feed control – Tracing banned bovine material in insect meal. Belghit I, Varunjikar M, Lecrenier MC, Steinhilber A, Niedzwiecka A, Wang YV, Dieu M, Azzollini D, Lie K, Lock EJ, Berntssen MHG, Renard P, Zagon J, Fumière O, van Loon JJA, Larsen T, Poetz O, Braeuning A, Palmblad M and Rasinger JD, *Food Control* **128**:108183, 2021, [doi.org/10.1016/j.foodcont.2021.108183](https://doi.org/10.1016/j.foodcont.2021.108183)

Species-Specific Discrimination of Insect Meals for Aquafeeds by Direct Comparison of Tandem Mass Spectra. Belghit I, Lock EJ, Fumière O, Lecrenier MC, Renard P, Dieu M, Berntssen MHG, Palmblad M and Rasinger JD, *Animals* **9(5)**:222, 2019 [doi.org/10.3390/ani9050222](https://doi.org/10.3390/ani9050222)

Palaeoproteomics of bird bones for taxonomic classification. Horn IR, Kenens Y, Palmblad M, van der Plas-Duivesteijn SJ, Langeveld BW, Meijer HJM, Dalebout H, Marissen RJ, Fischer A, Vincent Florens FB, Niemann J, Rijsdijk KF, Schulp AS, Laros JFJ and Gravendeel B, *Zoological Journal of the Linnean Society* **186(3)**:650–665, 2019, [doi.org/10.1093/zoolinnean/zlz012](https://doi.org/10.1093/zoolinnean/zlz012)

Species and tissues specific differentiation of processed animal proteins in aquafeeds using proteomics tools. Rasinger JD, Marbaix H, Dieu M, Fumière O, Mauro S, Palmblad M, Raes M and Berntssen MHG, *J. Proteomics* **147**:125-131, 2016, [doi.org/10.1016/j.jprot.2016.05.036](https://doi.org/10.1016/j.jprot.2016.05.036)   

Authentication of closely related fish and derived fish products using tandem mass spectrometry and spectral library matching. Nessen M, van der Zwaan D, Greevers S, Dalebout H, Staats M, Kok E and Palmblad M, *J. Agric. Food Chem.* **64(18)**:3669-3677, 2016, [doi.org/10.1021/acs.jafc.5b05322](https://doi.org/10.1021/acs.jafc.5b05322)

Identification of meat products by shotgun spectral matching. Ohana D, Dalebout H, Marissen RJ, Wulff J, Bergquist J, Deelder AM and Palmblad M, *Food Chem.* **203**:28-34, 2016, [doi.org/10.1016/j.foodchem.2016.01.138](https://doi.org/10.1016/j.foodchem.2016.01.138)   

Differentiating samples and experimental protocols by direct comparison of tandem mass spectra. van der Plas-Duivesteijn SJ, Wulff T, Klychnikov O, Ohana D, Dalebout H, van Veelen PA, de Keijzer J, Nessen MA, van der Burgt YEM, Deelder AM and Palmblad M, *Rapid Commun. Mass Spectrom.* **30**:731-738, 2016, [doi.org/10.1002/rcm.7494](https://doi.org/10.1002/rcm.7494) 

Identifying Proteins in Zebrafish Embryos Using Spectral Libraries Generated from Dissected Adult Organs and Tissues. van der Plas-Duivesteijn SJ, Mohammed Y, Dalebout H, Meijer A, Botermans A, Hoogendijk JL, Henneman AA, Deelder AM, Spaink HP and Palmblad M, *J. Proteome Res.* **13(3)**:1537-1544, 2014, [doi.org/10.1021/pr4010585](https://doi.org/10.1021/pr4010585)

Authentication of Fish Products by Large-Scale Comparison of Tandem Mass Spectra. Wulff  T, Nielsen ME, Deelder AM, Jessen F and Palmblad M, *J. Proteome Res.* **12(11)**:5253-5259, 2013, [doi.org/10.1021/pr4006525](https://doi.org/10.1021/pr4006525)

Molecular phylogenetics by direct comparison of tandem mass spectra. Palmblad M and Deelder AM, *Rapid Commun. Mass Spectrom.* **26(7)**:728-732, 2012, [doi.org/10.1002/rcm.6162](https://doi.org/10.1002/rcm.6162)
