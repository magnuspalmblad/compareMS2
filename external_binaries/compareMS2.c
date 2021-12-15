/*                                                                                                                   */
/* compareMS2 - comparison of tandem mass spectra datasets, e.g. from LC-MS/MS                                       */
/*                                                                                                                   */
/* MIT License                                                                                                       */
/*                                                                                                                   */
/* Copyright (c) 2021 Magnus Palmblad                                                                                */
/*                                                                                                                   */
/* Permission is hereby granted, free of charge, to any person obtaining a copy                                      */
/* of this software and associated documentation files (the "Software"), to deal                                     */
/* in the Software without restriction, including without limitation the rights                                      */
/* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell                                         */
/* copies of the Software, and to permit persons to whom the Software is                                             */
/* furnished to do so, subject to the following conditions:                                                          */
/*                                                                                                                   */
/* The above copyright notice and this permission notice shall be included in all                                    */
/* copies or substantial portions of the Software.                                                                   */
/*                                                                                                                   */
/* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR                                        */
/* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,                                          */
/* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE                                       */
/* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER                                            */
/* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,                                     */
/* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE                                     */
/* SOFTWARE.                                                                                                         */
/*                                                                                                                   */
/* Contact information: n.m.palmblad@lumc.nl                                                                         */
/*                                                                                                                   */
/* Compile with gcc -o compareMS2 compareMS2.c or x86_64-w64-mingw32-gcc compareMS2.c -o compareMS2                  */
/*                                                                                                                   */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

#define MAX_LEN 8192
#define HISTOGRAM_BINS 200
#define	DEFAULT_MIN_BASEPEAK_INTENSITY 0
#define	DEFAULT_MIN_TOTAL_ION_CURRENT 0
#define	DEFAULT_MAX_SCAN_NUMBER_DIFFERENCE 10000
#define	DEFAULT_MAX_RT_DIFFERENCE 60
#define	DEFAULT_MAX_PRECURSOR_DIFFERENCE 2.05
#define	DEFAULT_START_SCAN 1
#define	DEFAULT_END_SCAN 1000000
#define	DEFAULT_START_RT 0
#define	DEFAULT_END_RT 100000
#define	DEFAULT_CUTOFF 0.8
#define	DEFAULT_SCALING 0.5
#define	DEFAULT_NOISE 0
#define	DEFAULT_METRIC 2
#define	DEFAULT_QC 0
#define	DEFAULT_BIN_SIZE 0.2
#define	DEFAULT_MIN_PEAKS 20
#define	DEFAULT_MAX_PEAKS 5000
#define	DEFAULT_MIN_MZ 109
#define	DEFAULT_MAX_MZ 2000
#define	DEFAULT_N_BINS 9455
#define	DEFAULT_TOP_N -1
#define DEFAULT_EXPERIMENTAL_FEATURES 0
#define MASS_DIFF_HISTOGRAM_BINS 320
#define DEFAULT_SCAN_NUMBERS_COULD_BE_READ 0
#define DEFAULT_RTS_COULD_BE_READ 0

/* atol0 acts the same as atol, but handles a null pointer without crashing */
static int atol0(const char *p) {
	if (p == 0) {
		return 0;
	}
	return atol(p);
}

/* atof0 acts the same as atof, but handles a null pointer without crashing */
static double atof0(const char *p) {
	if (p == 0) {
		return 0;
	}
	return atof(p);
}

/* an implementation of quickselect for finding the intensity of the k:th most intense spectra (for filtering) */
#define SWAP(x,y) {double temp=x; x=y; y=temp;}

long partition(double a[], long left, long right, long pIndex) {
	double pivot = a[pIndex];
	SWAP(a[pIndex], a[right]);
	pIndex = left;
	for (long i = left; i < right; i++)
		if (a[i] >= pivot) {
			SWAP(a[i], a[pIndex]);
			pIndex++;
		}
	SWAP(a[pIndex], a[right]);
	return pIndex;
}

double quickSelect(double A[], long left, long right, long k) {
	if (left == right)
		return A[left];
	long pIndex = left + rand() % (right - left + 1);
	pIndex = partition(A, left, right, pIndex);
	if (k == pIndex)
		return A[k];
	else if (k < pIndex)
		return quickSelect(A, left, pIndex - 1, k);
	else
		return quickSelect(A, pIndex + 1, right, k);
}

/* compareMS2 main function */

int main(int argc, char *argv[]) {
	FILE *datasetA, *datasetB, *output;
	char datasetAFilename[MAX_LEN], datasetBFilename[MAX_LEN],
			outputFilename[MAX_LEN], experimentalOutputFilename[MAX_LEN],
			temp[MAX_LEN], line[MAX_LEN], *p, metric, qc, experimentalFeatures,
			datasetAScanNumbersCouldBeRead, datasetBScanNumbersCouldBeRead,
			datasetARTsCouldBeRead, datasetBRTsCouldBeRead;
	long i, j, k, datasetASize, datasetBSize, startScan, endScan, nComparisons,
			minPeaks, maxPeaks, nBins, nPeaks, topN, histogram[HISTOGRAM_BINS],
			massDiffHistogram[HISTOGRAM_BINS], **massDiffDotProductHistogram,
			greaterThanCutoff, sAB, sBA, datasetAActualCompared,
			datasetBActualCompared;
	double minBasepeakIntensity, minTotalIonCurrent, maxScanNumberDifference,
			maxRTDifference, startRT, endRT, maxPrecursorDifference, cutoff,
			scaling, noise, minMz, maxMz, binSize, *datasetAIntensities,
			*datasetBIntensities, datasetACutoff, datasetBCutoff, dotProd,
			maxDotProd, dotProdSum, squareSum, rootSquareSum;

	typedef struct {
		long scan; /* scan number */
		double rt; /* retention time in seconds */
		double *mz; /* measured m/z */
		double *intensity; /* measured intensities */
		char charge; /* deconvoluted charge (in MGF file) */
		double *bin; /* binned spectra */
		double precursorMz; /* precursor m/z */
		int nPeaks; /* number of peaks in spectrum */
		double basepeakIntensity; /* basepeak intensity */
		double totalIonCurrent; /* total ion current for spectrum */
	} DatasetType;

	DatasetType *A, *B;

	/* parsing command line parameters */

	if ((argc == 2)
			&& ((strcmp(argv[1], "--help") == 0)
					|| (strcmp(argv[1], "-help") == 0)
					|| (strcmp(argv[1], "-h") == 0))) /* want help? */
			{
		printf(
				"compareMS2 - (c) Magnus Palmblad 2010-2021\n\ncompareMS2 is developed to compare, globally, all MS/MS spectra between two datasets in MGF acquired under similar conditions, or aligned so that they are comparable. This may be useful for molecular phylogenetics based on shared peptide sequences quantified by the share of highly similar tandem mass spectra. The similarity between a pair of tandem mass spectra is calculated essentially as in SpectraST [see Lam et al. Proteomics 2007, 7, 655-667 (2007)].\n\nusage: compareMS2 -A <first dataset filename> -B <second dataset filename> [-W <first scan number>,<last scan number> -R <first retention time>,<last retention time> -c <score cutoff, default=0.8> -o <output filename> -m<minimum base peak signal in MS/MS spectrum for comparison>,<minimum total ion signal in MS/MS spectrum for comparison> -w <maximum scan number difference> -r <maximum retention time difference> -p <maximum difference in precursor mass> -e <maximum mass measurement error in MS/MS> -s <scaling power> -n <noise threshold> -d <distance metric (0, 1 or 2)> -q <QC measure (0)>]\n");
		return 0;
	}

	/* test for correct number of parameters */

	if (argc < 3) {
		printf(
				"usage: compareMS2 -A <first dataset filename> -B <second dataset filename> [-W <first scan number>,<last scan number> -R <first retention time>,<last retention time> -c <score cutoff> -o <output filename> -m <minimum base peak signal in MS/MS spectrum for comparison>,<minimum total ion signal in MS/MS spectrum for comparison> -w <maximum scan number difference> -r <maximum retention time difference> -p <maximum difference in precursor mass> -e <maximum mass measurement error> -s <scaling power> -n <noise threshold> -d <distance metric (0, 1 or 2)> -q <QC measure (0)>] (type compareMS2 --help for more information)\n");
		return -1;
	}

	/* assign default values */

	strcpy(outputFilename, "output.txt");
	minBasepeakIntensity = DEFAULT_MIN_BASEPEAK_INTENSITY;
	minTotalIonCurrent = DEFAULT_MIN_TOTAL_ION_CURRENT;
	maxScanNumberDifference = DEFAULT_MAX_SCAN_NUMBER_DIFFERENCE;
	maxRTDifference = DEFAULT_MAX_RT_DIFFERENCE;
	maxPrecursorDifference = DEFAULT_MAX_PRECURSOR_DIFFERENCE;
	startScan = DEFAULT_START_SCAN;
	endScan = DEFAULT_END_SCAN;
	startRT = DEFAULT_START_RT;
	endRT = DEFAULT_END_RT;
	cutoff = DEFAULT_CUTOFF;
	scaling = DEFAULT_SCALING;
	noise = DEFAULT_NOISE;
	metric = DEFAULT_METRIC;
	qc = DEFAULT_QC;
	binSize = DEFAULT_BIN_SIZE;
	minPeaks = DEFAULT_MIN_PEAKS;
	maxPeaks = DEFAULT_MAX_PEAKS;
	minMz = DEFAULT_MIN_MZ;
	maxMz = DEFAULT_MAX_MZ;
	nBins = DEFAULT_N_BINS;
	topN = DEFAULT_TOP_N;
	experimentalFeatures = DEFAULT_EXPERIMENTAL_FEATURES;
	datasetAScanNumbersCouldBeRead = DEFAULT_SCAN_NUMBERS_COULD_BE_READ;
	datasetBScanNumbersCouldBeRead = DEFAULT_SCAN_NUMBERS_COULD_BE_READ;
	datasetARTsCouldBeRead = DEFAULT_RTS_COULD_BE_READ;
	datasetBRTsCouldBeRead = DEFAULT_RTS_COULD_BE_READ;
	strcpy(experimentalOutputFilename, "experimental_output.txt");

	/* read and replace parameter values */

	for (i = 1; i < argc; i++) {
		if ((argv[i][0] == '-') && (argv[i][1] == 'A')) /* dataset A filename */
			strcpy(datasetAFilename,
					&argv[strlen(argv[i]) > 2 ? i : i + 1][
							strlen(argv[i]) > 2 ? 2 : 0]);
		if ((argv[i][0] == '-') && (argv[i][1] == 'B')) /* dataset B filename */
			strcpy(datasetBFilename,
					&argv[strlen(argv[i]) > 2 ? i : i + 1][
							strlen(argv[i]) > 2 ? 2 : 0]);
		if ((argv[i][0] == '-') && (argv[i][1] == 'W')) { /* range of spectra (in scans) */
			strcpy(temp,
					&argv[strlen(argv[i]) > 2 ? i : i + 1][
							strlen(argv[i]) > 2 ? 2 : 0]);
			p = strtok(temp, ",");
			startScan = atol0(p);
			p = strtok('\0', ",");
			endScan = atol0(p);
		}
		if ((argv[i][0] == '-') && (argv[i][1] == 'R')) { /* range of spectra (in RT) */
			strcpy(temp,
					&argv[strlen(argv[i]) > 2 ? i : i + 1][
							strlen(argv[i]) > 2 ? 2 : 0]);
			p = strtok(temp, ",");
			startRT = atof0(p);
			p = strtok('\0', ",");
			endRT = atof0(p);
		}
		if ((argv[i][0] == '-') && (argv[i][1] == 'o')) /* output filename */
			strcpy(outputFilename,
					&argv[strlen(argv[i]) > 2 ? i : i + 1][
							strlen(argv[i]) > 2 ? 2 : 0]);
		if ((argv[i][0] == '-') && (argv[i][1] == 'p')) /* maximum precursor m/z difference */
			maxPrecursorDifference = atof(
					&argv[strlen(argv[i]) > 2 ? i : i + 1][
							strlen(argv[i]) > 2 ? 2 : 0]);
		if ((argv[i][0] == '-') && (argv[i][1] == 'm')) { /* minimum basepeak and total intensity */
			strcpy(temp,
					&argv[strlen(argv[i]) > 2 ? i : i + 1][
							strlen(argv[i]) > 2 ? 2 : 0]);
			p = strtok(temp, ",");
			minBasepeakIntensity = atof(p);
			p = strtok('\0', ",");
			minTotalIonCurrent = atof(p);
		}
		if ((argv[i][0] == '-') && (argv[i][1] == 'w')) /* maximum scan number difference */
			maxScanNumberDifference = atof0(
					&argv[strlen(argv[i]) > 2 ? i : i + 1][
							strlen(argv[i]) > 2 ? 2 : 0]);
		if ((argv[i][0] == '-') && (argv[i][1] == 'r')) /* maximum RT difference */
			maxRTDifference = atof0(
					&argv[strlen(argv[i]) > 2 ? i : i + 1][
							strlen(argv[i]) > 2 ? 2 : 0]);
		if ((argv[i][0] == '-') && (argv[i][1] == 'c')) /* cutoff for spectral similarity */
			cutoff = atof(
					&argv[strlen(argv[i]) > 2 ? i : i + 1][
							strlen(argv[i]) > 2 ? 2 : 0]);
		if ((argv[i][0] == '-') && (argv[i][1] == 's')) /* intensity scaling for dot product */
			scaling = atof(
					&argv[strlen(argv[i]) > 2 ? i : i + 1][
							strlen(argv[i]) > 2 ? 2 : 0]);
		if ((argv[i][0] == '-') && (argv[i][1] == 'n')) /* noise threshold for dot product */
			noise = atof(
					&argv[strlen(argv[i]) > 2 ? i : i + 1][
							strlen(argv[i]) > 2 ? 2 : 0]);
		if ((argv[i][0] == '-') && (argv[i][1] == 'd')) /* version of set distance metric */
			metric = atoi(
					&argv[strlen(argv[i]) > 2 ? i : i + 1][
							strlen(argv[i]) > 2 ? 2 : 0]);
		if ((argv[i][0] == '-') && (argv[i][1] == 'q')) /* version of QC metric */
			qc = atoi(
					&argv[strlen(argv[i]) > 2 ? i : i + 1][
							strlen(argv[i]) > 2 ? 2 : 0]);
		if ((argv[i][0] == '-') && (argv[i][1] == 'N')) /* compare only the N most intense spectra */
			topN = atoi(
					&argv[strlen(argv[i]) > 2 ? i : i + 1][
							strlen(argv[i]) > 2 ? 2 : 0]);
		if ((argv[i][0] == '-') && (argv[i][1] == 'b')) /* bin size (advanced parameter) */
			binSize = atof(
					&argv[strlen(argv[i]) > 2 ? i : i + 1][
							strlen(argv[i]) > 2 ? 2 : 0]);
		if ((argv[i][0] == '-') && (argv[i][1] == 'I')) /* minimum number of peaks (advanced parameter) */
			minPeaks = atoi(
					&argv[strlen(argv[i]) > 2 ? i : i + 1][
							strlen(argv[i]) > 2 ? 2 : 0]);
		if ((argv[i][0] == '-') && (argv[i][1] == 'L')) /* minimum m/z for dot product (advanced parameter) */
			minMz = atof(
					&argv[strlen(argv[i]) > 2 ? i : i + 1][
							strlen(argv[i]) > 2 ? 2 : 0]);
		if ((argv[i][0] == '-') && (argv[i][1] == 'U')) /* maximum m/z for dot product (advanced parameter) */
			maxMz = atof(
					&argv[strlen(argv[i]) > 2 ? i : i + 1][
							strlen(argv[i]) > 2 ? 2 : 0]);
		if ((argv[i][0] == '-') && (argv[i][1] == 'x')) /* level of experimental features enabled */
			experimentalFeatures = atoi(
					&argv[strlen(argv[i]) > 2 ? i : i + 1][
							strlen(argv[i]) > 2 ? 2 : 0]);
	}

	if (((maxMz - minMz) / binSize) == floor((maxMz - minMz) / binSize))
		nBins = (int) ((maxMz - minMz) / binSize);
	else
		nBins = (int) ((maxMz - minMz) / binSize) + 1;
	printf("spectrum bin size %1.3f Th -> %ld bins in [%.3f,%.3f]\n", binSize,
			nBins, minMz, maxMz);
	fflush(stdout);

	printf(
			"scan range=[%ld,%ld]\nmax scan difference=%.2f\nmax m/z difference=%.4f\nscaling=^%.2f\nnoise=%.1f\nmin basepeak intensity=%.2f\nmin total ion current=%.2f\n",
			startScan, endScan, maxScanNumberDifference, maxPrecursorDifference,
			scaling, noise, minBasepeakIntensity, minTotalIonCurrent);
	fflush(stdout);

	/* check MGF dataset A for number of MS/MS spectra */
	if ((datasetA = fopen(datasetAFilename, "r")) == NULL) {
		printf("error opening dataset A %s for reading", datasetAFilename);
		return -1;
	}
	printf("checking dataset A (\"%s\")...", datasetAFilename);
	fflush(stdout);
	datasetASize = 0;
	nPeaks = 0;
	maxPeaks = 0;
	while (fgets(line, MAX_LEN, datasetA) != NULL) {
		if (strcmp(line, "\n") == 0)
			continue;
		p = strtok(line, "=");
		if (strcmp("TITLE", p) == 0) {
			datasetASize++;
			if (nPeaks > maxPeaks)
				maxPeaks = nPeaks;
			nPeaks = 0;
		}
		if (isdigit(p[0]))
			nPeaks++;
	}
	printf("done (contains %ld MS2 spectra)\n", datasetASize);
	fflush(stdout);
	fclose(datasetA);

	if ((topN > -1) && (topN < datasetASize)) {
		printf("filtering top-%ld spectra...", topN);
		fflush(stdout);
		datasetAIntensities = malloc(datasetASize * sizeof(double));
		if ((datasetA = fopen(datasetAFilename, "r")) == NULL) {
			printf("error opening dataset A %s for reading", datasetAFilename);
			return -1;
		}
		i = -1;
		while (fgets(line, MAX_LEN, datasetA) != NULL) {
			if (strcmp(line, "\n") == 0)
				continue;
			p = strtok(line, " \t");
			if (strcmp("BEGIN", p) == 0) {
				i++;
				datasetAIntensities[i] = 0.0;
			}
			if (isdigit(p[0])) {
				p = strtok('\0', " \t");
				datasetAIntensities[i] += atof0(p);
			}
		}

		datasetACutoff = quickSelect(datasetAIntensities, 0, i, topN); /* quickselect top-Nth intensity */
		printf("done (ion intensity threshold %.3f)\n", datasetACutoff);
		fflush(stdout);
		fclose(datasetA);
	}

	if ((datasetB = fopen(datasetBFilename, "r")) == NULL) {
		printf("error opening dataset B %s for reading", datasetBFilename);
		return -1;
	}
	printf("checking dataset B (\"%s\")...", datasetBFilename);
	fflush(stdout);
	datasetBSize = 0;
	nPeaks = 0;
	while (fgets(line, MAX_LEN, datasetB) != NULL) {
		if (strcmp(line, "\n") == 0)
			continue;
		p = strtok(line, "=");
		if (strcmp("TITLE", p) == 0) {
			datasetBSize++;
			if (nPeaks > maxPeaks)
				maxPeaks = nPeaks;
			nPeaks = 0;
		}
	}
	printf("done (contains %ld MS2 spectra)\n", datasetBSize);
	fflush(stdout);
	fclose(datasetB);

	if ((topN > -1) && (topN < datasetBSize)) {
		printf("filtering top-%ld spectra...", topN);
		fflush(stdout);
		datasetBIntensities = malloc(datasetBSize * sizeof(double));
		if ((datasetB = fopen(datasetBFilename, "r")) == NULL) {
			printf("error opening dataset B %s for reading", datasetBFilename);
			return -1;
		}
		i = -1;
		while (fgets(line, MAX_LEN, datasetB) != NULL) {
			if (strcmp(line, "\n") == 0)
				continue;
			p = strtok(line, " \t");
			if (strcmp("BEGIN", p) == 0) {
				i++;
				datasetBIntensities[i] = 0.0;
			}
			if (isdigit(p[0])) {
				p = strtok('\0', " \t");
				datasetBIntensities[i] += atof0(p);
			}
		}

		datasetBCutoff = quickSelect(datasetBIntensities, 0, i, topN); /* quickselect top-Nth intensity */
		printf("done (ion intensity threshold %.3f)\n", datasetBCutoff);
		fflush(stdout);
		fclose(datasetB);
	}

	/* allocate memory */

	printf("allocating memory...");
	fflush(stdout);
	A = (DatasetType*) malloc(datasetASize * sizeof(DatasetType));
	B = (DatasetType*) malloc(datasetBSize * sizeof(DatasetType));
	if (experimentalFeatures == 1) {
		massDiffDotProductHistogram = (long**) malloc(
		MASS_DIFF_HISTOGRAM_BINS * sizeof(long*));
		for (i = 0; i < MASS_DIFF_HISTOGRAM_BINS; i++)
			massDiffDotProductHistogram[i] = malloc(
			HISTOGRAM_BINS * sizeof(long));
	}

	/* read in tandem mass spectra from MGF files */

	printf("done\nreading %ld MS2 spectra from %s...", datasetASize,
			datasetAFilename);
	fflush(stdout);
	if ((datasetA = fopen(datasetAFilename, "r")) == NULL) {
		printf("error opening dataset A %s for reading", datasetAFilename);
		return -1;
	}
	i = 0;
	j = 0;
	while (fgets(line, MAX_LEN, datasetA) != NULL) {
		if (strcmp(line, "\n") == 0)
			continue;
		p = strtok(line, " \t");
		if (strspn("PEPMASS", p) > 6) {
			A[i].precursorMz = atof(strpbrk(p, "0123456789"));
			A[i].mz = (double*) malloc(maxPeaks * sizeof(double));
			A[i].intensity = (double*) malloc(maxPeaks * sizeof(double));
			A[i].bin = (double*) malloc(nBins * sizeof(double));
		}
		if (strspn("CHARGE", p) > 5)
			A[i].charge = (char) atoi(strpbrk(p, "0123456789"));
		if (isdigit(p[0])) {
			A[i].mz[j] = atof(p);
			p = strtok('\0', " \t");
			if (j < maxPeaks) {
				A[i].intensity[j] = atof(p);
				j++;
			}
		}
		A[i].scan = startScan; /* default if no scan information is available */
		if (strspn("SCANS", p) > 4) { /* MGFs with SCANS attributes */
			A[i].scan = (long) atol0(strpbrk(p, "0123456789"));
			// printf("A[%ld].scan = %ld\n", i, A[i].scan); fflush(stdout);
			datasetAScanNumbersCouldBeRead = 1;
			continue;
		}
		if (strncmp("###MSMS:", p, 8) == 0) { /* Bruker-style MGFs */
			p = strtok('\0', " \t");
			A[i].scan = (long) atol0(strpbrk(p, "0123456789"));
			// printf("A[%ld].scan = %ld\n", i, A[i].scan); fflush(stdout);
			datasetAScanNumbersCouldBeRead = 1;
			continue;
		}
		if (strspn("TITLE", p) > 4) { /* msconvert-style MGFs with NativeID and scan= */
			while (p != NULL) {
				if (strstr(p, "scan=") != NULL) {
					A[i].scan = (long) atol0(strpbrk(p, "0123456789"));
					// printf("A[%ld].scan = %ld\n", i, A[i].scan); fflush(stdout);
					datasetAScanNumbersCouldBeRead = 1;
				}
				p = strtok('\0', " \t");
			}
			continue;
		}
		A[i].rt = startRT; /* default if no scan information is available */
			if (strspn("RTINSECONDS", p) > 10) { /* MGFs with RTINSECONDS attributes */
				A[i].rt = (double) atof0(strpbrk(p, "0123456789"));
				// printf("A[%ld].rt = %ld\n", i, A[i].scan); fflush(stdout);
				datasetARTsCouldBeRead = 1;
				continue;
			}
		if (strcmp("END", p) == 0) {
			A[i].nPeaks = j;
			A[i].basepeakIntensity = 0;
			A[i].totalIonCurrent = 0;
			for (k = 1; k <= j; k++) {
				if (A[i].basepeakIntensity < A[i].intensity[k])
					A[i].basepeakIntensity = A[i].intensity[k];
				A[i].totalIonCurrent = A[i].totalIonCurrent + A[i].intensity[k];
			}
			i++;
			j = 0;
		}
	}

	printf("done\nreading %ld MS2 spectra from %s...", datasetBSize,
			datasetBFilename);
	fflush(stdout);
	if ((datasetB = fopen(datasetBFilename, "r")) == NULL) {
		printf("error opening dataset B %s for reading", datasetBFilename);
		return -1;
	}
	i = 0;
	j = 0;
	while (fgets(line, MAX_LEN, datasetB) != NULL) {
		if (strcmp(line, "\n") == 0)
			continue;
		p = strtok(line, " \t");
		if (strspn("PEPMASS", p) > 6) {
			B[i].precursorMz = atof(strpbrk(p, "0123456789"));
			B[i].mz = (double*) malloc(maxPeaks * sizeof(double));
			B[i].intensity = (double*) malloc(maxPeaks * sizeof(double));
			B[i].bin = (double*) malloc(nBins * sizeof(double));
		}
		if (strspn("CHARGE", p) > 5)
			B[i].charge = (char) atoi(strpbrk(p, "0123456789"));
		if (isdigit(p[0])) {
			B[i].mz[j] = atof(p);
			p = strtok('\0', " \t");
			if (j < maxPeaks) {
				B[i].intensity[j] = atof(p);
				j++;
			}
		}
		B[i].scan = startScan; /* default if no scan information is available */
		if (strspn("SCANS", p) > 4) { /* MGFs with SCANS attributes */
			B[i].scan = (long) atol0(strpbrk(p, "0123456789"));
			// printf("B[%ld].scan = %ld\n", i, B[i].scan); fflush(stdout);
			datasetBScanNumbersCouldBeRead = 1;
			continue;
		}
		if (strncmp("###MSMS:", p, 8) == 0) { /* Bruker-style MGFs */
			p = strtok('\0', " \t");
			B[i].scan = (long) atol0(strpbrk(p, "0123456789"));
			// printf("B[%ld].scan = %ld\n", i, B[i].scan); fflush(stdout);
			datasetBScanNumbersCouldBeRead = 1;
			continue;
		}
		if (strspn("TITLE", p) > 4) { /* msconvert-style MGFs with NativeID and scan= */
			while (p != NULL) {
				if (strstr(p, "scan=") != NULL) {
					B[i].scan = (long) atol0(strpbrk(p, "0123456789"));
					// printf("B[%ld].scan = %ld\n", i, B[i].scan); fflush(stdout);
					datasetBScanNumbersCouldBeRead = 1;
				}
				p = strtok('\0', " \t");
			}
			continue;
		}
		B[i].rt = startRT; /* default if no scan information is available */
			if (strspn("RTINSECONDS", p) > 10) { /* MGFs with RTINSECONDS attributes */
				B[i].rt = (double) atof0(strpbrk(p, "0123456789"));
				// printf("B[%ld].rt = %ld\n", i, B[i].scan); fflush(stdout);
				datasetBRTsCouldBeRead = 1;
				continue;
			}
		if (strcmp("END", p) == 0) {
			B[i].nPeaks = j;
			B[i].basepeakIntensity = 0;
			B[i].totalIonCurrent = 0;
			for (k = 1; k <= j; k++) {
				if (B[i].basepeakIntensity < B[i].intensity[k])
					B[i].basepeakIntensity = B[i].intensity[k];
				B[i].totalIonCurrent = B[i].totalIonCurrent + B[i].intensity[k];
			}
			i++;
			j = 0;
		}
	}
	printf("done\n");

	if (datasetAScanNumbersCouldBeRead == 0) {
		printf("warning: scan numbers could not be read from dataset A (%s)\n",
				datasetAFilename);
	}
	if (datasetBScanNumbersCouldBeRead == 0) {
		printf("warning: scan numbers could not be read from dataset B (%s)\n",
				datasetBFilename);
	}
	if ((datasetAScanNumbersCouldBeRead == 0)
			|| (datasetBScanNumbersCouldBeRead == 0)) {
		maxScanNumberDifference = DEFAULT_MAX_SCAN_NUMBER_DIFFERENCE;
		printf("scan filters will be ignored\n");
	}
	fflush(stdout);

	if (datasetARTsCouldBeRead == 0) {
		printf("warning: retention times could not be read from dataset A (%s)\n",
				datasetAFilename);
	}
	if (datasetBRTsCouldBeRead == 0) {
		printf("warning: retention times could not be read from dataset B (%s)\n",
				datasetBFilename);
	}
	if ((datasetARTsCouldBeRead == 0)
			|| (datasetBRTsCouldBeRead == 0)) {
		maxRTDifference = DEFAULT_MAX_RT_DIFFERENCE;
		printf("retention time filters will be ignored\n");
	}
	fflush(stdout);

	printf("scaling, normalizing and binning %ld MS2 spectra from %s..",
			datasetASize, datasetAFilename);
	fflush(stdout);
	for (j = 0; j < datasetASize; j++) {
		for (k = 0; k < A[j].nPeaks; k++)
			A[j].intensity[k] =
					A[j].intensity[k] > noise ?
							pow(A[j].intensity[k], scaling) : 0; /* nth root scaling and noise removal */
		squareSum = 0;
		for (k = 0; k < A[j].nPeaks; k++)
			squareSum += A[j].intensity[k] * A[j].intensity[k];
		rootSquareSum = sqrt(squareSum); /* normalization against spectral vector magnitude */
		for (k = 0; k < nBins; k++)
			A[j].bin[k] = 0; /* set all bins to zero */
		for (k = 0; k < A[j].nPeaks; k++) {
			A[j].intensity[k] = A[j].intensity[k] / rootSquareSum;
			if ((A[j].mz[k] >= minMz) && (A[j].mz[k] < maxMz))
				A[j].bin[(long) floor(
						binSize * (A[j].mz[k] - minMz) + binSize / 2)] +=
						A[j].intensity[k];
		}
		squareSum = 0;
		for (k = 0; k < nBins; k++)
			squareSum += A[j].bin[k] * A[j].bin[k];
		rootSquareSum = sqrt(squareSum);
		for (k = 0; k < nBins; k++)
			A[j].bin[k] = A[j].bin[k] / rootSquareSum; /* normalize binned spectra to binned vector magnitude */
	}

	printf(".done\nscaling, normalizing and binning %ld MS2 spectra from %s..",
			datasetBSize, datasetBFilename);
	fflush(stdout);
	for (j = 0; j < datasetBSize; j++) {
		for (k = 0; k < B[j].nPeaks; k++)
			B[j].intensity[k] =
					B[j].intensity[k] > noise ?
							pow(B[j].intensity[k], scaling) : 0; /* nth root scaling and noise removal */
		squareSum = 0;
		for (k = 0; k < B[j].nPeaks; k++)
			squareSum += B[j].intensity[k] * B[j].intensity[k];
		rootSquareSum = sqrt(squareSum); /* normalization against spectral vector magnitude */
		for (k = 0; k < nBins; k++)
			B[j].bin[k] = 0; /* set all bins to zero */
		for (k = 0; k < B[j].nPeaks; k++) {
			if ((B[j].mz[k] >= minMz) && (B[j].mz[k] < maxMz))
				B[j].bin[(long) floor(
						binSize * (B[j].mz[k] - minMz) + binSize / 2)] +=
						B[j].intensity[k];
		} /* populate bins */
		squareSum = 0;
		for (k = 0; k < nBins; k++)
			squareSum += B[j].bin[k] * B[j].bin[k];
		rootSquareSum = sqrt(squareSum);
		for (k = 0; k < nBins; k++)
			B[j].bin[k] = B[j].bin[k] / rootSquareSum; /* normalize binned spectra to binned vector magnitude */
	}

	/* go through spectra (entries) in dataset A and compare with those in dataset B and vice versa */

	printf(".done\nmatching spectra and computing set distance.");
	fflush(stdout);
	dotProdSum = 0.0;
	nComparisons = 0;
	greaterThanCutoff = 0;
	sAB = 0;
	sBA = 0;
	datasetAActualCompared = 0;
	datasetBActualCompared = 0;
	for (i = 0; i < HISTOGRAM_BINS; i++) {
		histogram[i] = 0;
		massDiffHistogram[i] = 0;
	}
	if (experimentalFeatures == 1) {
		for (i = 0; i < HISTOGRAM_BINS; i++)
			for (j = 0; j < MASS_DIFF_HISTOGRAM_BINS; j++)
				massDiffDotProductHistogram[j][i] = 0;
	}

	for (i = 0; i < datasetASize; i++) {
		if (topN > -1)
			if (topN < datasetASize)
				if (datasetAIntensities[i] <= datasetACutoff)
					continue; /* compare only top-N spectra (including spectra of same intensity as Nth spectrum) */
		/* if(A[i].scan<startScan) continue; */
		/* if(A[i].scan>endScan) continue; */
		if (A[i].basepeakIntensity < minBasepeakIntensity)
			continue;
		if (A[i].totalIonCurrent < minTotalIonCurrent)
			continue;

		maxDotProd = 0.0;
		datasetAActualCompared++;

		for (j = 0; j < datasetBSize; j++) {
			if (topN > -1)
				if (topN < datasetBSize)
					if (datasetBIntensities[j] <= datasetBCutoff)
						continue;
			if (B[j].basepeakIntensity < minBasepeakIntensity)
				continue;
			if (B[j].totalIonCurrent < minTotalIonCurrent)
				continue;
			if ((A[i].scan - B[j].scan) > maxScanNumberDifference)
				continue;
			if ((B[j].scan - A[i].scan) > maxScanNumberDifference)
				break;
			if ((A[i].rt - B[j].rt) > maxRTDifference)
				continue;
			if ((B[j].rt - A[i].rt) > maxRTDifference)
				break;
			if (fabs(B[j].precursorMz - A[i].precursorMz)
					< maxPrecursorDifference) {
				dotProd = 0;
				for (k = 0; k < nBins; k++)
					dotProd += A[i].bin[k] * B[j].bin[k];
				if (fabs(dotProd) <= 1.00) {
					histogram[(int) (HISTOGRAM_BINS / 2)
							+ (int) floor(dotProd * (HISTOGRAM_BINS / 2 - 1E-9))]++;
					if (experimentalFeatures == 1)
						massDiffDotProductHistogram[(int) (MASS_DIFF_HISTOGRAM_BINS
								/ 2)
								+ (int) floor(
										(B[j].precursorMz - A[i].precursorMz)
												* 99.99999999999999999)][(int) (HISTOGRAM_BINS
								/ 2) /* constant scaling 1 bin = 0.01 m/z units */
						+ (int) floor(dotProd * (HISTOGRAM_BINS / 2 - 1E-9))]++;
				}
				nComparisons++;
				if (dotProd > maxDotProd) {
					maxDotProd = dotProd;
				}
			}
		}
		if (maxDotProd > cutoff)
			greaterThanCutoff++;
		if (maxDotProd > cutoff)
			sAB++;
	}

	printf(".");
	fflush(stdout);

	for (i = 0; i < datasetBSize; i++) {
		if (topN > -1)
			if (topN < datasetBSize)
				if (datasetBIntensities[i] <= datasetBCutoff)
					continue;
		/* if(B[i].scan<startScan) continue; */
		/* if(B[i].scan>endScan) continue; */
		if (B[i].basepeakIntensity < minBasepeakIntensity)
			continue;
		if (B[i].totalIonCurrent < minTotalIonCurrent)
			continue;

		maxDotProd = 0.0;
		datasetBActualCompared++;

		for (j = 0; j < datasetASize; j++) {
			if (topN > -1)
				if (topN < datasetASize)
					if (datasetAIntensities[j] <= datasetACutoff)
						continue;
			if (A[j].basepeakIntensity < minBasepeakIntensity)
				continue;
			if (A[j].totalIonCurrent < minTotalIonCurrent)
				continue;
			if ((B[i].scan - A[j].scan) > maxScanNumberDifference)
				continue;
			if ((A[j].scan - B[i].scan) > maxScanNumberDifference)
				break;
			if ((B[i].rt - A[j].rt) > maxRTDifference)
				continue;
			if ((A[j].rt - B[i].rt) > maxRTDifference)
				break;
			if (fabs(A[j].precursorMz - B[i].precursorMz)
					< maxPrecursorDifference) {
				//printf("2: A[%i] vs B[%i]\n", j, i);
				dotProd = 0;
				for (k = 0; k < nBins; k++)
					dotProd += B[i].bin[k] * A[j].bin[k];
				if (fabs(dotProd) <= 1.00)
					histogram[(int) (HISTOGRAM_BINS / 2)
							+ (int) floor(dotProd * (HISTOGRAM_BINS / 2 - 1E-9))]++;
				nComparisons++;
				if (dotProd > maxDotProd) {
					maxDotProd = dotProd;
				}
			}
		}
		if (maxDotProd > cutoff)
			greaterThanCutoff++; /* counting shared spectra from both datasets */
		if (maxDotProd > cutoff)
			sBA++;
	}

	printf(
			".done (compared %ld (|S_AB|=%ld) spectra from dataset A with %ld (|S_BA|=%ld) spectra from dataset B)\nwriting results to file...",
			datasetAActualCompared, sAB, datasetBActualCompared, sBA);

	/* print output to file */

	if ((output = fopen(outputFilename, "w")) == NULL) {
		printf("error opening output file %s for writing", outputFilename);
		return -1;
	}
	fprintf(output, "dataset_A\t%s\n", datasetAFilename);
	fprintf(output, "dataset_B\t%s\n", datasetBFilename);
	if (metric == 0) /* original metric */
	{
		if (greaterThanCutoff > 0)
			fprintf(output, "set_distance\t%1.10f\n",
					(double) nComparisons / 2.0 / greaterThanCutoff);
		if (greaterThanCutoff == 0)
			fprintf(output, "set_distance\tINF\n");
	}
	if (metric == 1) /* symmetric metric */
	{
		if (greaterThanCutoff > 0)
			fprintf(output, "set_distance\t%1.10f\n",
					(double) (datasetASize + datasetBSize) / greaterThanCutoff);
		if (greaterThanCutoff == 0)
			fprintf(output, "set_distance\tINF\n");
	}
	if (metric == 2) { /* compareMS2 2.0 symmetric metric */
		if ((sAB + sBA) > 0) {
			fprintf(output, "set_distance\t%1.10f\n",
					1.0
							/ ((double) sAB / (2 * datasetASize)
									+ (double) sBA / (2 * datasetBSize)) - 1.0);
		}
		if ((sAB + sBA) == 0) { /* distance between sets with no similar spectra */
			fprintf(output, "set_distance\t%1.10f\n",
					(4.0 * (double) datasetASize * datasetBSize)
							/ (datasetASize + datasetBSize) - 1.0);
		}
	}
	fprintf(output, "set_metric\t%i\n", metric);
	fprintf(output,
			"scan_range\t%ld\t%ld\nmax_scan_diff\t%.5f\nmax_m/z_diff\t%.5f\nscaling_power\t%.5f\nnoise_threshold\t%.5f\nmin_basepeak_intensity\t%.2f\nmin_total_ion_current\t%.2f\n",
			startScan, endScan, maxScanNumberDifference, maxPrecursorDifference,
			scaling, noise, minBasepeakIntensity, minTotalIonCurrent);
	if (qc == 0)
		fprintf(output, "dataset_A_QC\t%.4f\n", (float) datasetASize);
	if (qc == 0)
		fprintf(output, "dataset_B_QC\t%.4f\n", (float) datasetBSize);
	fprintf(output, "n_gt_cutoff\t%ld\n", greaterThanCutoff);
	fprintf(output, "n_comparisons\t%ld\n", nComparisons);
	fprintf(output, "min_peaks\t%ld\n", minPeaks);
	fprintf(output, "max_peaks\t%ld\n", maxPeaks);
	fprintf(output, "m/z_range\t%.4f\t%.4f\n", minMz, maxMz);
	fprintf(output, "m/z_bin_size\t%.4f\n", binSize);
	fprintf(output, "n_m/z_bins\t%ld\n", nBins);
	/* fprintf(output,"histogram (interval, midpoint, comparisons)\n"); */
	for (i = 0; i < HISTOGRAM_BINS; i++)
		fprintf(output, "histogram\t%1.3f\t%1.3f\t%1.3f\t%ld\t%ld\n",
				(double) (i - 100) / 100, (double) (i + 1 - 100) / 100,
				(double) (i + 0.5 - 100) / 100, histogram[i],
				massDiffHistogram[i]);
	fflush(output);
	fclose(output);

	if (experimentalFeatures == 1) {
		if ((output = fopen(experimentalOutputFilename, "w")) == NULL) {
			printf("error opening experimental output file %s for writing",
					experimentalOutputFilename);
			return -1;
		}
		for (i = 0; i < HISTOGRAM_BINS; i++) {
			for (j = 0; j < MASS_DIFF_HISTOGRAM_BINS - 1; j++)
				fprintf(output, "%ld\t", massDiffDotProductHistogram[j][i]);
			fprintf(output, "%ld\n",
					massDiffDotProductHistogram[MASS_DIFF_HISTOGRAM_BINS - 1][i]);
		}
		fflush(output);
		fclose(output);
	}

	/* free memory */
	printf("done\nfreeing memory...");
	fflush(stdout);
	free(A);
	free(B);
	printf("done\n");
	fflush(stdout);

	/* return from main */

	return 0;
}
