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
#define	DEFAULT_MIN_BASEPEAK_INTENSITY 10000
#define	DEFAULT_MAX_SCAN_NUMBER_DIFFERENCE  1500
#define	DEFAULT_MAX_PRECURSOR_DIFFERENCE  2.05
#define	DEFAULT_START_SCAN  1
#define	DEFAULT_END_SCAN  1000000
#define	DEFAULT_CUTOFF  0.8
#define	DEFAULT_SCALING 0.5
#define	DEFAULT_NOISE 10
#define	DEFAULT_METRIC 2
#define	DEFAULT_QC 0
#define	DEFAULT_BIN_SIZE 0.2
#define	DEFAULT_MIN_PEAKS 20
#define	DEFAULT_MAX_PEAKS 5000
#define	DEFAULT_MIN_MZ 109
#define	DEFAULT_MAX_MZ 2000
#define	DEFAULT_N_BINS 9455
#define	DEFAULT_TOP_N -1

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
			outputFilename[MAX_LEN], temp[MAX_LEN], line[MAX_LEN], *p, metric,
			qc;
	long i, j, k, datasetAsize, datasetBsize, startScan, endScan, nComparisons,
			minPeaks, maxPeaks, nBins, nPeaks, topN, histogram[HISTOGRAM_BINS],
			greaterThanCutoff, setSAB, setSBA, datasetAActualCompared,
			datasetBActualCompared;
	double minBasepeakIntensity, minTotalIntensity, maxScanNumberDifference,
			maxPrecursorDifference, cutoff, scaling, noise, minMz, maxMz,
			binSize, *datasetAIntensities, *datasetBIntensities, datasetACutoff,
			datasetBCutoff, dotProd, maxDotProd, dotProdSum, squareSum,
			rootSquareSum;

	typedef struct {
		long scan; /* scan number */
		double *mz; /* measured m/z */
		double *intensity; /* measured intensities */
		char charge; /* deconvoluted charge (in MGF file) */
		double *bin; /* binned spectra */
		double precursorMz; /* precursor m/z */
		int nPeaks; /* number of peaks in spectrum */
		double basePeakIntensity; /* basepeak intensity */
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
				"compareMS2 - (c) Magnus Palmblad 2010-2021\n\ncompareMS2 is developed to compare, globally, all MS/MS spectra between two datasets in MGF acquired under similar conditions, or aligned so that they are comparable. This may be useful for molecular phylogenetics based on shared peptide sequences quantified by the share of highly similar tandem mass spectra. The similarity between a pair of tandem mass spectra is calculated essentially as in SpectraST [see Lam et al. Proteomics 2007, 7, 655-667 (2007)].\n\nusage: compareMS2 -A <first dataset filename> -B <second dataset filename> [-R <first scan number>,<last scan number> -c <score cutoff, default=0.8> -o <output filename> -m<minimum base peak signal in MS/MS spectrum for comparison>,<minimum total ion signal in MS/MS spectrum for comparison> -a <alignment piecewise linear function filename> -w <maximum scan number difference> -p <maximum difference in precursor mass> -e <maximum mass measurement error in MS/MS> -s <scaling power> -n <noise threshold> -d <distance metric (0, 1 or 2)> -q <QC measure (0)>]\n");
		return 0;
	}

	/* test for correct number of parameters */

	if (argc < 3 || argc > 17) {
		printf(
				"usage: compareMS2 -A <first dataset filename> -B <second dataset filename> -R <first scan number>,<last scan number> [-c <score cutoff> -o <output filename> -m <minimum base peak signal in MS/MS spectrum for comparison>,<minimum total ion signal in MS/MS spectrum for comparison> -a <alignment piecewise linear function filename> -w <maximum scan number difference> -p <maximum difference in precursor mass> -e <maximum mass measurement error> -s <scaling power> -n <noise threshold> -d <distance metric (0, 1 or 2)> -q <QC measure (0)>] (type compareMS2 --help for more information)\n");
		return -1;
	}

	/* assign default values */

	strcpy(outputFilename, "output.txt");
	minBasepeakIntensity = DEFAULT_MIN_BASEPEAK_INTENSITY;
	maxScanNumberDifference = DEFAULT_MAX_PRECURSOR_DIFFERENCE;
	maxPrecursorDifference = 2.05;
	startScan = DEFAULT_START_SCAN;
	endScan = DEFAULT_END_SCAN;
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
		if ((argv[i][0] == '-') && (argv[i][1] == 'R')) { /* range of spectra (scans) */
			strcpy(temp,
					&argv[strlen(argv[i]) > 2 ? i : i + 1][
							strlen(argv[i]) > 2 ? 2 : 0]);
			p = strtok(temp, ",");
			startScan = atol0(p);
			p = strtok('\0', ",");
			endScan = atol0(p);
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
			minTotalIntensity = atof(p);
		}
		if ((argv[i][0] == '-') && (argv[i][1] == 'w')) /* maximum scan number difference */
			maxScanNumberDifference = atof(
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
	}

	if (((maxMz - minMz) / binSize) == floor((maxMz - minMz) / binSize))
		nBins = (int) ((maxMz - minMz) / binSize);
	else
		nBins = (int) ((maxMz - minMz) / binSize) + 1;
	printf("spectrum bin size %1.3f Th -> %i bins in [%.3f,%.3f]\n", binSize,
			nBins, minMz, maxMz);
	fflush(stdout);

	printf(
			"scan range=[%i,%i], max scan difference=%.2f, max m/z difference=%.4f, scaling=^%.2f, noise=%.1f\n",
			startScan, endScan, maxScanNumberDifference, maxPrecursorDifference,
			scaling, noise);
	fflush(stdout);

	/* check MGF dataset A for number of MS/MS spectra */
	if ((datasetA = fopen(datasetAFilename, "r")) == NULL) {
		printf("error opening dataset A %s for reading", datasetAFilename);
		return -1;
	}
	printf("checking dataset A (\"%s\")...", datasetAFilename);
	fflush(stdout);
	datasetAsize = 0;
	nPeaks = 0;
	maxPeaks = 0;
	while (fgets(line, 8192, datasetA) != NULL) {
		if (strcmp(line, "\n") == 0)
			continue;
		p = strtok(line, "=");
		if (strcmp("TITLE", p) == 0) {
			datasetAsize++;
			if (nPeaks > maxPeaks)
				maxPeaks = nPeaks;
			nPeaks = 0;
		}
		if (isdigit(p[0]))
			nPeaks++;
	}
	printf("done (contains %i MS2 spectra)\n", datasetAsize);
	fflush(stdout);
	fclose(datasetA);

	if ((topN > -1) && (topN < datasetAsize)) {
		printf("filtering top-%i spectra...", topN);
		fflush(stdout);
		datasetAIntensities = malloc(datasetAsize * sizeof(double));
		if ((datasetA = fopen(datasetAFilename, "r")) == NULL) {
			printf("error opening dataset A %s for reading", datasetAFilename);
			return -1;
		}
		i = -1;
		while (fgets(line, 8192, datasetA) != NULL) {
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

		datasetACutoff = quickSelect(datasetAIntensities, 0, i, topN); // quickselect top-Nth intensity
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
	datasetBsize = 0;
	nPeaks = 0;
	while (fgets(line, 8192, datasetB) != NULL) {
		if (strcmp(line, "\n") == 0)
			continue;
		p = strtok(line, "=");
		if (strcmp("TITLE", p) == 0) {
			datasetBsize++;
			if (nPeaks > maxPeaks)
				maxPeaks = nPeaks;
			nPeaks = 0;
		}
	}
	printf("done (contains %i MS2 spectra)\n", datasetBsize);
	fflush(stdout);
	fclose(datasetB);

	if ((topN > -1) && (topN < datasetBsize)) {
		printf("filtering top-%i spectra...", topN);
		fflush(stdout);
		datasetBIntensities = malloc(datasetBsize * sizeof(double));
		if ((datasetB = fopen(datasetBFilename, "r")) == NULL) {
			printf("error opening dataset B %s for reading", datasetBFilename);
			return -1;
		}
		i = -1;
		while (fgets(line, 8192, datasetB) != NULL) {
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
	A = (DatasetType*) malloc(datasetAsize * sizeof(DatasetType));
	B = (DatasetType*) malloc(datasetBsize * sizeof(DatasetType));

	/* read in tandem mass spectra from MGF files */

	printf("done\nreading %i MS2 spectra from %s...", datasetAsize,
			datasetAFilename);
	fflush(stdout);
	if ((datasetA = fopen(datasetAFilename, "r")) == NULL) {
		printf("error opening dataset A %s for reading", datasetAFilename);
		return -1;
	}
	i = 0;
	j = 0;
	while (fgets(line, 8192, datasetA) != NULL) {
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
		if (strspn("SCANS", p) > 4) {
			A[i].scan = (long) atol0(strpbrk(p, "0123456789"));
			continue;
		}
		if (strcmp("END", p) == 0) {
			A[i].nPeaks = j;
			i++;
			j = 0;
		}
	}

	printf("done\nreading %i MS2 spectra from %s...", datasetBsize,
			datasetBFilename);
	fflush(stdout);
	if ((datasetB = fopen(datasetBFilename, "r")) == NULL) {
		printf("error opening dataset B %s for reading", datasetBFilename);
		return -1;
	}
	i = 0;
	j = 0;
	while (fgets(line, 8192, datasetB) != NULL) {
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
		if (strspn("SCANS", p) > 4) {
			B[i].scan = (long) atol0(strpbrk(p, "0123456789"));
			continue;
		}
		if (strcmp("END", p) == 0) {
			B[i].nPeaks = j;
			i++;
			j = 0;
		}
	}
	printf("done\n");
	fflush(stdout);

	printf("scaling, normalizing and binning %i MS2 spectra from %s..",
			datasetAsize, datasetAFilename);
	fflush(stdout);
	for (j = 0; j < datasetAsize; j++) {
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

	printf(".done\nscaling, normalizing and binning %i MS2 spectra from %s..",
			datasetBsize, datasetBFilename);
	fflush(stdout);
	for (j = 0; j < datasetBsize; j++) {
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

	/* go through spectra (entries) in MGF file 1 and compare with those in MGF file 2 and vice versa */

	printf(".done\nmatching spectra and computing set distance.");
	fflush(stdout);
	dotProdSum = 0.0;
	nComparisons = 0;
	greaterThanCutoff = 0;
	setSAB = 0;
	setSBA = 0;
	datasetAActualCompared = 0;
	datasetBActualCompared = 0;
	for (i = 0; i < 200; i++)
		histogram[i] = 0;

	for (i = 0; i < datasetAsize; i++) {
		if (topN > -1)
			if (topN < datasetAsize)
				if (datasetAIntensities[i] <= datasetACutoff)
					continue; /* compare only top-N spectra (including spectra of same intensity as Nth spectrum) */
		/* if(A[i].scan<start_scan) continue; */
		/* if(A[i].scan>end_scan) continue; */

		maxDotProd = 0.0;
		datasetAActualCompared++;

		for (j = 0; j < datasetBsize; j++) {
			if (topN > -1)
				if (topN < datasetBsize)
					if (datasetBIntensities[j] <= datasetBCutoff)
						continue;
			if ((A[i].scan - B[j].scan) > maxScanNumberDifference)
				continue;
			if ((B[j].scan - A[i].scan) > maxScanNumberDifference)
				break;
			if (fabs(B[j].precursorMz - A[i].precursorMz)
					< maxPrecursorDifference) {
				dotProd = 0;
				for (k = 0; k < nBins; k++)
					dotProd += A[i].bin[k] * B[j].bin[k];
				if (fabs(dotProd) <= 1.00)
					histogram[100 + (int) floor(dotProd * 99.999999999999)]++;
				nComparisons++;
				if (dotProd > maxDotProd) {
					maxDotProd = dotProd;
				}
			}
		}
		if (maxDotProd > cutoff)
			greaterThanCutoff++;
		if (maxDotProd > cutoff)
			setSAB++;
	}

	printf(".");
	fflush(stdout);

	for (i = 0; i < datasetBsize; i++) {
		if (topN > -1)
			if (topN < datasetBsize)
				if (datasetBIntensities[i] <= datasetBCutoff)
					continue;
		/* if(B[i].scan<startScan) continue; */
		/* if(B[i].scan>end_scan) continue; */
		maxDotProd = 0.0;
		datasetBActualCompared++;
		for (j = 0; j < datasetAsize; j++) {
			if (topN > -1)
				if (topN < datasetAsize)
					if (datasetAIntensities[j] <= datasetACutoff)
						continue;
			if ((B[i].scan - A[j].scan) > maxScanNumberDifference)
				continue;
			if ((A[j].scan - B[i].scan) > maxScanNumberDifference)
				break;
			if (fabs(A[j].precursorMz - B[i].precursorMz)
					< maxPrecursorDifference) {
				dotProd = 0;
				for (k = 0; k < nBins; k++)
					dotProd += B[i].bin[k] * A[j].bin[k];
				if (fabs(dotProd) <= 1.00)
					histogram[100 + (int) floor(dotProd * 99.999999999999)]++;
				/* nComparisons++; */
				if (dotProd > maxDotProd) {
					maxDotProd = dotProd;
				}
			}
		}
		if (maxDotProd > cutoff)
			greaterThanCutoff++; /* counting shared spectra from both datasets */
		if (maxDotProd > cutoff)
			setSBA++;
	}

	printf(
			".done (compared %i spectra from dataset A with %i spectra from dataset B)\nwriting results to file...",
			datasetAActualCompared, datasetBActualCompared);

	/* print output to file */

	if ((output = fopen(outputFilename, "w")) == NULL) {
		printf("error opening output file %s for writing", outputFilename);
		return -1;
	}
	fprintf(output, "dataset_A\t%s\n", datasetAFilename);
	fprintf(output, "dataset_B\t%s\n", datasetBFilename);
	if (metric == 0)
		fprintf(output, "set_distance\t%1.10f\n",
				(double) greaterThanCutoff / nComparisons); /* original metric */
	if (metric == 1)
		fprintf(output, "set_distance\t%1.10f\n",
				(double) greaterThanCutoff / (datasetAsize + datasetBsize)); /* symmetric metric */
	if (metric == 2)
		fprintf(output, "set_distance\t%1.10f\n",
				(double) setSAB / (2 * datasetAsize)
						+ (double) setSBA / (2 * datasetBsize)); /* compareMS2 2.0 symmetric metric */
	fprintf(output, "set_metric\t%i\n", metric);
	fprintf(output,
			"scan_range\t%i\t%i\nmax_scan_diff\t%.5f\nmax_m/z_diff\t%.5f\nscaling_power\t%.5f\nnoise_threshold\t%.5f\n",
			startScan, endScan, maxScanNumberDifference, maxPrecursorDifference,
			scaling, noise);
	if (qc == 0)
		fprintf(output, "dataset_A_QC\t%.4f\n", (float) datasetAsize);
	if (qc == 0)
		fprintf(output, "dataset_B_QC\t%.4f\n", (float) datasetBsize);
	fprintf(output, "n_gt_cutoff\t%i\n", greaterThanCutoff);
	fprintf(output, "n_comparisons\t%i\n", nComparisons);
	fprintf(output, "min_peaks\t%i\n", minPeaks);
	fprintf(output, "max_peaks\t%i\n", maxPeaks);
	fprintf(output, "m/z_range\t%.4f\t%.4f\n", minMz, maxMz);
	fprintf(output, "m/z_bin_size\t%.4f\n", binSize);
	fprintf(output, "n_m/z_bins\t%i\n", nBins);
	/* fprintf(output,"histogram (interval, midpoint, comparisons)\n"); */
	for (i = 0; i < 200; i++)
		fprintf(output, "histogram\t%1.3f\t%1.3f\t%1.3f\t%i\n",
				(double) (i - 100) / 100, (double) (i + 1 - 100) / 100,
				(double) (i + 0.5 - 100) / 100, histogram[i]);
	fflush(output);

	/* close files and free memory */

	fclose(output);
	printf("done\nfreeing memory...");
	fflush(stdout);
	free(A);
	free(B);
	printf("done\n");
	fflush(stdout);

	/* return from main */

	return 0;
}
