#include <stdio.h>
#include <string.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>

typedef struct moleculeIndex
{
	int molID, start, end;
} MOLECULEINDEX;

typedef struct trajectory
{
	int atomID, atomType, molType, molID, ix, iy, iz;
	float x, y, z;
	int isEndGroup;
	int clusterID; // this cluster ID is from DBSCAN
	bool core, edge;
} TRAJECTORY;

typedef struct vector
{
	float x1, y1, z1;
	float x2, y2, z2;
	float xc, yc, zc;
} VECTOR;

typedef struct datafileInfo
{
	int nAtoms, nBonds, nAngles, nDihedrals, nImpropers;
	int nAtomTypes, nBondTypes, nAngleTypes, nDihedralTypes, nImproperTypes;
} DATAFILE_INFO;

typedef struct datafile_atoms
{
	int resNumber;
	char resName[6], atomName[6], atomType2[6], molName[6];

	int id, molType, atomType;
	float charge, x, y, z;
} DATA_ATOMS;

typedef struct datafile_bonds
{
	int id, bondType, atom1, atom2;
} DATA_BONDS;

typedef struct datafile_angles
{
	int id, angleType, atom1, atom2, atom3;
} DATA_ANGLES;

typedef struct datafile_dihedrals
{
	int id, dihedralType, atom1, atom2, atom3, atom4;
} DATA_DIHEDRALS;

typedef struct datafile_impropers
{
	int id, improperType, atom1, atom2, atom3, atom4;
} DATA_IMPROPERS;

typedef struct simulationBoundary
{
	float xlo, xhi, ylo, yhi, zlo, zhi;
	float xLength, yLength, zLength;
} SIMULATION_BOUNDARY;

typedef struct rdf
{
	float rlo, rhi, gofr;
} RDF;

typedef struct stats
{
	float average, standardDeviation;
} STATS;

typedef struct orderParameterBins
{
	float orderParameter, rlo, rhi, count;
} ORDERPARAMETER_BINS;

TRAJECTORY *readTimestep (FILE *file_dump, TRAJECTORY *atoms, int nAtomEntries, SIMULATION_BOUNDARY *boundary)
{
	char lineString[2000];
	int currentAtomID = 1;

	for (int i = 0; i < 5; ++i) {
		fgets (lineString, 2000, file_dump); }

	fgets (lineString, 2000, file_dump); sscanf (lineString, "%f %f\n", &(*boundary).xlo, &(*boundary).xhi);
	fgets (lineString, 2000, file_dump); sscanf (lineString, "%f %f\n", &(*boundary).ylo, &(*boundary).yhi);
	fgets (lineString, 2000, file_dump); sscanf (lineString, "%f %f\n", &(*boundary).zlo, &(*boundary).zhi);
	fgets (lineString, 2000, file_dump);

	(*boundary).xLength = (*boundary).xhi - (*boundary).xlo;
	(*boundary).yLength = (*boundary).yhi - (*boundary).ylo;
	(*boundary).zLength = (*boundary).zhi - (*boundary).zlo;

	for (int i = 0; i < nAtomEntries; ++i)
	{
		fgets (lineString, 2000, file_dump);
		sscanf (lineString, "%d\n", &currentAtomID);

		if (currentAtomID > 0)
		{
			sscanf (lineString, "%d %d %d %d %f %f %f %d %d %d\n", &atoms[currentAtomID - 1].atomID, &atoms[currentAtomID - 1].atomType, &atoms[currentAtomID - 1].molID, &atoms[currentAtomID - 1].molType, &atoms[currentAtomID - 1].x, &atoms[currentAtomID - 1].y, &atoms[currentAtomID - 1].z, &atoms[currentAtomID - 1].ix, &atoms[currentAtomID - 1].iy, &atoms[currentAtomID - 1].iz);
			// sscanf (lineString, "%d %d %f %f %f %d %d %d\n", &atoms[currentAtomID - 1].atomID, &atoms[currentAtomID - 1].atomType, &atoms[currentAtomID - 1].x, &atoms[currentAtomID - 1].y, &atoms[currentAtomID - 1].z, &atoms[currentAtomID - 1].ix, &atoms[currentAtomID - 1].iy, &atoms[currentAtomID - 1].iz);
		}

		atoms[currentAtomID - 1].isEndGroup = 0;
	}

	return atoms;
}

SIMULATION_BOUNDARY readDumpBoundary (FILE *file_dump, SIMULATION_BOUNDARY boundary)
{
	rewind (file_dump);
	char lineString[2000];

	for (int i = 0; i < 5; ++i)
	{
		fgets (lineString, 2000, file_dump);
	}

	fgets (lineString, 2000, file_dump);
	sscanf (lineString, "%f %f\n", &boundary.xlo, &boundary.xhi);
	fgets (lineString, 2000, file_dump);
	sscanf (lineString, "%f %f\n", &boundary.ylo, &boundary.yhi);
	fgets (lineString, 2000, file_dump);
	sscanf (lineString, "%f %f\n", &boundary.zlo, &boundary.zhi);
	rewind (file_dump);

	boundary.xLength = boundary.xhi - boundary.xlo;
	boundary.yLength = boundary.yhi - boundary.ylo;
	boundary.zLength = boundary.zhi - boundary.zlo;

	return boundary;
}

DATAFILE_INFO readData (const char *dataFileName, DATA_ATOMS **atoms, DATA_BONDS **bonds, DATA_ANGLES **angles, DATA_DIHEDRALS **dihedrals, DATA_IMPROPERS **impropers)
{
	printf("Reading LAMMPS data file...\n");
	FILE *input;
	input = fopen (dataFileName, "r");

	int isAtomLine = 0, /*nAtoms = -1,*/ nAtomLine = 0;
	int isBondLine = 0, /*nBonds = -1,*/ nBondLine = 0;
	int isAngleLine = 0, /*nAngles = -1,*/ nAngleLine = 0;
	int isDihedralLine = 0, /*nDihedrals = -1,*/ nDihedralLine = 0;
	int isImproperLine = 0, /*nImpropers = -1,*/ nImproperLine = 0;
	int printHeaderInfo = 1;

	DATAFILE_INFO datafile;
	datafile.nAtoms = -1;
	datafile.nBonds = -1;
	datafile.nAngles = -1;
	datafile.nDihedrals = -1;
	datafile.nImpropers = -1;

	char lineString[1000];

	// DATA_ATOMS *atoms;
	// DATA_BONDS *bonds;
	// DATA_ANGLES *angles;
	// DATA_DIHEDRALS *dihedrals;
	// DATA_IMPROPERS *impropers;
	*atoms = NULL;
	*bonds = NULL;
	*angles = NULL;
	*dihedrals = NULL;
	*impropers = NULL;

	while ((fgets (lineString, 1000, input) != NULL))
	{
		if (strstr (lineString, "atoms"))
		{
			sscanf (lineString, "%d \n", &datafile.nAtoms);
			(*atoms) = (DATA_ATOMS *) malloc (datafile.nAtoms * sizeof (DATA_ATOMS));
		}

		if (strstr (lineString, "bonds"))
		{
			sscanf (lineString, "%d \n", &datafile.nBonds);
			(*bonds) = (DATA_BONDS *) malloc (datafile.nBonds * sizeof (DATA_BONDS));
		}

		if (strstr (lineString, "angles"))
		{
			sscanf (lineString, "%d \n", &datafile.nAngles);
			(*angles) = (DATA_ANGLES *) malloc (datafile.nAngles * sizeof (DATA_ANGLES));
		}

		if (strstr (lineString, "dihedrals"))
		{
			sscanf (lineString, "%d \n", &datafile.nDihedrals);
			(*dihedrals) = (DATA_DIHEDRALS *) malloc (datafile.nDihedrals * sizeof (DATA_DIHEDRALS));
		}

		if (strstr (lineString, "impropers"))
		{
			sscanf (lineString, "%d \n", &datafile.nImpropers);
			(*impropers) = (DATA_IMPROPERS *) malloc (datafile.nImpropers * sizeof (DATA_IMPROPERS));
		}

		if (strstr (lineString, "atom types"))
			sscanf (lineString, "%d \n", &datafile.nAtomTypes);

		if (strstr (lineString, "bond types"))
			sscanf (lineString, "%d \n", &datafile.nBondTypes);

		if (strstr (lineString, "angle types"))
			sscanf (lineString, "%d \n", &datafile.nAngleTypes);

		if (strstr (lineString, "dihedral types"))
			sscanf (lineString, "%d \n", &datafile.nDihedralTypes);

		if (strstr (lineString, "improper types"))
			sscanf (lineString, "%d \n", &datafile.nImproperTypes);

		if ((datafile.nAtoms >= 0) && (datafile.nBonds >= 0) && (datafile.nAngles >= 0) && (datafile.nDihedrals >= 0) && (datafile.nImpropers >= 0) && (printHeaderInfo))
			printHeaderInfo = 0;

		if (strstr (lineString, "Atoms"))
		{
			isAtomLine = 1;
			fgets (lineString, 1000, input);
			fgets (lineString, 1000, input);
		}

		if (strstr (lineString, "Bonds"))
		{
			isBondLine = 1;
			fgets (lineString, 1000, input);
			fgets (lineString, 1000, input);
		}

		if (strstr (lineString, "Angles"))
		{
			isAngleLine = 1;
			fgets (lineString, 1000, input);
			fgets (lineString, 1000, input);
		}

		if (strstr (lineString, "Dihedrals"))
		{
			isDihedralLine = 1;
			fgets (lineString, 1000, input);
			fgets (lineString, 1000, input);
		}

		if (strstr (lineString, "Impropers"))
		{
			isImproperLine = 1;
			fgets (lineString, 1000, input);
			fgets (lineString, 1000, input);
		}

		if (isAtomLine)
		{
			sscanf (lineString, "%d %d %d %f %f %f %f\n", 
				&(*atoms)[nAtomLine].id, 
				&(*atoms)[nAtomLine].molType, 
				&(*atoms)[nAtomLine].atomType, 
				&(*atoms)[nAtomLine].charge, 
				&(*atoms)[nAtomLine].x, 
				&(*atoms)[nAtomLine].y, 
				&(*atoms)[nAtomLine].z);
			nAtomLine++;
			if (nAtomLine == datafile.nAtoms)
				isAtomLine = 0;
		}

		if (isBondLine)
		{
			sscanf (lineString, "%d %d %d %d\n", 
				&(*bonds)[nBondLine].id, 
				&(*bonds)[nBondLine].bondType, 
				&(*bonds)[nBondLine].atom1, 
				&(*bonds)[nBondLine].atom2);
			nBondLine++;
			if (nBondLine == datafile.nBonds)
				isBondLine = 0;
		}

		if (isAngleLine)
		{
			sscanf (lineString, "%d %d %d %d %d\n", 
				&(*angles)[nAngleLine].id, 
				&(*angles)[nAngleLine].angleType, 
				&(*angles)[nAngleLine].atom1, 
				&(*angles)[nAngleLine].atom2, 
				&(*angles)[nAngleLine].atom3);
			nAngleLine++;
			if (nAngleLine == datafile.nAngles)
				isAngleLine = 0;
		}

		if (isDihedralLine)
		{
			sscanf (lineString, "%d %d %d %d %d %d\n", 
				&(*dihedrals)[nDihedralLine].id, 
				&(*dihedrals)[nDihedralLine].dihedralType, 
				&(*dihedrals)[nDihedralLine].atom1, 
				&(*dihedrals)[nDihedralLine].atom2, 
				&(*dihedrals)[nDihedralLine].atom3, 
				&(*dihedrals)[nDihedralLine].atom4);
			nDihedralLine++;
			if (nDihedralLine == datafile.nDihedrals)
				isDihedralLine = 0;
		}

		if (isImproperLine)
		{
			sscanf (lineString, "%d %d %d %d %d %d\n", 
				&(*impropers)[nImproperLine].id, 
				&(*impropers)[nImproperLine].improperType, 
				&(*impropers)[nImproperLine].atom1, 
				&(*impropers)[nImproperLine].atom2, 
				&(*impropers)[nImproperLine].atom3, 
				&(*impropers)[nImproperLine].atom4);
			nImproperLine++;
			if (nImproperLine == datafile.nImpropers)
				isImproperLine = 0;
		}
	}

	printf("\nFrom input data file:\n\n nAtoms: %d\n nBonds: %d\n nAngles: %d\n nDihedrals: %d\n nImpropers: %d\n\n", datafile.nAtoms, datafile.nBonds, datafile.nAngles, datafile.nDihedrals, datafile.nImpropers);

	return datafile;
}

TRAJECTORY *initializeAtoms (TRAJECTORY *atoms, int nAtoms)
{
	for (int i = 0; i < nAtoms; ++i)
	{
		atoms[i].atomID = 0;
		atoms[i].atomType = 0;
		atoms[i].molType = 0;
		atoms[i].molID = 0;
		atoms[i].ix = 0;
		atoms[i].iy = 0;
		atoms[i].iz = 0;
		atoms[i].x = 0;
		atoms[i].y = 0;
		atoms[i].z = 0;
		atoms[i].isEndGroup = 0;
		atoms[i].clusterID = 0;
		atoms[i].core = false;
		atoms[i].edge = false;
	}

	return atoms;
}

TRAJECTORY *initializeMolID (TRAJECTORY *atoms, int nAtoms)
{
	for (int i = 0; i < nAtoms + 1; ++i)
	{
		atoms[i].molType = 0;
	}
	return atoms;
}

TRAJECTORY *assignMolID (TRAJECTORY *atoms, int nAtoms)
{
	int nAtomsInbetween = 0, IDofPreviousBr = 0;
	int nCTAB = 0, nDDAB = 0, nSurfactants = 0;

	for (int i = 0; i < nAtoms; ++i)
	{
		if (atoms[i].atomType == 5)
		{
			nAtomsInbetween = atoms[i].atomID - IDofPreviousBr;
			IDofPreviousBr = atoms[i].atomID;
			nSurfactants++;

			if (nAtomsInbetween == 63)
			{
				nCTAB++;
				atoms[i].molType = 1;
				for (int j = i; j > (i - 63); --j)
				{
					if (atoms[j].molType == 0) {
						atoms[j].molType = 1; }
				}
			}
			else if (nAtomsInbetween == 84)
			{
				nDDAB++;
				atoms[i].molType = 2;
				for (int j = i; j > (i - 84); --j)
				{
					if (atoms[j].molType == 0) {
						atoms[j].molType = 2; }
				}
			}

			nAtomsInbetween = 0;
		}
	}

	return atoms;
}

int countNAtoms (FILE *file_dump, int *nAtomEntries)
{
	int nAtoms, currentAtomID, nAtomsFixed;
	char lineString[2000];
	rewind (file_dump);

	for (int i = 0; i < 4; ++i) {
		fgets (lineString, 2000, file_dump); }

	sscanf (lineString, "%d\n", &nAtoms);
	(*nAtomEntries) = nAtoms;
	rewind (file_dump);
	nAtomsFixed = nAtoms;

	for (int i = 0; i < 9; ++i) {
		fgets (lineString, 2000, file_dump); }

	for (int i = 0; i < nAtoms; ++i)
	{
		fgets (lineString, 2000, file_dump);
		sscanf (lineString, "%d\n", &currentAtomID);

		if (currentAtomID > nAtoms) {
			nAtomsFixed = currentAtomID; }
	}

	return nAtomsFixed;
}

int countNMolecules (TRAJECTORY *atoms, int nMolecules, int nAtoms)
{
	nMolecules = 0;

	for (int i = 0; i < nAtoms; ++i)
	{
		if (nMolecules < atoms[i].molID) {
			nMolecules = atoms[i].molID; }
	}

	return nMolecules;
}

MOLECULEINDEX *getIndexOfMolecules (TRAJECTORY *atoms, int nAtoms, int nMolecules, MOLECULEINDEX *indexOfMolecules)
{
	int previousMolID = 0, currentMolID;

	for (int i = 0; i < nAtoms; ++i)
	{
		currentMolID = atoms[i].molID;

		if (currentMolID != previousMolID)
		{
			indexOfMolecules[currentMolID - 1].molID = currentMolID;
			indexOfMolecules[currentMolID - 1].start = i;

			if (currentMolID > 1) {
				indexOfMolecules[currentMolID - 2].end = (i - 1); }
		}

		previousMolID = currentMolID;
	}

	return indexOfMolecules;
}

float translatePeriodicDistance (float x1, float x2, float simBoxLength, float newR)
{
	if (fabs (x1 - x2) > (simBoxLength / (float)2))
	{
		if (x1 >= x2) {
			newR = x1 - simBoxLength; }
		else if (x1 < x2) {
			newR = x1 + simBoxLength; }

		return newR;
	}
	else
	{
		return x1;
	}
}

int **getNeighbours (TRAJECTORY *atoms, int nAtomEntries, int **neighbourIDs, int maxNeighbors, SIMULATION_BOUNDARY boundary, float thresholdDistance, int atomType)
{
	float newX, newY, newZ;
	float distance;

	for (int i = 0; i < nAtomEntries; ++i)
	{
		if (atoms[i].atomID > 0 && atoms[i].atomType == atomType)
		{
			for (int j = 0; j < nAtomEntries; ++j)
			{
				if (atoms[j].atomID > 0 && i != j && atoms[i].atomType == atomType && atoms[i].molID != atoms[j].molID)
				{
					newX = translatePeriodicDistance (atoms[i].x, atoms[j].x, boundary.xLength, newX);
					newY = translatePeriodicDistance (atoms[i].y, atoms[j].y, boundary.yLength, newY);
					newZ = translatePeriodicDistance (atoms[i].z, atoms[j].z, boundary.zLength, newZ);

					distance = sqrt (
						(newX - atoms[j].x) * (newX - atoms[j].x) +
						(newY - atoms[j].y) * (newY - atoms[j].y) +
						(newZ - atoms[j].z) * (newZ - atoms[j].z)
						);

					if (distance < thresholdDistance)
					{
						for (int k = 0; k < maxNeighbors; ++k)
						{
							if (neighbourIDs[i][k] == -1)
							{
								neighbourIDs[i][k] = (j + 1);
								break;
							}
						}
					}
				}
			}
		}
	}

	return neighbourIDs;
}

// TRAJECTORY *markAtoms (TRAJECTORY *atoms, int nAtomEntries, float thresholdNeighbours, float thresholdDistance)
// {
// 	int maxNeighbors = ceil (thresholdNeighbours);

// 	for (int i = 0; i < nAtomEntries; ++i)
// 	{
		
// 	}

// 	return atoms;
// }

int **initNeighIDs (int **neighbourIDs, int nAtomEntries, int maxNeighbors)
{
	for (int i = 0; i < nAtomEntries; ++i)
	{
		for (int j = 0; j < maxNeighbors; ++j)
		{
			neighbourIDs[i][j] = -1;
		}
	}

	return neighbourIDs;
}

TRAJECTORY *markCoreAtoms (TRAJECTORY *atoms, int **neighbourIDs, int nAtomEntries, float thresholdNeighbours, float thresholdDistance, int maxNeighbors, int atomType, int currentTimeframe)
{
	int nNeighbors;

	for (int i = 0; i < nAtomEntries; ++i)
	{
		if (atoms[i].atomType == atomType)
		{
			nNeighbors = 0;

			for (int j = 0; j < maxNeighbors; ++j)
			{
				if (neighbourIDs[i][j] != -1) {
					nNeighbors++; }
			}

			if (nNeighbors >= thresholdNeighbours)
			{
				atoms[i].core = true;

				for (int j = 0; j < maxNeighbors; ++j)
				{
					if (neighbourIDs[i][j] > -1)
					{
						if (atoms[neighbourIDs[i][j] - 1].atomType == atomType)
						{
							atoms[neighbourIDs[i][j] - 1].edge = true;
						}
					}
				}
			}
		}
	}

	return atoms;
}

int countNCores (TRAJECTORY *atoms, int nAtomEntries, int nCore)
{
	nCore = 0;

	for (int i = 0; i < nAtomEntries; ++i)
	{
		if (atoms[i].core == true) {
			nCore++; }
	}

	return nCore;
}

int countNEdge (TRAJECTORY *atoms, int nAtomEntries, int nEdge)
{
	nEdge = 0;

	for (int i = 0; i < nAtomEntries; ++i)
	{
		if (atoms[i].edge == true && atoms[i].core == false) {
			nEdge++; }
	}

	return nEdge;
}

int *shuffleInt (int *array, int n)
{
	if (n > 1) 
	{
		int i;
		for (i = 0; i < n - 1; i++) 
		{
			int j = i + rand() / (RAND_MAX / (n - i) + 1);
			int t = array[j];
			array[j] = array[i];
			array[i] = t;
		}
	}

	return array;
}

TRAJECTORY *assignClusterID (TRAJECTORY *atoms, int nAtomEntries, int **neighbourIDs, int maxNeighbors, int atomType, MOLECULEINDEX *indexOfMolecules)
{
	int *randAtomIDs, currentID;
	randAtomIDs = (int *) malloc (nAtomEntries * sizeof (int));

	for (int i = 0; i < nAtomEntries; ++i)
	{
		randAtomIDs[i] = i;
	}

	// creating random number array to have random sequence for the clustering algorithm
	randAtomIDs = shuffleInt (randAtomIDs, nAtomEntries);

	// assign clusters sequentially, unique cluster ID to each atom
	for (int i = 0; i < nAtomEntries; ++i)
	{
		// currentID = randAtomIDs[i];

		if (atoms[i].atomID > 0 && atoms[i].atomType == atomType)
		{
			atoms[i].clusterID = i + 1;
			// printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\n", atoms[i].atomID, atoms[i].atomType, atoms[i].molID, atoms[i].molType, atoms[i].core, atoms[i].edge, atoms[i].clusterID);
			// usleep (100000);
		}
	}

	// reassign cluster ID based on molecule ID.

	// reassign cluster ID based on the nearest neighbour IDs.
	int nChanged = 0;

	do
	{
		nChanged = 0;

		for (int i = 0; i < nAtomEntries; ++i)
		{
			if (atoms[i].atomType == atomType)
			{
				for (int j = 0; j < maxNeighbors; ++j)
				{
					if (neighbourIDs[i][j] != -1 && atoms[neighbourIDs[i][j] - 1].clusterID < atoms[i].clusterID)
					{
						atoms[i].clusterID = atoms[neighbourIDs[i][j] - 1].clusterID;
						nChanged++;
					}
				}
			}
		}

	} while (nChanged > 0)

	return atoms;
}

int main(int argc, char const *argv[])
{
	if (argc != 7)
	{
		printf("INSUFFICIENT ARGUMENTS PASSED:\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n {~} argv[0] = program\n {~} argv[1] = input dump file\n {~} argv[2] = input data file\n {~} argv[3] = threshold distance\n {~} argv[4] = threshold neighbours\n {~} argv[5] = atom type to consider\n {~} argv[6] = maximum neighbours possible for each atom\n\n");
		exit (1);
	}

	FILE *file_dump, *file_data;
	file_dump = fopen (argv[1], "r");
	file_data = fopen (argv[2], "w");

	float thresholdDistance = atof (argv[3]), thresholdNeighbours = atof (argv[4]);
	int atomType = atoi (argv[5]);
	int maxNeighbors = atoi (argv[6]);

	int file_status, nAtoms, currentTimeframe = 0, nAtomEntries;
	nAtoms = countNAtoms (file_dump, &nAtomEntries);
	TRAJECTORY *atoms;
	atoms = (TRAJECTORY *) malloc (nAtoms * sizeof (TRAJECTORY));
	printf("Number of atoms in the trajectory file: %d\n", nAtoms);

	DATA_ATOMS *dataAtoms;
	DATA_BONDS *dataBonds;
	DATA_ANGLES *dataAngles;
	DATA_DIHEDRALS *dataDihedrals;
	DATA_IMPROPERS *dataImpropers;
	DATAFILE_INFO datafile;

	SIMULATION_BOUNDARY boundary;
	boundary = readDumpBoundary (file_dump, boundary);

	datafile = readData (argv[2], &dataAtoms, &dataBonds, &dataAngles, &dataDihedrals, &dataImpropers);

	atoms = initializeAtoms (atoms, nAtoms);

	int nMolecules = 0;
	MOLECULEINDEX *indexOfMolecules;

	rewind (file_dump);
	file_status = 1;

	int **neighbourIDs;
	neighbourIDs = (int **) malloc (nAtomEntries * sizeof (int *));

	for (int i = 0; i < nAtomEntries; ++i) {
		neighbourIDs[i] = (int *) malloc (maxNeighbors * sizeof (int)); }

	neighbourIDs = initNeighIDs (neighbourIDs, nAtomEntries, maxNeighbors);
	int nCore, nEdge;

	while (file_status != EOF)
	{
		// fprintf(stdout, "computing timestep: %d...                                                            \n", currentTimeframe);
		fflush (stdout);

		atoms = initializeAtoms (atoms, nAtoms);
		atoms = readTimestep (file_dump, atoms, nAtomEntries, &boundary);
		neighbourIDs = initNeighIDs (neighbourIDs, nAtomEntries, maxNeighbors);
		neighbourIDs = getNeighbours (atoms, nAtomEntries, neighbourIDs, maxNeighbors, boundary, thresholdDistance, atomType);

		// for (int i = 0; i < nAtomEntries; ++i)
		// {
		// 	printf("%d => ", (i + 1));
		// 	for (int j = 0; j < maxNeighbors; ++j)
		// 	{
		// 		printf("%d ", neighbourIDs[i][j]);
		// 	}
		// 	printf("\n");
		// 	usleep (100000);
		// }

		atoms = markCoreAtoms (atoms, neighbourIDs, nAtomEntries, thresholdNeighbours, thresholdDistance, maxNeighbors, atomType, currentTimeframe);
		nCore = countNCores (atoms, nAtomEntries, nCore);
		nEdge = countNEdge (atoms, nAtomEntries, nEdge);
		atoms = assignClusterID (atoms, nAtomEntries, neighbourIDs, maxNeighbors, atomType, indexOfMolecules);

		printf("edge: %d; core: %d\n", nEdge, nCore);
		usleep (100000);

		if (currentTimeframe == 0)
		{
			nMolecules = countNMolecules (atoms, nMolecules, nAtoms);
			indexOfMolecules = (MOLECULEINDEX *) malloc (nMolecules * sizeof (MOLECULEINDEX));
			indexOfMolecules = getIndexOfMolecules (atoms, nAtoms, nMolecules, indexOfMolecules);
		}

		// DBSCAN
		// Use the end of first peak in nearest neighbors and
		// use the end of first peak in rdf
		// for dbscan analysis.

		file_status = fgetc (file_dump);
		currentTimeframe++;
	}

	fclose (file_dump);
	fclose (file_data);
	return 0;
}