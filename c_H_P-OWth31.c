/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009,2010,2012, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#include <gromacs/copyrite.h>
#include <gromacs/filenm.h>
#include <gromacs/macros.h>
#include <gromacs/pbc.h>
#include <gromacs/smalloc.h>
#include <gromacs/statutil.h>
#include <gromacs/xvgr.h>
#include <gromacs/gmx_fatal.h>
#include <gromacs/rmpbc.h>
#include <gromacs/nbsearch.h>
#include <gromacs/trajana.h>
#include <gromacs/typedefs.h>
#include <gromacs/vec.h>
#include <gromacs/tpxio.h>
#include <math.h>
#include <gromacs/index.h>
#include <gromacs/typedefs.h>
/*! \brief
 * Template analysis data structure.
 */
#define delr   0.002
/* Threshold values are in nm */
#define thT_OW 0.600  /* Triad Heavy Atoms - OW */
#define thT_T  0.700  /* Triad Heavy Atoms - Triad Heavy Atoms */
#define thP_OW 0.310  /* Protein HA - OW. 1st solvation Minima */

#define thT_OWsq thT_OW * thT_OW
#define thT_Tsq  thT_T  * thT_T
#define thP_OWsq thP_OW * thP_OW


float distancesq(rvec [], matrix, int *, int *, int, int, float *, int *);
float distance(float *, float *, int);
float min_max(float *, int, float *, float *);
int  histogram(float *, int, float *, int );

int main(int argc, char *argv[]) {
   const char *desc[] = {
       "this is a small test program meant to serve as a template ", 
	   "when writing your own analysis tools. The advantage of this",
	   "is you can have access to the gromcas parameters in the",
	   "topology and also you can read the binary trajectories" 
	};
	
	int n = 1;
	static gmx_bool bPBC = FALSE;
	t_pargs pa[] = {
    	{ "-n", FALSE, etINT, {&n},
      	"Plot data for atom number n (starting on 1)" },
		{ "-pbc",      FALSE, etBOOL, {&bPBC},
		"Use periodic boundary conditions for computing distances" }
	};
	
	t_topology 		top;
	int        		ePBC;
	char       		title[STRLEN];
	t_trxframe 		fr;
	rvec       		*xtop;
	matrix     		box, box_pbc;
	int        		flags = TRX_READ_X;
	t_trxstatus 	*status;
	t_pbc       	pbc;
    t_atoms    		*atoms;
    int         	natoms;
	char       		*grpnm, *grpnmp;
    atom_id    		*index, *indexp;
    int         	i, nidx, nidxp;
    int        		v;
    int         	j, k;
    long         	***bin = (long ***)NULL;
	output_env_t    oenv;
	gmx_rmpbc_t     gpbc=NULL;

	t_filenm fnm[] = {
	    { efTPS,  NULL,  NULL, ffREAD },   /* this is for the topology */
    	{ efTRX, "-f", NULL, ffREAD }      /* and this for the trajectory */
  	};

	#define NFILE asize(fnm)

	parse_common_args(&argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME | PCA_TIME_UNIT | PCA_BE_NICE, NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv);


	read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&ePBC,&xtop,NULL,box,TRUE);
	sfree(xtop);

	n=n-1; /* Our enumeration started on 1, but C starts from 0 */
  	/* check that this atom exists */
  	if(n<0 || n>(top.atoms.nr)) {
	    printf("Error: Atom number %d is out of range.\n",n);
    	exit(1);
  	}
	
 /* RESINFO are meant for pdb;
 	printf("Atom name: %s\n",*(top.atoms.atomname[n]));
  	printf("Atom charge: %f\n",top.atoms.atom[n].q);
  	printf("Atom Mass: %f\n",top.atoms.atom[n].m);
  	printf("Atom resindex: %d\n",top.atoms.atom[n].resind); //Starts from 0
  	printf("Atom resname: %s\n",*(top.atoms.resinfo[top.atoms.atom[n].resind].name));
  	printf("Atom resno: %d %d\n",i, (top.atoms.resinfo[top.atoms.atom[i].resind].nr));
  	printf("Atom chain id: %d %c\n",i, top.atoms.resinfo[top.atoms.atom[i].resind].chainid); //Represents the pdb numbering
  	printf("Atom chain no: %d %d\n",i,top.atoms.resinfo[top.atoms.atom[i].resind].chainnum);
  	printf("Atom resno: %d %d\n",i, top.atoms.atom[i].resind);
*/
	atoms = &(top.atoms);

	/* The first time we read data is a little special */
	natoms = read_first_frame(oenv, &status, ftp2fn(efTRX, NFILE, fnm), &fr, flags);


	/*My Code Starts Here*/

	FILE  *P_OW;
    P_OW          = fopen("H_P-OWth31_0-250.xvg","wb");

    int tot_no_Amino_Acids = 0, d = 0, size = 99999;
    int pairsP_OW = 0, ProtMax = 0, ProteinH[50000] = {0};
    int uniqueP_OW = 0;
    int *storeP_OW, count = 0;
    int ProteinBegin = 0, ProteinEnd = 0, SOLBegin = 0, SOLEnd = 0;

    float distsq = 0.0;

    rvec B;

    char *AminoAcids[] = {"GLY", "ALA", "VAL", "LEU", "ILE", "PRO", "PHE", "TYR", "TRP",\
                          "SER", "THR", "CYS", "MET", "ASN", "GLN", "LYS", "ARG", "HIS",\
                          "ASP", "GLU", "HSD", "LYE"};

    storeP_OW = malloc(sizeof(int) * size);
    memset(storeP_OW, 0, sizeof(int) * size);

    k = 0;
    ProteinBegin = 0; ProteinEnd   = 3988;
    for (i = ProteinBegin; i < ProteinEnd; i++) {
        tot_no_Amino_Acids = (sizeof(AminoAcids)/sizeof(AminoAcids[0]));
        for( ; tot_no_Amino_Acids > 0; ) {
            if(strcmp(*(top.atoms.resinfo[top.atoms.atom[i].resind].name), AminoAcids[--tot_no_Amino_Acids]) == 0) {
                    ProteinH[k] = i; ProtMax = k; k++;
            }
        }
    }

    SOLBegin = 37818; SOLEnd   = 42234;

	do {
        /* Protein Water Contacts */
        memset(storeP_OW, 0, sizeof(int) * size);
        pairsP_OW = 0; uniqueP_OW = 0;

        for(i = 0; i <= ProtMax; i++) {
            for(j = SOLBegin; j <= SOLEnd; j++) {
                distsq = 0.0;
                for(d = 0; d < DIM; d++) {
                    /* Implementing the Minimum Image Distance Algorithm */
                    B[d] = fr.x[j][d];

                    if((fr.x[ProteinH[i]][d] - fr.x[j][d]) > (0.50 * fr.box[d][d]))
                        B[d] += fr.box[d][d];
                    if((fr.x[ProteinH[i]][d] - fr.x[j][d]) < -(0.50 * fr.box[d][d]))
                        B[d] -= fr.box[d][d];

                   distsq += (pow((fr.x[ProteinH[i]][d] - B[d]),2));
                }
                if(distsq < thP_OWsq) {
                    storeP_OW[pairsP_OW] = j;
                    pairsP_OW++;
                }
            }
        }
        /**************************************************************************/
        count = 0;
        for(j = 0; j < pairsP_OW; j++) {
            count = 0;
            for(i = j; i >= 0; i--) { 
                if(storeP_OW[j] == storeP_OW[i])
                    count++;
            }
            if(count < 2){
                uniqueP_OW++;
            }
        }
        fprintf(P_OW,"%8.3f %4d\n",fr.time, uniqueP_OW);
        /**************************************************************************/
        /**************************************************************************/
    }
    while (read_next_frame(oenv, status, &fr));	
    fclose(P_OW);
}

/*************************************** User Defined Functions ************************/
int histogram(float *dist, int pairs, float *hist, int maxbin) {
	int i = 0;
	for(i = 0; i < pairs; i++) {
		int bin = 0;
		bin = (int)(dist[i]/delr) + 1;

		if(bin < maxbin) 
			hist[bin] += 2;
	}
	return (0);
}

/*int jc = 0;
		for(i = 0; i < chainnr; i++) {
			for(j = 0; j < 40; j++) {
				printf("ATOM   %4d  Br  COM %d%4d%12.3f%8.3f%8.3f\n", count+1, i+1, j+1, 10*comg[i][jc][0],10* comg[i][jc][1],10* comg[i][jc][2]);
				count++; jc++;
			}
		}
*/
