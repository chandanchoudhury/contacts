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
#include <stdbool.h>
//##include <gromacs/types/block.h>
/*! \brief
 * Template analysis data structure.
 */
#define delr   0.002
/* Threshold values are in nm */
#define thP_P 0.300  /* Protein HA - OW. 1st solvation Minima */

#define thP_Psq  thP_P  * thP_P


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
    atom_id    		*index, *indexp, ix, jx, **pairs;
    int         	i, nidx, nidxp;
    int        		v, kx;
    int         	j, k, *npairs;
    long         	***bin = (long ***)NULL;
	output_env_t    oenv;
	gmx_rmpbc_t     gpbc=NULL;
    t_blocka        *excl;
    bool            *bExcl, bNonSelfExcl;
    bool            *bAtomList;

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
    printf("\n");
    FILE  *P_P;
    P_P = fopen("P-Pth30_0-250.xvg","wb");


    int ProteinEnd = 3988, tot_no_Amino_Acids = 0, tot_conj_res = 0, d = 0;
    int *P_Parray, pairs_count = 0;

    float distsq = 0.0, L_dx = 0.0, half_Lx = 0.0;
    real dist_dx = 0.0;
    rvec B;
    char *AminoAcids[] = {"GLY", "ALA", "VAL", "LEU", "ILE", "PRO", "PHE", "TYR", "TRP",\
                          "SER", "THR", "CYS", "MET", "ASN", "GLN", "LYS", "ARG", "HIS",\
                          "ASP", "GLU", "HSD", "LYE"};
    char *Conjugated[] = {"PEGM", "MET2", "MAT", "MHT"};

    tot_conj_res = (sizeof(Conjugated)/sizeof(Conjugated[0]));

    natoms = ProteinEnd; excl = NULL;
    snew(bExcl, natoms); snew(index, natoms); snew(pairs, natoms);
    snew(bAtomList, natoms); snew(npairs, natoms);
    snew(P_Parray, natoms);

    for(j = 0; j < natoms; j++) {
        bExcl[j] = FALSE; bAtomList[j] = FALSE;
        index[j] = j;
    }

    /* get exclusions from topology */
    excl = &(top.excls);
    for(j = 0; j < natoms; j++)
        bExcl[j] = FALSE;
    
/* Indices of the protein atoms are stored as TRUE, when it belongs to the four
 * residues as listed in the Conjugated Array. Here TRUE sets the ix value to 0.
 * Rest other indices are FALSE i.e., stored as 1.*/
/* Store a value of 1 (TRUE) if the atom is 'Hydrogen' or belongs to 'Conjugated
 * Residues'.
 * #define FALSE (0)
 * #define TRUE  (1)
 * Above definition in gmx*/
    for(ix = 0; ix < natoms; ix++) {
        bAtomList[ix] = FALSE;
        if(strncmp(*(top.atoms.atomname[ix]), "H",1) == 0)
            bAtomList[ix] = TRUE;
        else {
            tot_conj_res = (sizeof(Conjugated)/sizeof(Conjugated[0]));
            for(; tot_conj_res > 0; )
                if(strcmp(*(top.atoms.resinfo[top.atoms.atom[ix].resind].name), Conjugated[--tot_conj_res]) == 0)
                    bAtomList[ix] = TRUE;
            }
    }
/*    for(ix = 0; ix < natoms; ix++)
    printf("Atom List :%6d %4s %4s %2d\n", ix, *(top.atoms.atomname[ix]), \
    *(top.atoms.resinfo[top.atoms.atom[ix].resind].name), bAtomList[ix]);
*/    
/*Pull out the Non-Bonded pairs.*/
    for(ix = 0; ix < natoms; ix++) {
        k = 0;
        if(!bAtomList[ix]) {
            k = 0;
            snew(pairs[ix], top.atoms.nr);
            for(j = 0; j < natoms; j++) bExcl[j] = FALSE;
            for(j = excl->index[ix]; j < excl->index[ix+1]; j++) {
                bExcl[excl->a[j]] = TRUE;
            }
            pairs[ix][0]=ix; 
            for(jx = ix; jx < natoms; jx++)
                if(!bAtomList[jx])
                    if(!bExcl[jx]) {
                        pairs[ix][++k]=jx; 
                    }
        }
        npairs[ix] = k;
    }
/*        for(ix = 0; ix < natoms; ix++) 
            for(jx = 0; jx < npairs[ix]; jx++) 
                printf("NBPairs :%6d %6d %6d\n", ix, jx, pairs[ix][jx]);
*//*        for(ix = 0; ix < natoms; ix++) 
           printf("npairs : %6d %6d\n", ix, npairs[ix]); 
*/
	do {
        snew(P_Parray, natoms); distsq = 0.0; L_dx = 0.0; half_Lx = 0.0;
        dist_dx = 0.0; pairs_count = 0;

        for(ix = 0; ix < natoms; ix++) {
            distsq = 0.0;
            for(jx = 1; jx < npairs[ix]; jx++) { /* [ix][0] is the first atom. This atom acts as reference to [ix][jx], jx != 0, atoms. So, we start from jx = 1, instead of jx=1.*/
                distsq = 0.0;
                for(d = 0; d < DIM; d++) {
                    B[d] = fr.x[pairs[ix][jx]][d];
                    L_dx = fr.x[pairs[ix][0]][d] - fr.x[pairs[ix][jx]][d];
                    half_Lx = 0.5 * fr.box[d][d];

                    if(L_dx > half_Lx)
                        B[d] += fr.box[d][d];
                    if(L_dx < (-1 * half_Lx))
                        B[d] -= fr.box[d][d];

                    dist_dx = fr.x[pairs[ix][0]][d] - B[d];
                    distsq += dist_dx * dist_dx;
                }

                if(distsq < thP_Psq) {
                    P_Parray[pairs[ix][0]]++;
                    P_Parray[pairs[ix][jx]]++;
//                printf("Within cutoff %4.1fnm Pairs: %6d %6d %6d %6d\n", 100*thP_P, pairs[ix][0], pairs[ix][jx], P_Parray[pairs[ix][0]], P_Parray[pairs[ix][jx]]);
                }
            }
        }
//        for(ix = 0; ix < natoms; ix++)
//            printf("Contacts(%4.1fnm): %6d %6d \n", 100*thP_P, ix, P_Parray[ix]);
        for(ix = 0; ix < natoms; ix++) {
            if(P_Parray[ix]) {
//                printf("%6d %6d\n", ix, P_Parray[ix]);
                pairs_count++;
            }
        }
        fprintf(P_P,"%8.3f %6d\n", fr.time, pairs_count);
    }
    while (read_next_frame(oenv, status, &fr));
    fclose(P_P);
}

/*************************************** User Defined Functions ************************/

/*int jc = 0;
		for(i = 0; i < chainnr; i++) {
			for(j = 0; j < 40; j++) {
				printf("ATOM   %4d  Br  COM %d%4d%12.3f%8.3f%8.3f\n", count+1, i+1, j+1, 10*comg[i][jc][0],10* comg[i][jc][1],10* comg[i][jc][2]);
				count++; jc++;
			}
		}
*/
