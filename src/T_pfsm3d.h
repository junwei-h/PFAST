#ifndef _FIM3D_FUNCTIONS
#define _FIM3D_FUNCTIONS

#ifndef MAX
    #define MAX(a,b) (((a) > (b)) ? (a) : (b))
#endif

#ifndef MIN
    #define MIN(a,b) (((a) < (b)) ? (a) : (b))
#endif

#ifndef INF
    #define INF 9999999
#endif

#ifndef DIM
    #define DIM 3
#endif

/*#ifndef ROUND
#define ROUND(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))
#endif
*/
#ifndef PI
	#define PI (3.141592653589793)
#endif

#include <math.h>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <string.h>
//#include <algorithm>
//#include <vector>

using namespace std;

void readmod(double *mod, char *file, int nrows, int ncols, int nz)
{
	FILE *fp;
	if ((fp=fopen(file,"rb"))==NULL)
		printf(" The file %s cannot be opened\n",file);
	printf(" \n...reading file:\n\t %s\n",file);
	fread(mod, sizeof(double), nrows*ncols*nz, fp);
	fclose(fp);
}

void writemod(double *mod, char *file, int nrows, int ncols, int nz)
{
	FILE *fp;
	if ((fp=fopen(file,"wb"))==NULL)
		printf(" The file %s cannot be opened\n",file);
	printf(" \n...writing file:\n\t %s\n",file);
	fwrite(mod, sizeof(double), nrows*ncols*nz, fp);
	fclose(fp);
}

bool notin(int *node,int *CurrentNodes,int Num_CurrentNodes){
//to test if one node is in CurrentNodes
	
	//it = find (myvector.begin(), myvector.end(), 30);

	clock_t tic=clock();
	for (int i=0;i<Num_CurrentNodes;i++)
		if (CurrentNodes[i*DIM+0]==node[0] && CurrentNodes[i*DIM+1]==node[1] && CurrentNodes[i*DIM+2]==node[2])
			{if ((double)(clock()-tic)/CLOCKS_PER_SEC > 1.0 ) printf("notin function takes %8.5f seconds\n", (double)(clock()-tic)/CLOCKS_PER_SEC);
			return false;}
	if ((double)(clock()-tic)/CLOCKS_PER_SEC > 1.0 )  printf("notin function takes %8.5f seconds\n", (double)(clock()-tic)/CLOCKS_PER_SEC);
	return true;		
}

double eikslv3d(double *U,double *vp,int *node,double dx,double dy,double dz,int nrows, int ncols, int nz)
{	//to solve the time field U at node
	double q,a,b,c,ta,tb,tc,da=dx,db=dy,dc=dz,s=1.0/vp[node[2]*nrows*ncols+node[0]*ncols+node[1]];
	
	//a=min(Ux),b=min(Uy),c=min(Uz),s=slowness at node
	if (node[1]-1<0)              //node is at the left edge
		a=U[node[2]*nrows*ncols+node[0]*ncols+node[1]+1];
	else if (node[1]+1>=ncols)   //node is at the right edge
		a=U[node[2]*nrows*ncols+node[0]*ncols+node[1]-1];
	else                        //node is in the model
		a=MIN(U[node[2]*nrows*ncols+node[0]*ncols+node[1]+1],U[node[2]*nrows*ncols+node[0]*ncols+node[1]-1]);

	if (node[0]-1<0)              //node is at the upper edge
		b=U[node[2]*nrows*ncols+(node[0]+1)*ncols+node[1]];
	else if (node[0]+1>=nrows)    //node is at the lower edge
		b=U[node[2]*nrows*ncols+(node[0]-1)*ncols+node[1]];
	else							//node is in the model
		b=MIN(U[node[2]*nrows*ncols+(node[0]-1)*ncols+node[1]],U[node[2]*nrows*ncols+(node[0]+1)*ncols+node[1]]);
	
	if (node[2]-1<0)              //node is at the upper edge
		c=U[(node[2]+1)*nrows*ncols+node[0]*ncols+node[1]];
	else if (node[2]+1>=nz)    //node is at the lower edge
		c=U[(node[2]-1)*nrows*ncols+node[0]*ncols+node[1]];
	else							//node is in the model
		c=MIN(U[(node[2]-1)*nrows*ncols+node[0]*ncols+node[1]],U[(node[2]+1)*nrows*ncols+node[0]*ncols+node[1]]);

	//sort xa,xb,xc into c<=b<=a
	if (b>a) {swap(a,b);swap(da,db);}
	if (c>a) {swap(a,c);swap(da,dc);}
	if (c>b) {swap(b,c);swap(db,dc);} //to ensure a>=b>=c
	//double check
	if (!(c<=b && b<=a)) {printf("Error in sorting(c,b,a)=(%f,%f,%f)\n",c,b,a); return 0;}
	//////////////////////////////
	q=dc*s+c;
	if (q<=b) return q;
	
	ta=dc*dc+db*db;tb=-2.*(b*dc*dc+c*db*db);tc=dc*dc*b*b+db*db*c*c-db*db*dc*dc*s*s;
	if ((tb*tb-4.*ta*tc)>=0)
		q=(-tb+sqrt(tb*tb-4.*ta*tc))/2./ta;//q=(b+c+sqrt(-b*b-c*c+2*b*c+2*s*s))/2.0;
	else q=INF;
	if (q<=a) return q;
	
	ta=1.0/da/da+1.0/db/db+1.0/dc/dc;tb=-2.0*(a/da/da+b/db/db+c/dc/dc);tc=a*a/da/da+b*b/db/db+c*c/dc/dc-s*s;
	if ((tb*tb-4.*ta*tc)>=0)
		q=(-tb+sqrt(tb*tb-4.*ta*tc))/2./ta;//q=(2*(a+b+c)+sqrt(4*(a+b+c)*(a+b+c)-12*(a*a+b*b+c*c-s*s)))/6.0;
	else q=INF;
	return q;
}

bool isconv(double *p,double *U, int nrows, int ncols, int nz){
	for (int k=0;k<nz;k++)
		for (int i=0;i<nrows;i++)
			for (int j=0;j<ncols;j++)
				if (fabs(p[k*nrows*ncols+i*ncols+j]-U[k*nrows*ncols+i*ncols+j])>1e-9)
					return false;			
	return true;
}

void fsm3d(double *vp, double *U,int *src, int srcid, int MAXSwp, double dx,double dy,double dz, int nrows, int ncols, int nz,int pe){
	//printf("|\tPE_%d Called Eikonal Equation Solver for source #%d (%d,%d,%d)...\n",pe,srcid+1,src[srcid*DIM+0]+1,src[srcid*DIM+1]+1,src[srcid*DIM+2]+1);
	int count=0,OneNode[DIM];
	double *p=new double [nrows*ncols*nz],tmpU;
	if (!p) {printf("Memory Shotage\n");exit (1);}
	clock_t tic=clock();//start time recording
	
	//============Initialize Time Field for Point Source==========
	for (int kk=0;kk<nz;kk++){
		for (int ii=0;ii<nrows;ii++){
			for (int jj=0;jj<ncols;jj++){
				if (ii==src[srcid*DIM+0] && jj==src[srcid*DIM+1] && kk==src[srcid*DIM+2]) 
					U[kk*nrows*ncols+ii*ncols+jj]=0.0;
				else U[kk*nrows*ncols+ii*ncols+jj]=INF;
				p[kk*nrows*ncols+ii*ncols+jj]=0.0;			
			}
		}
	} 
	
	//////////////////////
	//Start Sweepting
	while (!isconv(p,U,nrows,ncols,nz) && count++<MAXSwp ){
		for (int kk=0;kk<nz;kk++)
			for (int ii=0;ii<nrows;ii++)
				for (int jj=0;jj<ncols;jj++)
					p[kk*nrows*ncols+ii*ncols+jj]=U[kk*nrows*ncols+ii*ncols+jj];		
			
		for (int s3=-1;s3<=1;s3+=2){
			for (int s1=-1;s1<=1;s1+=2){
				for (int s2=-1;s2<=1;s2+=2){
					//printf("PE_%d: Sweep (%d,%d,%d)...",pe,s3,s1,s2);clock_t tmp=clock();
					//printf("Sweep (%d,%d,%d)\n",s1,s2,s3);
					for (int k=(s3<0?nz-1:0);(s3<0?k>=0:k<nz);k+=s3){
						for (int i=(s1<0?nrows-1:0);(s1<0?i>=0:i<nrows);i+=s1){
							for (int j=(s2<0?ncols-1:0);(s2<0?j>=0:j<ncols);j+=s2){
								OneNode[0]=i;OneNode[1]=j;OneNode[2]=k;
								//if (notin(OneNode,src ,NumofSrc)){ //not a fixed node
								if (!(i==src[srcid*DIM+0] && j==src[srcid*DIM+1] && k==src[srcid*DIM+2])){
									tmpU=eikslv3d(U,vp,OneNode,dx,dy,dz,nrows,ncols,nz);
									U[k*nrows*ncols+i*ncols+j]=MIN(U[k*nrows*ncols+i*ncols+j],tmpU);
								}
							}
						}
					}
				//printf("takes %8.5f seconds\n", (double)(clock()-tmp)/CLOCKS_PER_SEC);
				}
			}
		}
	}
	delete [] p;p=NULL;
	//if (count<MAXSwp)
		//printf("|\tPE_%d: Eikonal Equation Solver Converged after %d sweeps, %8.5f seconds\n\n",pe, count,(double)(clock()-tic)/CLOCKS_PER_SEC);
	//else printf("|\tPE_%d: ....Warning! Eikonal Equation Solver  NOT Converged after %d sweeps, %8.5f seconds...\n",pe,MAXSwp,(double)(clock()-tic)/CLOCKS_PER_SEC);
	if (count>=MAXSwp) printf("|\tPE_%d: ....Warning! Eikonal Equation NOT Converged after %d sweeps, %8.5f seconds...\n",pe,MAXSwp,(double)(clock()-tic)/CLOCKS_PER_SEC);	
}

void usage()
{
	printf(" -----------------------------------------------------------------------------------------------------------------------\n");
	printf(" Wrong usage. See correct usage below:\n");
	printf(" mpirun -np number_of_processors T_PFAST3D input_file_name [Maximum Memory Number=5] [Tomography Scheme=3, Joint] [>output file]\"\n");
	printf("      Inverse Modeling usage 1: \"mpirun -np 10 ../bin/T_PFAST3D ./FSM3d01.inp\"\n");
	printf("      Inverse Modeling usage 2: \"mpirun -np 10 ../bin/T_PFAST3D ./FSM3d01.inp 5 3 > FSM3D01.out\"\n");
	printf("      Forward Modeling usage 3: \"mpirun -np 10 ../bin/T_PFAST3D_FW ./FSM3d01_FW.inp\"\n");
	printf("      Forward Modeling usage 4: \"mpirun -np 10 ../bin/T_PFAST3D_FW_NoS ./FSM3d01_FW.inp > FSM3D01.out\"\n");
	printf(" -----------------------------------------------------------------------------------------------------------------------\n");

}

void info() 
{
	printf(" *******************************\n");
	printf(" * This is Traveltime Tomography program: Version 1.0            \n");
	printf(" * Parallel 3-D Traveltime Tomography based on Fast Sweeping Method and Adjoint-state Technique      \n");
	printf(" *                                                           \n");
	printf(" * written by  J.W. Huang since April 2010                    \n");
	printf(" * See COPYING file for copying and redistribution conditions.\n");
	printf(" *******************************\n");
	printf(" \n");
}
#endif
