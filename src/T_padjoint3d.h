#ifndef _ADJOINT3D_METHODS
#define _ADJOINT3D_METHODS

#include "T_pfsm3d.h"
#include "fftn.h"
#include <cmath>

#ifndef POS
	#define POS(x) (((x)+fabs(x))/2.0)
#endif

#ifndef NEG
	#define NEG(x) (((x)-fabs(x))/2.0) 
#endif

#ifndef MAXV
	#define MAXV 8000 
#endif

#ifndef MINV
	#define MINV 1400  //speed of sound in 0 degree C water
#endif

/*--------FFT or IFFT Switch----------*/
#define FORWARD_SCALE	0.0
#define INVERSE_SCALE	-1.0

using namespace std;

void swap(double &a, double &b)
	{double dum=a; a=b; b=dum;}
	
int findin(int *node,int *CurrentNodes,int Num_CurrentNodes){
//to find the index of one node is in CurrentNodes
	//clock_t tic=clock();
	for (int i=0;i<Num_CurrentNodes;i++)
		if (CurrentNodes[i*DIM+0]==node[0] && CurrentNodes[i*DIM+1]==node[1] && CurrentNodes[i*DIM+2]==node[2])
			{//printf("findin function takes %8.5f seconds\n", (double)(clock()-tic)/CLOCKS_PER_SEC);
			return i;}
}

double adjoint_slv3d(double *U,double *T,int *OneNode,double dx,double dy,double dz, int nrows, int ncols, int nz, int *src, int srcid){
	double a1,a2,b1,b2,c1,c2,La1,La2,Lb1,Lb2,Lc1,Lc2;
	
	if (T[OneNode[2]*nrows*ncols+OneNode[0]*ncols+OneNode[1]]>=INF \
	|| T[OneNode[2]*nrows*ncols+OneNode[0]*ncols+OneNode[1]+1]>=INF || T[OneNode[2]*nrows*ncols+OneNode[0]*ncols+OneNode[1]-1]>=INF \
	|| T[OneNode[2]*nrows*ncols+(OneNode[0]+1)*ncols+OneNode[1]]>=INF || T[OneNode[2]*nrows*ncols+(OneNode[0]-1)*ncols+OneNode[1]]>=INF \
	|| T[(OneNode[2]+1)*nrows*ncols+OneNode[0]*ncols+OneNode[1]]>=INF || T[(OneNode[2]-1)*nrows*ncols+OneNode[0]*ncols+OneNode[1]]>=INF)
		return 0.0;
		
	a1=POS(-(T[OneNode[2]*nrows*ncols+OneNode[0]*ncols+OneNode[1]+1]-T[OneNode[2]*nrows*ncols+OneNode[0]*ncols+OneNode[1]])/dx);
	a2=NEG(-(T[OneNode[2]*nrows*ncols+OneNode[0]*ncols+OneNode[1]]-T[OneNode[2]*nrows*ncols+OneNode[0]*ncols+OneNode[1]-1])/dx);
	b1=POS(-(T[OneNode[2]*nrows*ncols+(OneNode[0]+1)*ncols+OneNode[1]]-T[OneNode[2]*nrows*ncols+OneNode[0]*ncols+OneNode[1]])/dy);
	b2=NEG(-(T[OneNode[2]*nrows*ncols+OneNode[0]*ncols+OneNode[1]]-T[OneNode[2]*nrows*ncols+(OneNode[0]-1)*ncols+OneNode[1]])/dy);
	c1=POS(-(T[(OneNode[2]+1)*nrows*ncols+OneNode[0]*ncols+OneNode[1]]-T[OneNode[2]*nrows*ncols+OneNode[0]*ncols+OneNode[1]])/dz);
	c2=NEG(-(T[OneNode[2]*nrows*ncols+OneNode[0]*ncols+OneNode[1]]-T[(OneNode[2]-1)*nrows*ncols+OneNode[0]*ncols+OneNode[1]])/dz);
	
	La1=POS(-(T[OneNode[2]*nrows*ncols+OneNode[0]*ncols+OneNode[1]]-T[OneNode[2]*nrows*ncols+OneNode[0]*ncols+OneNode[1]-1])/dx)*U[OneNode[2]*nrows*ncols+OneNode[0]*ncols+OneNode[1]-1];
	La2=NEG(-(T[OneNode[2]*nrows*ncols+OneNode[0]*ncols+OneNode[1]+1]-T[OneNode[2]*nrows*ncols+OneNode[0]*ncols+OneNode[1]])/dx)*U[OneNode[2]*nrows*ncols+OneNode[0]*ncols+OneNode[1]+1];
	Lb1=POS(-(T[OneNode[2]*nrows*ncols+OneNode[0]*ncols+OneNode[1]]-T[OneNode[2]*nrows*ncols+(OneNode[0]-1)*ncols+OneNode[1]])/dy)*U[OneNode[2]*nrows*ncols+(OneNode[0]-1)*ncols+OneNode[1]];
	Lb2=NEG(-(T[OneNode[2]*nrows*ncols+(OneNode[0]+1)*ncols+OneNode[1]]-T[OneNode[2]*nrows*ncols+OneNode[0]*ncols+OneNode[1]])/dy)*U[OneNode[2]*nrows*ncols+(OneNode[0]+1)*ncols+OneNode[1]];
	Lc1=POS(-(T[OneNode[2]*nrows*ncols+OneNode[0]*ncols+OneNode[1]]-T[(OneNode[2]-1)*nrows*ncols+OneNode[0]*ncols+OneNode[1]])/dz)*U[(OneNode[2]-1)*nrows*ncols+OneNode[0]*ncols+OneNode[1]];
	Lc2=NEG(-(T[(OneNode[2]+1)*nrows*ncols+OneNode[0]*ncols+OneNode[1]]-T[OneNode[2]*nrows*ncols+OneNode[0]*ncols+OneNode[1]])/dz)*U[(OneNode[2]+1)*nrows*ncols+OneNode[0]*ncols+OneNode[1]];
	
	if (fabs((a1-a2)/dx+(b1-b2)/dy+(c1-c2)/dz)>1e-9)
		return ((La1-La2)/dx+(Lb1-Lb2)/dy+(Lc1-Lc2)/dz)/((a1-a2)/dx+(b1-b2)/dy+(c1-c2)/dz);
	else //if (OneNode[0]==src[srcid*DIM+0] && OneNode[1]==src[srcid*DIM+1] && OneNode[2]==src[srcid*DIM+2]) //at source location
		return 0.0;
	//else printf("(%d,%d,%d):a1=%f,a2=%f,b1=%f,b2=%f,c1=%f,c2=%f (a1-a2+b1-b2+c1-c2)=%f\n",OneNode[0],OneNode[1],OneNode[2],a1,a2,b1,b2,c1,c2,a1-a2+b1-b2+c1-c2);
		
}

double adjoint_fsm3d(double *U,double *T,double *T_obs,int *rec, double *nxyz, int *src, int *srnum,int srcid,int MAXSwp,double dx,double dy,double dz,int nrows,int ncols,int nz,int pe){
	//printf("|\tPE_%d Called Adjoint Equation Solver for source #%d (%d,%d,%d)...\n",pe,srcid+1,src[srcid*DIM+0]+1,src[srcid*DIM+1]+1,src[srcid*DIM+2]+1);
	int count=0,OneNode[DIM],ch0=0,i,j,k;
	double nnx,nny,nnz,Tx,Ty,Tz,EE=0.0;
	double *p=new double[nrows*ncols*nz]; 
	int *rectag=new (nothrow)int[nrows*ncols*nz]; //tag 1 if it is a receiver location, and 0 otherwise
	
	if (!p) {printf("PE_%d: Memory shortage in adjoint_fsm3d\n",pe);exit (1);}
	clock_t toc=clock();//start time recording
	////////////////////////////
	//Calculate Normal Vector at the surface defined by receiver locations
	//surfnorm(rec,NumofRec,nxyz,simds);
	
	///////////////////////////
	//Initialize lambda (U) boundary condition
	//clock_t timeini=clock();
	//printf("PE_%d: Starts Boundary Initialization...\n",pe);
	for (int k=0;k<nrows*ncols*nz;k++){	U[k]=0.0;p[k]=0.0;	}
	for (int i=0;i<srcid;i++) ch0+=srnum[i];
	for (int k=0;k<nrows*ncols*nz;k++) rectag[k]=0;
			
	for (int idx=0;idx<srnum[srcid];idx++) {
		if (T_obs[idx+ch0]<100) { //ONLY when transmission time picked	
			i=rec[(idx+ch0)*DIM+0];j=rec[(idx+ch0)*DIM+1];k=rec[(idx+ch0)*DIM+2];
			rectag[k*nrows*ncols+i*ncols+j]=1;	
			//printf("PE_%d:(srcid,ch0,rec(i,j,k))=(%d,%d,rec(%d,%d,%d))\n",pe,srcid,ch0,i,j,k);
			EE+=(T_obs[idx+ch0]-T[k*nrows*ncols+i*ncols+j])*(T_obs[idx+ch0]-T[k*nrows*ncols+i*ncols+j])/2.0;
			if ( ! (rec[(idx+ch0)*DIM+0]==src[srcid*DIM+0] && rec[(idx+ch0)*DIM+1]==src[srcid*DIM+1] && rec[(idx+ch0)*DIM+2]==src[srcid*DIM+2]) )// not at a source location
			{
				nny=nxyz[(idx+ch0)*DIM+0];nnx=nxyz[(idx+ch0)*DIM+1];nnz=nxyz[(idx+ch0)*DIM+2];
				//One sided difference for travel time at the boundary
				if (k+1>=nz) 
					Tz=(T[k*nrows*ncols+i*ncols+j]-T[(k-1)*nrows*ncols+i*ncols+j])/dz;
				else
					Tz=(T[(k+1)*nrows*ncols+i*ncols+j]-T[k*nrows*ncols+i*ncols+j])/dz;
					
				if (i+1>=nrows)
					Ty=(T[k*nrows*ncols+i*ncols+j]-T[k*nrows*ncols+(i-1)*ncols+j])/dy;
				else
					Ty=(T[k*nrows*ncols+(i+1)*ncols+j]-T[k*nrows*ncols+i*ncols+j])/dy;
						
				if (j+1>=ncols)
					Tx=(T[k*nrows*ncols+i*ncols+j]-T[k*nrows*ncols+i*ncols+j-1])/dx;
				else
					Tx=(T[k*nrows*ncols+i*ncols+j+1]-T[k*nrows*ncols+i*ncols+j])/dx;
							
				if (fabs(T_obs[idx+ch0]-T[k*nrows*ncols+i*ncols+j])>1e-9 && fabs(nnx*Tx+nny*Ty+nnz*Tz)>1e-6)
					U[k*nrows*ncols+i*ncols+j]=(T_obs[idx+ch0]-T[k*nrows*ncols+i*ncols+j])/(nnx*Tx+nny*Ty+nnz*Tz);
				else
					U[k*nrows*ncols+i*ncols+j]=0.0;
				//if (pe==0) printf("%d rectag(%d %d %d)=%d\n",idx+1,i+1,j+1,k+1,rectag[k*nrows*ncols+i*ncols+j]);
				//if (pe==0) printf("%d %d %d %d %15.9f\t%15.9f\t%15.9f\t%15.9f\t%15.9f\t%15.9f\t%15.9f\t%15.9f\t%15.9f\n",idx+1,i+1,j+1,k+1,U[k*nrows*ncols+i*ncols+j],T_obs[idx+ch0],T[k*nrows*ncols+i*ncols+j],nnx,nny,nnz,Tx,Ty,Tz);
			}
		}
	}
	//printf("PE_%d: %8.5f seconds for Boundary Initialization\n",pe,(double)(clock()-timeini)/CLOCKS_PER_SEC);
	/////////////////Done Initialization////////////////////
	
	///////////////////////
	//Solve adjoint state	
	//while ( !ConvFlag && count++<10){
				
	while (count++<MAXSwp && !isconv(p,U,nrows,ncols,nz)  ){
		//clock_t tic=clock();
		for (int kk=0;kk<nrows*ncols*nz;kk++) p[kk]=U[kk];
		//printf("here\n");		
		for (int s3=-1;s3<=1;s3+=2){
			for (int s1=-1;s1<=1;s1+=2){
				for (int s2=-1;s2<=1;s2+=2){
					//printf("PE_%d: Sweep (%d,%d,%d)...",pe,s3,s1,s2);clock_t tmp=clock();
					for ( k=(s3<0?nz-2:1);(s3<0?k>=1:k<nz-1);k+=s3){
						for ( i=(s1<0?nrows-2:1);(s1<0?i>=1:i<nrows-1);i+=s1){
							for ( j=(s2<0?ncols-2:1);(s2<0?j>=1:j<ncols-1);j+=s2){
								OneNode[0]=i;OneNode[1]=j;OneNode[2]=k;
								//it=search (vector_rec.begin(), vector_rec.end(), OneNode, OneNode+DIM);
								//if (notin(OneNode,rec,NumofRec)){ //not a fixed node 	
								if (fabs((double)rectag[k*nrows*ncols+i*ncols+j])<1e-9) {
									//printf("[s1 s2]=[%d %d]\t[i j]=[%d %d]\t[p q]=[%f %f]\n",s1,s2,i,j,p,q);
									//U[i*ncols+j]=q;//MAX(p,q);
									U[k*nrows*ncols+i*ncols+j]=adjoint_slv3d(U,T,OneNode,dx,dy,dz,nrows,ncols,nz,src,srcid);

								}
							}
						}
					}
				//printf("takes %8.5f seconds\n", (double)(clock()-tmp)/CLOCKS_PER_SEC);
				}
			}
		}
		//printf("Sweep #%d...takes time %8.5f seconds\n",count,(double)(clock()-tic)/CLOCKS_PER_SEC);	
	}
	delete [] p;p=NULL;
	delete [] rectag;rectag=NULL;
	//if (count<MAXSwp)
	//	printf("|\tPE_%d: Adjoint Equation Converged after %d sweeps, %8.5f seconds\n\n",pe,count,(double)(clock()-toc)/CLOCKS_PER_SEC);
	//else printf("|\tPE_%d: ....Warning! Adjoint Equation NOT Converged after %d sweeps, %8.5f seconds...\n",pe,MAXSwp,(double)(clock()-toc)/CLOCKS_PER_SEC);
	if (count>=MAXSwp) printf("|\tPE_%d: ....Warning! Adjoint Equation NOT Converged after %d sweeps, %8.5f seconds...\n",pe,MAXSwp,(double)(clock()-toc)/CLOCKS_PER_SEC);
	
	return EE; //time residual
}

bool check_conv(double *dvp,double stepsize,double *vp,double eps,double *dEE, int Iterk, int nrows, int ncols,int nz){
	double mynorm=0,chgE;
	
	if (Iterk>=2) chgE=(dEE[Iterk-2]-dEE[Iterk-1])/dEE[Iterk-2];
	else chgE=1;
	
	for (int k=0;k<nrows*ncols*nz;k++){
		mynorm=(stepsize*stepsize*dvp[k]*dvp[k]/vp[k]/vp[k]>=mynorm)?\
		(stepsize*stepsize*dvp[k]*dvp[k]/vp[k]/vp[k]):mynorm;
				//mynorm=mynorm+dm[k]*dm[k]/vp[k]*vp[k];
	}
	mynorm=sqrt(mynorm);
	//if (mynorm<=eps || chgE <=1e-2) //Total update amount < eps or residual error less than dt ms
	if (chgE <=1e-2)
		return true;
	else if ((mynorm/nrows/ncols/nz)>=1) // average update is as twice as the previous velocity is considered as over corrected
		{printf("Over corrected..abort...\n");exit(1);}
	else
		return false;
}

void updatevp(double *vp,double *dm, double stepsize, int nrows, int ncols, int nz){
	for (int k=0;k<nrows*ncols*nz;k++){
		if (fabs(vp[k]-330)>1e-9) {
			vp[k]=vp[k]+stepsize*dm[k];
			if (vp[k]>MAXV) vp[k]=MAXV;
			if (vp[k]<MINV)	vp[k]=MINV;
		}
	}
}

void regu3d(double *Dj,double *vp,double nux, double nuy, double nuz, double dx,double dy,double dz, int nrows, int ncols, int nz){
	int ret,dims3[]={ncols,nrows,nz};
	double *rDj=new double [nrows*ncols*nz],*iDj=new double [nrows*ncols*nz];
	double *kx=new double [ncols],*ky=new double [nrows],*kz=new double [nz];
	double kxmax,kymax,kzmax,kxs,kys,kzs;
	
	for (int i=0;i<nrows*ncols*nz;i++) {
			rDj[i]=Dj[i];iDj[i]=0.0;
	}

	ret = fftn(3, dims3, rDj, iDj,  -1, FORWARD_SCALE); //combination of (-1, 0.0)==fftn.m in matlab
	if (ret) exit(1);
	
	kxmax=PI/dx;kymax=PI/dy;kzmax=PI/dz;kxs=2*PI/ncols/dx;kys=2*PI/nrows/dy;kzs=2*PI/nz/dz;
	
	for(int j=0;j<ncols;j++)
		kx[j]=-kxmax+kxs*j;
	for(int i=0;i<nrows;i++)	
		ky[i]=-kymax+kys*i;
	for(int k=0;k<nz;k++)	
		kz[k]=-kzmax+kzs*k;
		
	for(int j=0;j<ncols/2;j++)
		swap(kx[j],kx[ncols/2+j]);
	for(int i=0;i<nrows/2;i++)	
		swap(ky[i],ky[nrows/2+i]);
	for(int k=0;k<nz/2;k++)	
		swap(kz[k],kz[nz/2+k]);
	
	for (int k=0;k<nz;k++){	
		for (int i=0;i<nrows;i++){
			for (int j=0;j<ncols;j++){
				rDj[k*nrows*ncols+i*ncols+j]=rDj[k*nrows*ncols+i*ncols+j]/(1+nux*nux*kx[j]*kx[j]+nuy*nuy*ky[i]*ky[i]+nuz*nuz*kz[k]*kz[k]);
				iDj[k*nrows*ncols+i*ncols+j]=iDj[k*nrows*ncols+i*ncols+j]/(1+nux*nux*kx[j]*kx[j]+nuy*nuy*ky[i]*ky[i]+nuz*nuz*kz[k]*kz[k]);	
			}
		}
	}
	ret = fftn(3, dims3, rDj, iDj,  1, INVERSE_SCALE); //combination of (1, -1.0)==ifftn.m in matlab
	if (ret) exit(1);
		
	for (int k=0;k<nz;k++){	
		for (int i=0;i<nrows;i++) {
			for (int j=0;j<ncols;j++) {
				//if (fabs(vp[k*nrows*ncols+i*ncols+j]-330)<1e-9 || i==0 || i==nrows-1 || j==0 || j==ncols-1 || k==0 || k==nz-1) 
				//if (fabs(rDj[k*nrows*ncols+i*ncols+j])<1e-10)
				if (fabs(vp[k*nrows*ncols+i*ncols+j]-330)<1e-9 )
					Dj[k*nrows*ncols+i*ncols+j]=0.0;
				else
					Dj[k*nrows*ncols+i*ncols+j]=rDj[k*nrows*ncols+i*ncols+j];
			}
		}
	}
	delete [] rDj;rDj=NULL;
	delete [] iDj;iDj=NULL;
	delete [] kx;kx=NULL;
	delete [] ky;ky=NULL;
	delete [] kz;kz=NULL;	
}

double beta_HS(double *gk1,double *gk0,double *dk, int nrows, int ncols, int nz){
	double beta=0.0,tmpa=0.0,tmpb=0.0;
	
	for (int i=0;i<nrows*ncols*nz;i++){
		tmpa+=gk1[i]*(gk1[i]-gk0[i]);
		tmpb+=-dk[i]*(gk1[i]-gk0[i]);
	}
	beta=tmpa/tmpb;
	
	return beta;
}

double beta_CGDESCENT(double *gk1,double *gk0,double *dk, int nrows, int ncols, int nz){
	double beta=0.0,dTy=0.0,ny=0.0;
	
	for (int i=0;i<nrows*ncols*nz;i++){
		dTy+=dk[i]*(-gk1[i]+gk0[i]);
		ny+=(-gk1[i]+gk0[i])*(-gk1[i]+gk0[i]);
	}
	for (int i=0;i<nrows*ncols*nz;i++){
		beta+=((-gk1[i]+gk0[i])-2*dk[i]*ny/dTy)*(-gk1[i])/dTy;
	}
	return beta;
}

double gd2sk(double *gk1,double *gk0,double *dk, double alpha1, double *rhok, double *sk, double *yk, int Iterk, int lm, int nrows, int ncols, int nz){
	double gamma=0.0,dTy=0.0,ny=0.0;
	
	for (int i=0;i<nrows*ncols*nz;i++){
		dTy+=alpha1*dk[i]*(-gk1[i]+gk0[i]);
		ny+=(-gk1[i]+gk0[i])*(-gk1[i]+gk0[i]);
	}
	gamma=dTy/ny;

	rhok[(Iterk-2)%lm]=1.0/dTy;
	for (int i=0;i<nrows*ncols*nz;i++){
		sk[((Iterk-2)%lm)*nrows*ncols*nz+i]=alpha1*dk[i];
		yk[((Iterk-2)%lm)*nrows*ncols*nz+i]=-gk1[i]+gk0[i];
	}
	return gamma;
}

double zoom(double alpha0,double alpha1,double phi0,double phi1,double dphi0,double dphi1){
	double a,b,alphac;
	
	a=(2*(phi0-phi1)+(alpha1-alpha0)*(dphi0-dphi1))/((alpha1-alpha0)*(alpha1-alpha0)*(alpha1-alpha0));
	b=(3*(phi1-phi0)-(alpha1-alpha0)*(dphi1+2*dphi0))/((alpha1-alpha0)*(alpha1-alpha0));
	
	if ((b*b-3*a*dphi0)>=0) alphac=alpha0+(-b+sqrt(b*b-3*a*dphi0))/3.0/a;
	else alphac=alpha1-(alpha1-alpha0)*dphi1/(dphi1-dphi0);
	return alphac;
}
double cubicInterp(double alpha0,double alpha1,double phi0,double phi1,double dphi0,double dphi1){
	double a,b,c,d;
	a=(-2.*(phi1-phi0)+(dphi0+dphi1)*(alpha1-alpha0))/((alpha1-alpha0)*(alpha1-alpha0)*(alpha1-alpha0));
	b=(3.*(phi1-phi0)-(2.*dphi0+dphi1)*(alpha1-alpha0))/((alpha1-alpha0)*(alpha1-alpha0));
	c=dphi0;d=phi0;
	
	if (fabs(a)>0 && (b*b-3*a*c)>0)
		return alpha0+(-b+sqrt(b*b-3*a*c))/(3.*a);
	if (a==0 && fabs(b)>0 && (b*b-3*a*c)<0)
		return -c/(2.*b);
}

double magnitudeOf(double *g,int nrows,int ncols,int nz) {
	double maxg=abs(g[0]);
	for (int kk=1;kk<nrows*ncols*nz;kk++) 
		if (abs(g[kk])>maxg) maxg=abs(g[kk]);
	
	return 1./maxg;
}
double meanOf(double *g,int nrows,int ncols,int nz) {
	double meang=0;
	for (int kk=0;kk<nrows*ncols*nz;kk++) 
		meang=meang+g[kk];
	meang=meang/(nrows*ncols*nz);
	
	return 1./fabs(meang);
}

#endif
