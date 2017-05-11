#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include "inc.h"
#include "rmn.h"
#include "fld.h"
#include "cll.h"
#include "eos.h"
#include "trancoeff.h"
#include "cornelius.h"
//#include "fast.h"

#define OUTPI

using namespace std ;

size_t width = 15;

// returns the velocities in cartesian coordinates, fireball rest frame. Y=longitudinal rapidity of fluid
void Fluid::getCMFvariables(Cell *c, double tau, double &e, double &nb, double &nq, double &ns, double &vx, double &vy, double &Y)
{
	double p, vz;
	c->getPrimVar(eos, tau, e, p, nb, nq, ns, vx, vy, vz) ;
	double eta = getZ(c->getZ()) ;
	Y = eta + 0.5*log((1.+vz)/(1.-vz)) ;
	double factor = cosh(Y-eta)/cosh(Y) ;
	vx = vx * factor;
	vy = vy * factor;
}

// unput: geom. rapidity + velocities in Bjorken frame, --> output: velocities in lab.frame
void transformToLab(double eta, double &vx, double &vy, double &vz)
{
	const double Y = eta + 0.5*log((1.+vz)/(1.-vz)) ;
	double factor = cosh(Y-eta)/cosh(Y) ;
	vx = vx * factor ;
	vy = vy * factor ;
  	vz = tanh(Y) ;
}

/////**** modified by Yingru a little bit, now return not Y, but vz *******////
void Fluid::getCMFvariables_new(Cell *c, double tau, double &e, double &nb, double &nq, double &ns, double &vx, double &vy, double &vz)
{
	double p;
	c->getPrimVar(eos, tau, e, p, nb, nq, ns, vx, vy, vz) ;
	double eta = getZ(c->getZ()) ;
	transformToLab(eta, vx, vy, vz);
}


Fluid::Fluid(EoS *_eos, TransportCoeff *_trcoeff, int _nx, int _ny, int _nz, double _minx, double _maxx, double _miny, double _maxy, double _minz, double _maxz, double _dt,  double eCrit)
{
	eos = _eos ;
	trcoeff = _trcoeff ;
	nx = _nx ;
	ny = _ny ;
	nz = _nz ;
	minx = _minx ;
	maxx = _maxx ;
	miny = _miny ;
	maxy = _maxy ;
	minz = _minz ;
	maxz = _maxz ;
	dx = (maxx-minx)/(nx-1) ;
	dy = (maxy-miny)/(ny-1) ;
	dz = (maxz-minz)/(nz-1) ;
 	dt = _dt ;

	cell = new Cell [nx*ny*nz] ;

	cell0 = new Cell ;
	cell0 ->setPrimVar(eos, 1.0, 0., 0., 0., 0., 0., 0., 0.) ; // tau is not important here, since *0
	cell0 ->setAllM(0.) ;

	for(int ix=0; ix<nx; ix++)
	for(int iy=0; iy<ny; iy++)
	for(int iz=0; iz<nz; iz++)
	{
		getCell(ix,iy,iz)->setPrev(X_, getCell(ix-1, iy, iz)) ;
		getCell(ix,iy,iz)->setNext(X_, getCell(ix+1, iy, iz)) ;
		getCell(ix,iy,iz)->setPrev(Y_, getCell(ix, iy-1, iz)) ;
		getCell(ix,iy,iz)->setNext(Y_, getCell(ix, iy+1, iz)) ;
		getCell(ix,iy,iz)->setPrev(Z_, getCell(ix, iy, iz-1)) ;
		getCell(ix,iy,iz)->setNext(Z_, getCell(ix, iy, iz+1)) ;
		getCell(ix,iy,iz)->setPos(ix, iy, iz) ;
	}

	output_nt = 0 ;
	output_nx = 0 ;
	output_ny = 0 ;
  
//---- Cornelius init
  	double arrayDx [4] = {dt, dx, dy, dz} ;
  	cornelius = new Cornelius ;
 	cornelius->init(4, eCrit, arrayDx) ;
       	ecrit = eCrit ;
  	vEff = 0. ;
  	EtotSurf = 0.0 ;
}


Fluid::~Fluid()
{
	delete [] cell ;
	delete cell0 ;
}


void Fluid::initOutput(char *dir, int maxstep, double tau0, int cmpr2dOut)
{
//    directory = dir ;
	compress2dOut = cmpr2dOut ;
	cout << "maxstep=" << maxstep << endl ;
	char command [255] ;
	sprintf(command, "mkdir -p %s",dir) ;
	system(command) ;

	string outfreeze = dir; 
	outfreeze.append("/freezeout.dat");
	ffreeze.open(outfreeze.c_str()) ;
	
	// freezeout header file
	ffreeze << "# t x y z st sx sy sz vx vy vz pi11 pi12 pi13 pi22 pi23 pi33" << endl;
}


void Fluid::correctImagCells(void)
{
	double Q [7] ;
// Z
	for(int ix=0; ix<nx; ix++)
	{
		for(int iy=0; iy<ny; iy++)
		{
			// left boundary
			getCell(ix,iy,2)->getQ(Q) ;
			getCell(ix,iy,1)->setQ(Q) ;
			getCell(ix,iy,0)->setQ(Q) ;
			// right boundary
			getCell(ix,iy,nz-3)->getQ(Q) ;
			getCell(ix,iy,nz-2)->setQ(Q) ;
			getCell(ix,iy,nz-1)->setQ(Q) ;
		}
	}
// Y
	for(int ix=0; ix<nx; ix++)
	{
		for(int iz=0; iz<nz; iz++)
		{
			// left boundary
			getCell(ix,2,iz)->getQ(Q) ;
			getCell(ix,1,iz)->setQ(Q) ;
			getCell(ix,0,iz)->setQ(Q) ;
			// right boundary
			getCell(ix,ny-3,iz)->getQ(Q) ;
			getCell(ix,ny-2,iz)->setQ(Q) ;
			getCell(ix,ny-1,iz)->setQ(Q) ;
		}
	}
// X
	for(int iy=0; iy<ny; iy++)
	{
		for(int iz=0; iz<nz; iz++)
		{
			// left boundary
			getCell(2,iy,iz)->getQ(Q) ;
			getCell(1,iy,iz)->setQ(Q) ;
			getCell(0,iy,iz)->setQ(Q) ;
			// right boundary
			getCell(nx-3,iy,iz)->getQ(Q) ;
			getCell(nx-2,iy,iz)->setQ(Q) ;
			getCell(nx-1,iy,iz)->setQ(Q) ;
		}
	}
}


void Fluid::correctImagCellsFull(void)
{
	double Q [7], _pi[4][4], _Pi ;
// Z
	for(int ix=0; ix<nx; ix++)
	{
		for(int iy=0; iy<ny; iy++)
		{
			// left boundary
			getCell(ix,iy,2)->getQ(Q) ;
			getCell(ix,iy,1)->setQ(Q) ;
			getCell(ix,iy,0)->setQ(Q) ;
			for(int i=0; i<4; i++)
			{
				for(int j=0; j<4; j++)
				{
					_pi[i][j]=getCell(ix,iy,2)->getpi(i,j) ;
				}
			}
			_Pi=getCell(ix,iy,2)->getPi() ;
		
			for(int i=0; i<4; i++)
			{
				for(int j=0; j<=i; j++)
				{
					getCell(ix,iy,0)->setpi(i,j,_pi[i][j]) ;
					getCell(ix,iy,1)->setpi(i,j,_pi[i][j]) ;
				}
			}
			getCell(ix,iy,0)->setPi(_Pi) ;
			getCell(ix,iy,1)->setPi(_Pi) ;
		
			// right boundary
			getCell(ix,iy,nz-3)->getQ(Q) ;
			getCell(ix,iy,nz-2)->setQ(Q) ;
			getCell(ix,iy,nz-1)->setQ(Q) ;
			for(int i=0; i<4; i++)
			{
				for(int j=0; j<4; j++)
				{
					_pi[i][j]=getCell(ix,iy,nz-3)->getpi(i,j) ;
				}
			}
			_Pi=getCell(ix,iy,nz-3)->getPi() ;
		
			for(int i=0; i<4; i++)
			{
				for(int j=0; j<=i; j++)
				{
					getCell(ix,iy,nz-2)->setpi(i,j,_pi[i][j]) ;
					getCell(ix,iy,nz-1)->setpi(i,j,_pi[i][j]) ;
				}
			}
			getCell(ix,iy,nz-2)->setPi(_Pi) ;
			getCell(ix,iy,nz-1)->setPi(_Pi) ;
		}
	}
// Y
	for(int ix=0; ix<nx; ix++)
	{	for(int iz=0; iz<nz; iz++)
		{
			// left boundary
			getCell(ix,2,iz)->getQ(Q) ;
			getCell(ix,1,iz)->setQ(Q) ;
			getCell(ix,0,iz)->setQ(Q) ;
			for(int i=0; i<4; i++)
			{
				for(int j=0; j<4; j++)
				{
					_pi[i][j]=getCell(ix,2,iz)->getpi(i,j) ;
				}
			}
			_Pi=getCell(ix,2,iz)->getPi() ;
		
			for(int i=0; i<4; i++)
			{
				for(int j=0; j<=i; j++)
				{
					getCell(ix,0,iz)->setpi(i,j,_pi[i][j]) ;
					getCell(ix,1,iz)->setpi(i,j,_pi[i][j]) ;
				}
			}
			getCell(ix,0,iz)->setPi(_Pi) ;
			getCell(ix,1,iz)->setPi(_Pi) ;
		
			// right boundary
			getCell(ix,ny-3,iz)->getQ(Q) ;
			getCell(ix,ny-2,iz)->setQ(Q) ;
			getCell(ix,ny-1,iz)->setQ(Q) ;
			for(int i=0; i<4; i++)
			{
				for(int j=0; j<4; j++)
				{
					_pi[i][j]=getCell(ix,ny-3,iz)->getpi(i,j) ;
				}
			}
			_Pi=getCell(ix,ny-3,iz)->getPi() ;
		
			for(int i=0; i<4; i++)
			{
				for(int j=0; j<=i; j++)
				{
					getCell(ix,ny-2,iz)->setpi(i,j,_pi[i][j]) ;
					getCell(ix,ny-1,iz)->setpi(i,j,_pi[i][j]) ;
				}
			}
			getCell(ix,ny-2,iz)->setPi(_Pi) ;
			getCell(ix,ny-1,iz)->setPi(_Pi) ;
		}
	}
// X
	for(int iy=0; iy<ny; iy++)
	{
		for(int iz=0; iz<nz; iz++)
		{
			// left boundary
			getCell(2,iy,iz)->getQ(Q) ;
			getCell(1,iy,iz)->setQ(Q) ;
			getCell(0,iy,iz)->setQ(Q) ;
			for(int i=0; i<4; i++)
				for(int j=0; j<4; j++)
					_pi[i][j]=getCell(2,iy,iz)->getpi(i,j) ;
			_Pi=getCell(2,iy,iz)->getPi() ;
		
			for(int i=0; i<4; i++)
			{
				for(int j=0; j<=i; j++)
				{
					getCell(0,iy,iz)->setpi(i,j,_pi[i][j]) ;
					getCell(1,iy,iz)->setpi(i,j,_pi[i][j]) ;
				}
			}
			getCell(0,iy,iz)->setPi(_Pi) ;
			getCell(1,iy,iz)->setPi(_Pi) ;
			// right boundary
			getCell(nx-3,iy,iz)->getQ(Q) ;
			getCell(nx-2,iy,iz)->setQ(Q) ;
			getCell(nx-1,iy,iz)->setQ(Q) ;
			for(int i=0; i<4; i++)
				for(int j=0; j<4; j++)
					_pi[i][j]=getCell(nx-3,iy,iz)->getpi(i,j) ;
			_Pi=getCell(nx-3,iy,iz)->getPi() ;
		
			for(int i=0; i<4; i++)
			{	
				for(int j=0; j<=i; j++)
				{
					getCell(nx-2,iy,iz)->setpi(i,j,_pi[i][j]) ;
					getCell(nx-1,iy,iz)->setpi(i,j,_pi[i][j]) ;
				}
			}
			getCell(nx-2,iy,iz)->setPi(_Pi) ;
			getCell(nx-1,iy,iz)->setPi(_Pi) ;
		}
	}
}
    

void Fluid::updateM(double tau, double dt)
{
	for(int ix=0; ix<getNX(); ix++)
	{
		for(int iy=0; iy<getNY(); iy++)
		{
			for(int iz=0; iz<getNZ(); iz++)
			{
				Cell* c = getCell(ix,iy,iz) ;
				c->setDM(X_, 0.) ;
				c->setDM(Y_, 0.) ;
				c->setDM(Z_, 0.) ;
				if(getCell(ix,iy,iz)->getLM()<1.)
				{
					if(getCell(ix+1,iy,iz)->getM(X_)>=1. || getCell(ix-1,iy,iz)->getM(X_)>=1.) c->setDM(X_, dt/dx) ;
					if(getCell(ix,iy+1,iz)->getM(Y_)>=1. || getCell(ix,iy-1,iz)->getM(Y_)>=1.) c->setDM(Y_, dt/dy) ;
					if(getCell(ix,iy,iz+1)->getM(Z_)>=1. || getCell(ix,iy,iz-1)->getM(Z_)>=1.) c->setDM(Z_, dt/dz/tau) ;

					if(c->getDM(X_)==0. && c->getDM(Y_)==0.)
					{
						if(getCell(ix+1,iy+1,iz)->getLM()>=1. || getCell(ix+1,iy-1,iz)->getLM()>=1. || getCell(ix-1,iy+1,iz)->getLM()>=1. || getCell(ix-1,iy-1,iz)->getLM()>=1.)
						{
							c->setDM(X_, 0.707*dt/dx) ;
							c->setDM(Y_, 0.707*dt/dy) ;
						}
					}
				} //if
			}
		}
	}

	for(int ix=0; ix<getNX(); ix++)
	{
		for(int iy=0; iy<getNY(); iy++)
		{
			for(int iz=0; iz<getNZ(); iz++)
			{
				Cell* c = getCell(ix,iy,iz) ;
				c->addM(X_,c->getDM(X_)) ;
				c->addM(Y_,c->getDM(Y_)) ;
				c->addM(Z_,c->getDM(Z_)) ;
			}
		}
	}
}


// Modified by Weiyao to reduce the hyper surface output
void Fluid::outputSurface(double tau)
{
	static double nbSurf = 0.0 ;
	double e, p, nb, nq, ns, t, mub, muq, mus, vx, vy, vz, Q[7] ;
	double 	E=0., Efull=0., S=0., Px=0., 
			vt_num=0., vt_den=0., vxvy_num=0.,vxvy_den=0., 
			pi0x_num=0., pi0x_den=0., txxyy_num=0., txxyy_den=0., Nb1 = 0., Nb2 = 0. ;
	double eta=0 ;
	int nelements=0, nsusp=0 ; // all surface emenents and suspicious ones
	int nCoreCells=0, nCoreCutCells=0 ; // cells with e>eCrit and cells with cut visc.corr.
	//-- Cornelius: allocating memory for corner points
	double ****ccube = new double***[2];
  	for (int i1=0; i1 < 2; i1++)
  	{
    	ccube[i1] = new double**[2];
    	for (int i2=0; i2 < 2; i2++)
    	{
      		ccube[i1][i2] = new double*[2];
      		for (int i3=0; i3 < 2; i3++)
      		{
        		ccube[i1][i2][i3] = new double[2];
      		}
    	}
	}
	//----end Cornelius
	
	double cosh_int, sinh_int, tanh_vz, tanh_vz_sq2, norm_v, product_v,
	  vx2, vy2, vt_sq2, norm_vt, eh, cosh_eta, sinh_eta, norm_v_sq2;
		
	for(int ix=2; ix<nx-2; ix++)
	{
		for(int iy=2; iy<ny-2; iy++)
		{
	 		for(int iz=2; iz<nz-2; iz++)
	 		{
  				Cell *c = getCell(ix,iy,iz) ;
  				getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz) ;
  				c->getQ(Q) ;
  				eos->eos(e, nb, nq, ns, t, mub, muq, mus, p);
  				double s = eos->s(e, nb, nq, ns);
  				eta=getZ(iz) ;
     				// >> aux-variable
				cosh_int = (sinh(eta+0.5*dz)-sinh(eta-0.5*dz))/dz ;
  				sinh_int = (cosh(eta+0.5*dz)-cosh(eta-0.5*dz))/dz ;
								
				tanh_vz = tanh(vz);
				tanh_vz_sq2 = tanh_vz*tanh_vz;
				norm_v = sqrt(1. - vx*vx - vy*vy - tanh_vz_sq2);
				product_v = cosh_int - tanh_vz*sinh_int;
				norm_v_sq2 = norm_v*norm_v;
  				// << aux-variable

				eh = e + p;
				E += tau * (eh/norm_v_sq2*product_v - p*cosh_int) ;
 		 		Nb1 += Q[NB_] ;
  				Nb2 += tau*nb*product_v/norm_v;
//---- inf check
  				if(isinf(E))
  				{
				    cout<<"EEinf"<<setw(14)<<e<<setw(14)<<p<<setw(14)<<vx<<setw(14)<<vy<<setw(14)<<vz<<endl ;
				    exit(1) ;
  				}
				//-------------
				// >> aux-variable
			        sinh_eta = sinh(eta);
			        cosh_eta = sqrt(1.+sinh_eta*sinh_eta);
				// << aux-variable
  				Efull += tau * (eh/norm_v_sq2*(cosh_eta-tanh_vz*sinh_eta) - p*cosh_eta);

  				if(trcoeff->isViscous()) Efull += tau*(c->getpi(0,0)*cosh_eta + c->getpi(0,3)*sinh_eta);
  				
				if(e>ecrit)
  				{
   					nCoreCells++ ;
   					if(c->getViscCorrCutFlag()<0.9) nCoreCutCells++ ;
  				}
  				// -- noneq. corrections to entropy flux
  				const double gmumu[4] = {1., -1., -1., -1.} ;
  				double deltas = 0. ;
  				if(trcoeff->isViscous())
  				{
  					for(int i=0; i<4; i++)
  					{
  						for(int j=0; j<4; j++)
  						{
  							deltas += pow(c->getpi(i,j),2)*gmumu[i]*gmumu[j] ;
  						}
  					}
  				}
  				if(t>0.02)
  				{
   					s += 1.5*deltas/(eh*t) ;
   					S += tau*s*product_v/norm_v ;
  				}
				// >> aux-variable
				vx2 = vx*vx;
				vy2 = vy*vy;
				vt_sq2 = vx2 + vy2;
				norm_vt = sqrt(1. - vt_sq2);
				// << aux-variable
  				Px += tau*eh*vx/norm_v ;
  				vt_num += e/norm_vt*sqrt(vt_sq2) ;
  				vt_den += e/norm_vt ;
  				vxvy_num += e*(fabs(vx)-fabs(vy)) ;
  				vxvy_den += e ;
  				txxyy_num += eh/norm_v_sq2*(vx2-vy2) ;
  				txxyy_den += eh/norm_v_sq2*(vt_sq2)+2.*p ;
  				pi0x_num += e/norm_v_sq2*fabs(c->getpi(0,1)) ;
  				pi0x_den += e/norm_v_sq2;

//----- Cornelius stuff
  				double QCube[2][2][2][2][7] ;
 				double piSquare[2][2][2][10], PiSquare[2][2][2] ;
  				for(int jx=0; jx<2; jx++)
  				{
  					for(int jy=0; jy<2; jy++)
  					{
  						for(int jz=0; jz<2; jz++)
  						{
						      double _p, _nb, _nq, _ns, _vx, _vy, _vz ;
						      Cell* cc = getCell(ix+jx,iy+jy,iz+jz) ;
						      cc->getPrimVar(eos,tau,e,_p,_nb,_nq,_ns,_vx,_vy,_vz) ;
						      
						      cc->getQ(QCube[1][jx][jy][jz]) ;
						      ccube[1][jx][jy][jz] = e ;

						      cc->getPrimVarPrev(eos,tau-dt,e,_p,_nb,_nq,_ns,_vx,_vy,_vz) ;
						      cc->getQprev(QCube[0][jx][jy][jz]) ;
						      ccube[0][jx][jy][jz] = e ;
    						// ---- get viscous tensor
						      for(int ii=0; ii<4; ii++)
	  					      {
							  for(int jj=0; jj<=ii; jj++)
							  {
							    piSquare[jx][jy][jz][index44(ii,jj)] = cc->getpi(ii,jj) ;
							  }
						      }		
						      PiSquare[jx][jy][jz] = cc->getPi() ;
						}
					}
				}

  				cornelius->find_surface_4d(ccube);
 				const int Nsegm = cornelius->get_Nelements() ;
  				for(int isegm=0; isegm<Nsegm; isegm++)
  				{
				    nelements++ ;
				    ffreeze.precision(6);
				    
					// Weiyao: from curvelinear cx^mu to cartesian x^mu
					double cx0 = tau+cornelius->get_centroid_elem(isegm,0), 
						   x1 = getX(ix)+cornelius->get_centroid_elem(isegm,1), 						   x2 = getY(iy)+cornelius->get_centroid_elem(isegm,2),
						   cx3 = getZ(iz)+cornelius->get_centroid_elem(isegm,3);
					double x0 = cx0*cosh(cx3), x3 = cx0*sinh(cx3);
				    ffreeze << setw(width) << x0 
							<< setw(width) << x1
					   		<< setw(width) << x2
							<< setw(width) << x3;
				    
				    // ---- interpolation procedure
				    double vxC=0., vyC=0., vzC=0., TC=0., mubC=0., muqC=0., musC=0., piC[10], PiC=0., nbC=0., nqC=0. ; // values at the centre, to be interpolated
				    double QC [7] = {0., 0., 0., 0., 0., 0., 0.} ;
				    double eC=0., pC=0. ;
				    for(int ii=0; ii<10; ii++) piC[ii] = 0.0 ;
				    double delta_x0 = cornelius->get_centroid_elem(isegm,0)/dt;
				    double delta_x1 = cornelius->get_centroid_elem(isegm,1)/dx;
				    double delta_x2 = cornelius->get_centroid_elem(isegm,2)/dy;
				    double delta_x3 = cornelius->get_centroid_elem(isegm,3)/dz;
				    double wCenT[2] = {1. - delta_x0, delta_x0} ;
				    double wCenX[2] = {1. - delta_x1, delta_x1} ;
				    double wCenY[2] = {1. - delta_x2, delta_x2} ;
				    double wCenZ[2] = {1. - delta_x3, delta_x3} ;
				    for(int jt=0; jt<2; jt++)
    					for(int jx=0; jx<2; jx++)
    						for(int jy=0; jy<2; jy++)
    							for(int jz=0; jz<2; jz++)
    								for(int i=0; i<7; i++)
    								{
    	  								QC[i] += QCube[jt][jx][jy][jz][i]*wCenT[jt]*wCenX[jx]*wCenY[jy]*wCenZ[jz] ;
    								}
   				    for(int i=0; i<7; i++) QC[i] = QC[i]/(tau+cornelius->get_centroid_elem(isegm,0)) ;
				    double _ns = 0.0 ;
				    transformPV(eos, QC, eC, pC, nbC, nqC, _ns, vxC, vyC, vzC) ;
				    eos->eos(eC,nbC,nqC,_ns,TC,mubC,muqC,musC,pC) ;
				    if(TC>0.4 || fabs(mubC)>0.85)
				    {
    	 				cout << "#### Error (surface): high T/mu_b ####\n" ;
				    }
				    if(eC>ecrit*2.0 || eC<ecrit*0.5) nsusp++ ;
				    
				    double dV_tmp;
				    for(int jx=0; jx<2; jx++)
				      for(int jy=0; jy<2; jy++)
				        for(int jz=0; jz<2; jz++)
					  {
					    dV_tmp = wCenX[jx]*wCenY[jy]*wCenZ[jz];
					    for(int ii=0; ii<10; ii++)
					      piC[ii] += piSquare[jx][jy][jz][ii]*dV_tmp ;
					    PiC += PiSquare[jx][jy][jz]*dV_tmp ;
					  }
				    double v2C = vxC*vxC+vyC*vyC+vzC*vzC ;
				    if(v2C>1.)
				      {
					double v2C_factor = sqrt(0.99/v2C);
					vxC *= v2C_factor ;
					vyC *= v2C_factor ;
					vzC *= v2C_factor ;
					v2C = 0.99 ;
				      }
				    double etaC = getZ(iz) + cornelius->get_centroid_elem(isegm,3) ;
				    transformToLab(etaC, vxC, vyC, vzC) ; // viC is now in lab.frame!
				    double gammaC = 1./sqrt(1.-vxC*vxC-vyC*vyC-vzC*vzC) ;
				    double uC[4] = {gammaC, gammaC*vxC, gammaC*vyC, gammaC*vzC} ;
				    const double tauC = tau+cornelius->get_centroid_elem(isegm,0) ;
				    double dsigma[4] ;
				    // ---- transform dsigma to lab.frame :
				    const double ch = cosh(etaC) ;
				    const double sh = sinh(etaC) ;
				    dsigma [0] = tauC*( ch*cornelius->get_normal_elem(0,0) - sh/tauC*cornelius->get_normal_elem(0,3) ) ;
				    dsigma [3] = tauC*(-sh*cornelius->get_normal_elem(0,0) + ch/tauC*cornelius->get_normal_elem(0,3) );
				    dsigma [1] = tauC*cornelius->get_normal_elem(0,1) ;
				    dsigma [2] = tauC*cornelius->get_normal_elem(0,2) ;
				    double dVEff = 0.0 ;
				    for(int ii=0; ii<4; ii++) dVEff += dsigma[ii]*uC[ii] ; // normalize for Delta eta=1
				    vEff += dVEff ;

					// Weiyao: output dsigma_mu
				    for(int ii=0; ii<4; ii++) ffreeze << setw(width) << dsigma[ii] ;
				    // Weiyao: output 3-velocity
					ffreeze << setw(width) << vxC
				    		<< setw(width) << vyC
							<< setw(width) << vzC;
#ifdef OUTPI
					// Weiyao: only output pi tensor spatial component 11, 12, 13, 22, 23, other components can be reconstructed by symmetry+traceless+orthogonality
				    double picart[5];
				    /*pi11*/ picart[0] = piC[index44(1,1)] ;
				    /*pi12*/ picart[1] = piC[index44(1,2)] ;
				    /*pi13*/ picart[2] = sh*piC[index44(0,1)]+ch*piC[index44(3,1)] ;
				    /*pi22*/ picart[3] = piC[index44(2,2)] ;
				    /*pi23*/ picart[4] = sh*piC[index44(0,2)]+ch*piC[index44(3,2)] ;
				    for(int ii=0; ii<5; ii++) ffreeze <<setw(width) << picart[ii] ;

					// no bulk output!
				    ffreeze << endl;
#else	
				    ffreeze << setw(width) << dVEff << endl ;
#endif	
				    double dEsurfVisc = 0. ;
				    for(int i=0; i<4; i++) dEsurfVisc += picart[index44(0,i)]*dsigma[i] ;
				    EtotSurf += (eC+pC)*uC[0]*dVEff - pC*dsigma[0] + dEsurfVisc ;
				    nbSurf += nbC*dVEff ;
				    
  				}
  	  			//----- end Cornelius
			}
		}
 	}
	double dV = dx*dy*dz;
 	E=E*dV ;
	Efull=Efull*dV ;
 	S=S*dV ;
 	Nb1 *= dV ; 
	Nb2 *= dV ;
 	cout << setw(10) << tau << setw(13) << E << setw(13) << Efull << setw(13) << nbSurf
	     << setw(13) << S << setw(10) << nelements << setw(10) << nsusp 
	     << setw(13) << (float)(nCoreCutCells)/(float)(nCoreCells) << endl ;
	//-- Cornelius: all done, let's free memory
  	for (int i1=0; i1 < 2; i1++)
  	{
	    for (int i2=0; i2 < 2; i2++)
	    {
    	  	for (int i3=0; i3 < 2; i3++)
    	  	{
    	  		delete[] ccube[i1][i2][i3];
    	  	}
    	  	delete[] ccube[i1][i2];
	    }
	    delete[] ccube[i1];
  	}
  	delete[] ccube;
	//----end Cornelius	
  	if(nelements==0) exit(0) ;
}


void Fluid::outputCorona(double tau)
{
    static double nbSurf = 0.0 ;
    double e, p, nb, nq, ns, t, mub, muq, mus, vx, vy, vz, Q[7] ;
    double E=0., Efull=0., S=0., Px=0., vt_num=0., vt_den=0., vxvy_num=0.,
    vxvy_den=0., pi0x_num=0., pi0x_den=0., txxyy_num=0., txxyy_den=0., Nb1 = 0., Nb2 = 0. ;
    double eta=0 ;
    int nelements=0 ;

    for(int ix=2; ix<nx-2; ix++)
        for(int iy=2; iy<ny-2; iy++)
	    for(int iz=2; iz<nz-2; iz++)
	    {
	        Cell *c = getCell(ix,iy,iz) ;
		getCMFvariables(c, tau, e, nb, nq, ns, vx, vy, vz) ;
		c->getQ(Q) ;
		eos->eos(e, nb, nq, ns, t, mub, muq, mus, p);
		double s = eos->s(e, nb, nq, ns) ;
		eta=getZ(iz) ;
		const double cosh_int = (sinh(eta+0.5*dz)-sinh(eta-0.5*dz))/dz ;
		const double sinh_int = (cosh(eta+0.5*dz)-cosh(eta-0.5*dz))/dz ;
		
		E += tau*(e+p)/(1.-vx*vx-vy*vy-tanh(vz)*tanh(vz))*(cosh_int-tanh(vz)*sinh_int) - tau*p*cosh_int ;
		Nb1 += Q[NB_] ;
		Nb2 += tau*nb*(cosh_int-tanh(vz)*sinh_int)/sqrt(1.-vx*vx-vy*vy-tanh(vz)*tanh(vz)) ;
//---- inf check
		if(isinf(E))
		{
		    cout<<"EEinf"<<setw(14)<<e<<setw(14)<<p<<setw(14)<<vx<<setw(14)<<vy<<setw(14)<<vz<<endl ;
		    exit(1) ;
		}
//--------------
		Efull += tau*(e+p)/(1.-vx*vx-vy*vy-tanh(vz)*tanh(vz))*(cosh(eta)-tanh(vz)*sinh(eta)) - tau*p*cosh(eta) ;
		if(trcoeff->isViscous()) Efull += tau*c->getpi(0,0)*cosh(eta)+tau*c->getpi(0,3)*sinh(eta);
		// -- noneq. corrections to entropy flux
		const double gmumu[4] = {1., -1., -1., -1.} ;
		double deltas = 0. ;
		if(trcoeff->isViscous())
		    for(int i=0; i<4; i++)
		        for(int j=0; j<4; j++)
			    deltas += pow(c->getpi(i,j),2)*gmumu[i]*gmumu[j] ;
		if(t>0.02)
		{
		    s += 1.5*deltas/((e+p)*t) ;
		    S += tau*s*(cosh_int-tanh(vz)*sinh_int)/sqrt(1.-vx*vx-vy*vy-tanh(vz)*tanh(vz)) ;
		}
		Px += tau*(e+p)*vx/(1.-vx*vx-vy*vy-tanh(vz)*tanh(vz)) ;
		vt_num += e/sqrt(1.-vx*vx-vy*vy)*sqrt(vx*vx+vy*vy) ;
		vt_den += e/sqrt(1.-vx*vx-vy*vy) ;
		vxvy_num += e*(fabs(vx)-fabs(vy)) ;
		vxvy_den += e ;
		txxyy_num += (e+p)/(1.-vx*vx-vy*vy-tanh(vz)*tanh(vz))*(vx*vx-vy*vy) ;
		txxyy_den += (e+p)/(1.-vx*vx-vy*vy-tanh(vz)*tanh(vz))*(vx*vx+vy*vy)+2.*p ;
		pi0x_num += e/(1.-vx*vx-vy*vy-tanh(vz)*tanh(vz))*fabs(c->getpi(0,1)) ;
		pi0x_den += e/(1.-vx*vx-vy*vy-tanh(vz)*tanh(vz)) ;

//----- Cornelius stuff
		bool isCorona = true, isTail = true ;
		double QCube[2][2][2][7] ;
		double piSquare[2][2][2][10], PiSquare[2][2][2] ;
		for(int jx=0; jx<2; jx++)
		    for(int jy=0; jy<2; jy++)
		        for(int jz=0; jz<2; jz++)
			{
			    double _p, _nb, _nq, _ns, _vx, _vy, _vz ;
			    Cell* cc = getCell(ix+jx,iy+jy,iz+jz) ;
			    cc->getPrimVar(eos,tau,e,_p,_nb,_nq,_ns,_vx,_vy,_vz) ;
			    cc->getQ(QCube[jx][jy][jz]) ;
			    if(e>ecrit) isCorona = false ;
			    if(e>0.01) isTail = false ;
			    // ---- get viscous tensor
			    for(int ii=0; ii<4; ii++)
			        for(int jj=0; jj<=ii; jj++)
				    piSquare[jx][jy][jz][index44(ii,jj)] = cc->getpi(ii,jj) ;
			    PiSquare[jx][jy][jz] = cc->getPi() ;
			}
		if(isCorona && !isTail)
		{
		    nelements++ ;
		    ffreeze.precision(6) ;
			double cx0 = tau, 
				   x1 = getX(ix)+0.5*dx,
				   x2 = getY(iy)+0.5*dy,
				   cx3 = getZ(iz)+0.5*dz;
			double x0 = cx0*cosh(cx3), x3 = cx0*sinh(cx3);
		    ffreeze << setw(width) << x0
					<< setw(width) << x1
					<< setw(width) << x2
					<< setw(width) << x3;
		    // ---- interpolation procedure
		    double vxC=0., vyC=0., vzC=0., TC=0., mubC=0., muqC=0., musC=0., piC[10], PiC=0., nbC=0., nqC=0. ; // values at the centre, to be interpolated
		    double QC [7] = {0., 0., 0., 0., 0., 0., 0.} ;
		    double eC=0., pC=0. ;
		    for(int ii=0; ii<10; ii++) piC[ii] = 0.0 ;
		    for(int jx=0; jx<2; jx++)
		        for(int jy=0; jy<2; jy++)
			    for(int jz=0; jz<2; jz++)
			        for(int i=0; i<7; i++)
				{
				    QC[i] += QCube[jx][jy][jz][i]*0.125 ;
				}
		    for(int i=0; i<7; i++) QC[i] = QC[i]/tau ;
		    double _ns = 0.0 ;
		    transformPV(eos, QC, eC, pC, nbC, nqC, _ns, vxC, vyC, vzC) ;
		    eos->eos(eC,nbC,nqC,_ns,TC,mubC,muqC,musC,pC) ;
		    if(TC>0.4 || fabs(mubC)>0.85)
		    {
		        cout << "#### Error (surface): high T/mu_b ####\n" ;
		    }
		    for(int jx=0; jx<2; jx++)
		        for(int jy=0; jy<2; jy++)
			    for(int jz=0; jz<2; jz++)
			    {
			        for(int ii=0; ii<10; ii++)
				    piC[ii] += piSquare[jx][jy][jz][ii]*0.125 ;
				PiC += PiSquare[jx][jy][jz]*0.125 ;
			    }
		    double v2C = vxC*vxC+vyC*vyC+vzC*vzC ;
		    if(v2C>1.)
		    {
		        vxC *= sqrt(0.99/v2C) ;
		        vyC *= sqrt(0.99/v2C) ;
			vzC *= sqrt(0.99/v2C) ;
			v2C = 0.99 ;
		    }
		    double etaC = getZ(iz) + 0.5*dz ;
		    transformToLab(etaC, vxC, vyC, vzC) ; // viC is now in lab.frame!
		    double gammaC = 1./sqrt(1.-vxC*vxC-vyC*vyC-vzC*vzC) ;
		    double uC [4] = {gammaC, gammaC*vxC, gammaC*vyC, gammaC*vzC} ;
		    const double tauC = tau ;
		    double dsigma [4] ;
		    // ---- transform dsigma to lab.frame :
		    const double ch = cosh(etaC) ;
		    const double sh = sinh(etaC) ;
		    dsigma [0] = tauC*( ch*dx*dy*dz ) ;
		    dsigma [3] = tauC*(-sh*dx*dy*dz );
		    dsigma [1] = 0.0 ;
		    dsigma [2] = 0.0 ;
		    double dVEff = 0.0 ;
		    for(int ii=0; ii<4; ii++) dVEff += dsigma[ii]*uC[ii] ; // normalize for Delta eta=1
		    vEff += dVEff ;
		    for(int ii=0; ii<4; ii++) ffreeze << setw(width) << dsigma[ii] ;
		    ffreeze << setw(width) << vxC 
				    << setw(width) << vyC 
					<< setw(width) << vzC;
#ifdef OUTPI
		    double picart[5] ;
		    /*pi11*/ picart[0] = piC[index44(1,1)] ;
		    /*pi12*/ picart[1] = piC[index44(1,2)] ;
		    /*pi13*/ picart[2] = sh*piC[index44(0,1)]+ch*piC[index44(3,1)] ;
		    /*pi22*/ picart[3] = piC[index44(2,2)] ;
		    /*pi23*/ picart[4] = sh*piC[index44(0,2)]+ch*piC[index44(3,2)] ;
		    for(int ii=0; ii<5; ii++) ffreeze << setw(width) << picart[ii] ;
		    ffreeze <<endl ;
#else
		    ffreeze <<setw(width)<<dVEff<<endl ;
#endif
		    double dEsurfVisc = 0. ;
		    for(int i=0; i<4; i++) dEsurfVisc += picart[index44(0,i)]*dsigma[i] ;
		    EtotSurf += (eC+pC)*uC[0]*dVEff - pC*dsigma[0] + dEsurfVisc ;
		    nbSurf += nbC*dVEff ;
		}
//----- end Cornelius
	    }
    E=E*dx*dy*dz ;
    Efull=Efull*dx*dy*dz ;
    S=S*dx*dy*dz ;
    Nb1 *= dx*dy*dz ; Nb2 *= dx*dy*dz ;
    cout << setw(10) << "tau" << setw(13) << "E" << setw(13) << "Efull" << setw(13) << "Nb"
	 << setw(13) << "Sfull" << setw(10) << "elements" << setw(10) << "susp." << setw(13) << "\%cut" << endl ;
    cout << setw(10) << tau << setw(13) << E << setw(13) << Efull << setw(13) << nbSurf
	 << setw(13) << S << endl ;
    cout << "corona elements : " << nelements << endl ;
}
