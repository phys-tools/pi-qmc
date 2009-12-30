// $Id$
/*  Copyright (C) 2004-2006 John B. Shumway, Jr.

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.  */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "ConductivityEstimator2D.h"
#include "SimulationInfo.h"
#include "Action.h"
#include "DoubleAction.h"
#include "SuperCell.h"
#include <blitz/tinyvec.h>
#include "stats/MPIManager.h"

ConductivityEstimator2D::ConductivityEstimator2D(const SimulationInfo&simInfo,
  const double xmin, const double xmax, const double ymin, const double ymax, const int nfreq, 
  const int nxbin, const int nybin, const int nxdbin, const int nydbin, const int nstride,
  MPIManager *mpi)
  : BlitzArrayBlkdEst<7>("conductivity2D",IVecN(2*nxdbin-1,nxbin,2*nydbin-1,nybin,2,2,nfreq),true),
    npart(simInfo.getNPart()), nslice(simInfo.getNSlice()), nfreq(nfreq), nstride(nstride),
    tau(simInfo.getTau()), tauinv(1./tau), massinv(1./simInfo.getSpecies(0).mass), q(npart),
    nxbin(nxbin), nybin(nybin), nxdbin(nxdbin), nydbin(nydbin),
    xmin(xmin), xmax(xmax), ymin(ymin), ymax(ymax),
    dx((xmax-xmin)/nxbin), dy((ymax-ymin)/nybin), dxinv(1./dx), dyinv(1./dy),
    temp(2*nxdbin-1,nxbin,2*nydbin-1,nybin,2,2,nslice/nstride), jx(2,nxbin,nybin,nslice/nstride),
    jy(2,nxbin,nybin,nslice/nstride),  mpi(mpi) {
//  nbin[0]=nxbin; nbin[1]=nybin; 
//  ndbin[0]=nxdbin; ndbin[1]=nydbin;
//  min[0]=xmin; min[1]=ymin;
//  max[0]=xmax; max[1]=xmax;
//  dx[0]=(xmax-xmin)/nxbin; dx[1]=(ymax-ymin)/nybin;
//  dxinv[0]=1./(dx[0]); dxinv[1]=1./(dx[1]);
  fftw_complex *ptrx = (fftw_complex*)jx.data();
  fftw_complex *ptry = (fftw_complex*)jy.data();
  int nsliceEff=nslice/nstride;
  fwdx = fftw_plan_many_dft(1, &nsliceEff, 2*nxbin*nybin, 
                            ptrx, 0, 1, nsliceEff,
                            ptrx, 0, 1, nsliceEff,
                            FFTW_FORWARD, FFTW_MEASURE);
  fwdy = fftw_plan_many_dft(1, &nsliceEff, 2*nxbin*nybin,
                            ptry, 0, 1, nsliceEff,
                            ptry, 0, 1, nsliceEff,
                            FFTW_FORWARD, FFTW_MEASURE);
  for (int i=0; i<npart; ++i) q(i)=simInfo.getPartSpecies(i).charge;
}

ConductivityEstimator2D::~ConductivityEstimator2D() {
  fftw_destroy_plan(fwdx);
  fftw_destroy_plan(fwdy);
}

void ConductivityEstimator2D::initCalc(const int lnslice, const int firstSlice) {
  temp=0;
  jx=0; jy=0;
}

void ConductivityEstimator2D::handleLink(const Vec& start, const Vec& end, 
                                    const int ipart, const int islice, const Paths& paths) {
  if ((start[0]<xmin & end[0]<xmin) | (start[0]>xmax & end[0]>xmax)) return;
  if ((start[1]<ymin & end[1]<ymin) | (start[1]>ymax & end[1]>ymax)) return;
  int ibin = ((int)((start[0] - xmin)*dxinv));
  if (ibin < 0) {ibin = 0;}
  else if (ibin > nxbin-1) {ibin = nxbin;}
  int jbin = ((int)((end[0] - xmin)*dxinv));
  if (jbin < 0) {jbin = 0;}
  else if (jbin > nxbin-1) {jbin = nxbin;}
  if (ibin!=jbin) {
    double k = (end[1] - start[1])/(end[0] - start[0]);
    int idir = (end[0]-start[0]>0)?1:(-1);
    if (idir < 0) {ibin-=1; jbin-=1;}
    for (int i=ibin; i*idir<jbin*idir; i+=idir) {
      double x = xmin + (i+1)*dx;
      // The y index in matrix here is not the reverse of y coordinate.
      int ybin = (int)((start[1] + k*(x - start[0]) - ymin)*dyinv);
      jx(0,i,ybin,islice/nstride)+=idir*q(ipart);
    }
  }
  ibin = ((int)((start[1] - ymin)*dyinv));
  if (ibin < 0) {ibin = 0;}
  else if (ibin > nybin-1) {ibin = nybin;}
  jbin = ((int)((end[1] - ymin)*dyinv));
  if (jbin < 0) {jbin = 0;}
  else if (jbin > nybin-1) {jbin = nybin;}
  if (ibin!=jbin) {
    double kinv = (end[0] - start[0])/(end[1] - start[1]);
    int idir = (end[1]-start[1]>0)?1:(-1);
    if (idir < 0) {ibin-=1; jbin-=1;}
    for (int i=ibin; i*idir<jbin*idir; i+=idir) {
      double y = ymin + (i+1)*dy;
      int xbin = (int)((start[0] + kinv*(y - start[1]) - xmin)*dxinv);
      jy(0,xbin,i,islice/nstride)+=idir*q(ipart);
    }
  }
} 

void ConductivityEstimator2D::endCalc(const int lnslice) {
  blitz::Range allSlice = blitz::Range::all();
  blitz::Range allBin = blitz::Range::all();
  blitz::Range alldir = blitz::Range::all();
  // First move all data to 1st worker.
  int workerID = (mpi)?mpi->getWorkerID():0;
#ifdef ENABLE_MPI
  if (mpi) {
    mpi->getWorkerComm().Reduce(&jx(0,0,0,0), &jx(1,0,0,0), 2*nxbin*nybin*nslice/nstride,
                                MPI::DOUBLE, MPI::SUM, 0);
    mpi->getWorkerComm().Reduce(&jy(0,0,0,0), &jy(1,0,0,0), 2*nxbin*nybin*nslice/nstride,
                                MPI::DOUBLE, MPI::SUM, 0);
    jx(0,allBin,allBin,allSlice) = jx(1,allBin,allBin,allSlice);
    jy(0,allBin,allBin,allSlice) = jy(1,allBin,allBin,allSlice);
  }
#endif
  // Calculate correlation function using FFT's.
  if (workerID==0) {
    fftw_execute(fwdx);
    fftw_execute(fwdy); 
    for (int ibinx=0; ibinx<nxbin; ++ibinx) {
      for (int jdbinx=1; jdbinx<nxdbin; ++jdbinx) {
        int jbinx = (ibinx + jdbinx) % nxbin;
        for(int ibiny=0; ibiny<nybin; ++ibiny) {
          for(int jdbiny=1; jdbiny<nydbin; ++jdbiny) {
            // jx-jx
            int jbiny = (ibiny + jdbiny) % nybin;
            temp(jdbinx, ibinx, jdbiny, ibiny, 0, 0, allSlice)
             = conj(jx(0,ibinx,ibiny,allSlice)) * jx(0,jbinx,jbiny,allSlice);
            jbiny = (ibiny - jdbiny + nybin) % nybin;
            temp(jdbinx, ibinx, 2*nydbin-1-jdbiny, ibiny, 0, 0, allSlice)
             = conj(jx(0,ibinx,ibiny,allSlice)) * jx(0,jbinx,jbiny,allSlice);
            // jy-jy
            jbiny = (ibiny + jdbiny) % nybin;
            temp(jdbinx, ibinx, jdbiny, ibiny, 1, 1, allSlice)
             = conj(jy(0,ibinx,ibiny,allSlice)) * jy(0,jbinx,jbiny,allSlice);
            jbiny = (ibiny - jdbiny + nybin) % nybin;
            temp(jdbinx, ibinx, 2*nydbin-1-jdbiny, ibiny, 1, 1, allSlice)
             = conj(jy(0,ibinx,ibiny,allSlice)) * jy(0,jbinx,jbiny,allSlice);
            // jx-jy
            jbiny = (ibiny + jdbiny) % nybin;
            temp(jdbinx, ibinx, jdbiny, ibiny, 0, 1, allSlice)
             = conj(jx(0,ibinx,ibiny,allSlice)) * jy(0,jbinx,jbiny,allSlice);
            jbiny = (ibiny - jdbiny + nybin) % nybin;
            temp(jdbinx, ibinx, 2*nydbin-1-jdbiny, ibiny, 0, 1, allSlice)
             = conj(jx(0,ibinx,ibiny,allSlice)) * jy(0,jbinx,jbiny,allSlice);
            // jy-jx
            jbiny = (ibiny + jdbiny) % nybin;
            temp(jdbinx, ibinx, jdbiny, ibiny, 1, 0, allSlice)
             = conj(jy(0,ibinx,ibiny,allSlice)) * jx(0,jbinx,jbiny,allSlice);
            jbiny = (ibiny - jdbiny + nybin) % nybin;
            temp(jdbinx, ibinx, 2*nydbin-1-jdbiny, ibiny, 1, 0, allSlice)
             = conj(jy(0,ibinx,ibiny,allSlice)) * jx(0,jbinx,jbiny,allSlice);
          }
        }
        jbinx = (ibinx - jdbinx + nxbin) % nxbin;
        for(int ibiny=0; ibiny<nybin; ++ibiny) {
          for(int jdbiny=1; jdbiny<nydbin; ++jdbiny) {
            // jx-jx
            int jbiny = (ibiny + jdbiny) % nybin;
            temp(2*nxdbin-1-jdbinx, ibinx, jdbiny, ibiny, 0, 0, allSlice)
             = conj(jx(0,ibinx,ibiny,allSlice)) * jx(0,jbinx,jbiny,allSlice);
            jbiny = (ibiny - jdbiny + nybin) % nybin;
            temp(2*nxdbin-1-jdbinx, ibinx, 2*nydbin-1-jdbiny, ibiny, 0, 0, allSlice)
             = conj(jx(0,ibinx,ibiny,allSlice)) * jx(0,jbinx,jbiny,allSlice);
            // jy-jy
            jbiny = (ibiny + jdbiny) % nybin;
            temp(2*nxdbin-1-jdbinx, ibinx, jdbiny, ibiny, 1, 1, allSlice)
             = conj(jy(0,ibinx,ibiny,allSlice)) * jy(0,jbinx,jbiny,allSlice);
            jbiny = (ibiny - jdbiny + nybin) % nybin;
            temp(2*nxdbin-1-jdbinx, ibinx, 2*nydbin-1-jdbiny, ibiny, 1, 1, allSlice)
             = conj(jy(0,ibinx,ibiny,allSlice)) * jy(0,jbinx,jbiny,allSlice);
            // jx-jy
            jbiny = (ibiny + jdbiny) % nybin;
            temp(2*nxdbin-1-jdbinx, ibinx, jdbiny, ibiny, 0, 1, allSlice)
             = conj(jx(0,ibinx,ibiny,allSlice)) * jy(0,jbinx,jbiny,allSlice);
            jbiny = (ibiny - jdbiny + nybin) % nybin;
            temp(2*nxdbin-1-jdbinx, ibinx, 2*nydbin-1-jdbiny, ibiny, 0, 1, allSlice)
             = conj(jx(0,ibinx,ibiny,allSlice)) * jy(0,jbinx,jbiny,allSlice);
            // jy-jx
            jbiny = (ibiny + jdbiny) % nybin;
            temp(2*nxdbin-1-jdbinx, ibinx, jdbiny, ibiny, 1, 0, allSlice)
             = conj(jy(0,ibinx,ibiny,allSlice)) * jx(0,jbinx,jbiny,allSlice);
            jbiny = (ibiny - jdbiny + nybin) % nybin;
            temp(2*nxdbin-1-jdbinx, ibinx, 2*nydbin-1-jdbiny, ibiny, 1, 0, allSlice)
             = conj(jy(0,ibinx,ibiny,allSlice)) * jx(0,jbinx,jbiny,allSlice);
          }
        }
      }
    }
    for (int ibinx=0; ibinx<nxbin; ++ibinx) {
      for (int ibiny=0; ibiny<nybin; ++ibiny) {
        // jx-jx
        temp(0, ibinx, 0, ibiny, 0, 0, allSlice)
          = conj(jx(0,ibinx,ibiny,allSlice)) * jx(0,ibinx,ibiny,allSlice);
        // jy-jy
        temp(0, ibinx, 0, ibiny, 1, 1, allSlice) 
          = conj(jy(0,ibinx,ibiny,allSlice)) * jy(0,ibinx,ibiny,allSlice);
        // jx-jy
        temp(0, ibinx, 0, ibiny, 0, 1, allSlice)
          = conj(jx(0,ibinx,ibiny,allSlice)) * jy(0,ibinx,ibiny,allSlice);
        // jy-jx
        temp(0, ibinx, 0, ibiny, 1, 0, allSlice)
          = conj(jy(0,ibinx,ibiny,allSlice)) * jx(0,ibinx,ibiny,allSlice);
      }
    }
    value += real(temp(allBin, allBin, allBin, allBin, alldir, alldir, 
                  blitz::Range(0,nfreq-1))) / (tau*nslice);
    norm+=1;
  }
}
