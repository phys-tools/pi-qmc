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
#include "StdoutReportBuilder.h"
#include "ScalarEstimator.h"
#include "AccRejEstimator.h"
#include "ArrayBlockedEstimator.h"
#include "EstimatorManager.h"
#include <ctime>
#include <string.h>
#include <iostream>

void StdoutReportBuilder::startWritingGroup(EstimatorManager& manager) {
  nstep=manager.nstep;
  istep=0;
  int n=manager.estimator.size();
  sum.resize(n); sum2.resize(n); norm.resize(n);
  sum=0; sum2=0; norm=0; 
}

void StdoutReportBuilder::writeStep(EstimatorManager& manager) {
  std::time_t rawtime;
  std::time ( &rawtime );
  char * myTime = std::ctime (&rawtime);
  myTime[strlen(myTime)-1] = '\0'; //Hackish, gets rid of a trailing newline
  std::cout << "********** Block " << istep+1 << " of " << nstep << " blocks **********" ; //<< std::endl;
  std::cout << " (" << myTime << ")" << std::endl; 
  iscalar=0;
  for (EstimatorManager::EstimatorIter est=manager.estimator.begin();
       est!=manager.estimator.end(); ++est) {
    (*est)->reportStep(*this);
  }
  std::cout << std::endl;
  if (++istep == nstep) {}
}

void StdoutReportBuilder::reportScalarStep(const ScalarEstimator& est) {
  double value=est.getValue();
  sum(iscalar)+=value; sum2(iscalar)+=value*value; norm(iscalar)+=1;
  std::cout << est.getName();
  if (est.getUnitName()!="") std::cout << " (" << est.getUnitName() << ")";
  std::cout << ": " << value << ", " 
            << "Av=" << sum(iscalar)/(norm(iscalar))
            << " +-" << sqrt(sum2(iscalar)-sum(iscalar)*sum(iscalar)
                            /(norm(iscalar)))/(norm(iscalar)-1)
            << std::endl;
  ++iscalar;
}

void StdoutReportBuilder::reportAccRejStep(const AccRejEstimator& est) {
  std::cout << est.getName() << std::endl;
  int nlevel=est.getNLevel();
  const AccRejEstimator::IArray& nacc(est.getNAccept());
  const AccRejEstimator::IArray& ntrl(est.getNTrial());
  for (int i=nlevel-1; i>=0; --i)  {
    std::cout << "Level " << i << ": " 
              << nacc(i) << "/" << ntrl(i) << " "
              << nacc(i)/(float)(ntrl(i)==0?1:ntrl(i)) << std::endl;
  }
  std::cout << "Total:   " << nacc(0)<< "/" << ntrl(nlevel-1) << " " 
         << nacc(0)/(float)(ntrl(nlevel-1)==0?1:ntrl(nlevel-1)) << std::endl;
}

void 

StdoutReportBuilder::reportArrayBlockedStep(const ArrayBlockedEstimator& est) {
  std::cout << "(measured " << est.getName() << ")" << std::endl;
}
