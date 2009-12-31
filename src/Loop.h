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
#ifndef __Loop_h_
#define __Loop_h_
#ifdef ENABLE_MPI 
#include <mpi.h>
#endif
#include "stats/MPIManager.h"

#include "CompositeAlgorithm.h"
#include <iostream> ////////// dak
#include <string>   ////////// sak
#include <time.h>
/** Algorithm class for looping.
 * @version $Revision$
 * @author John Shumway */
class Loop : public CompositeAlgorithm {
public:
  Loop(const int nrepeat, const int totalSimTime, const std::string timer,
       const MPIManager* mpi, const int nsteps=0) 
    : CompositeAlgorithm(nsteps), nrepeat(nrepeat), totalSimTime(totalSimTime),
      timer(timer), mpi(mpi) {
  }
    
    virtual ~Loop() {}

    virtual void run() {

      // std :: cout << mpi->getCloneID()<<" cid. "<<mpi->getWorkerID()<<" iw. "<<timer<<" "<<totalSimTime<<std :: endl;

      //mpi->getWorkerComm().Barrier();
      
      // print out timeinfo for Main or sub loops. totalsimtime default value is 12hrs in case it is not set.
      if (timer =="Main" || totalSimTime>0){
	time_t startSim, elapsedTime;
	time (&startSim);
	double dif=0;
	double oldTime=0;
	double dt = 0;
	for(int i=0; (dif+dt) < totalSimTime && i<nrepeat; ++i ){

#ifdef ENABLE_MPI
	  if (mpi) {
	    mpi->getWorkerComm().Barrier();
	    if (mpi->isCloneMain())  {
	      mpi->getCloneComm().Barrier();
	    }
	  }
#endif
	  CompositeAlgorithm::run();
	  time (&elapsedTime);
	  dif = difftime (elapsedTime,startSim);
	  dt = dif-oldTime;
	  oldTime=dif;

#ifdef ENABLE_MPI
	  if ( mpi->isMain()) {
	    printElapsedTime(dif);	  
	    printAlgorithmTime(dt);
	  }
	
#else
	printElapsedTime(dif);	  
	printAlgorithmTime(dt);
#endif

	//	std :: cout << mpi->getCloneID()<<" cid. "<<mpi->getWorkerID()<<" iw. "<<dt<<" "<<  dif<<std :: endl<<std::flush;
	}
      } else if (timer!="" && timer !="Main"){
	time_t startSim, elapsedTime;
	time (&startSim);
	double dif=0;
	double oldTime=0;
	double dt = 0;
#ifdef ENABLE_MPI
	  //db sak
	  if (mpi) {
	    mpi->getWorkerComm().Barrier();
	    if (mpi->isCloneMain())  {
	      mpi->getCloneComm().Barrier();
	    }
	  }
#endif
	for(int i=0; i<nrepeat; ++i )  CompositeAlgorithm::run();
	time (&elapsedTime);
	dif = difftime (elapsedTime,startSim);
	dt = dif-oldTime ;
	oldTime=dif;

#ifdef ENABLE_MPI
	if ( mpi->isMain()) printAlgorithmTime(dt);
#else
	printAlgorithmTime(dt);
#endif
      }else {
 	//No timeinfo printed. timer not used. Just use nrepeat
      	if (timer=="" && totalSimTime==0){
#ifdef ENABLE_MPI
	    //db sak
	    if (mpi) {
	      mpi->getWorkerComm().Barrier();
	      if (mpi->isCloneMain())  {
		mpi->getCloneComm().Barrier();
	      }
	    }
#endif
	  for(int i=0;  i<nrepeat; ++i )	 {
	    CompositeAlgorithm::run();

	  }
	}
      }
    }
    
    void printElapsedTime(double dif){
      int d = int(dif);
      int hr = d/3600;
      int mn = ((d)%3600)/60;
      int sc = ((d)%3600)%60;
      std :: cout << "Time elapsed from start of simulation :: "<<hr<<" hour(s), "<< mn <<" min(s), "<<sc<<" sec(s), or total time in seconds :: "<<dif<<" secs."<<std :: endl<<std::flush;
    }
    
    void printAlgorithmTime(double dt){
      int d = int(dt);
      int mn = d/60;
      int sc = (d)%60;
      std :: cout << "Time spent inside loop ["<<timer<< "] ..... :: "<<mn<<" min(s), "<<sc<<" sec(s), or total time in seconds :: "<<dt<<" secs."<<std :: endl<<std::flush;
    }
    
 private:
    const int nrepeat;
    const int totalSimTime;
    const std::string timer;
    const MPIManager* mpi;
};
#endif
