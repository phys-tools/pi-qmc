//$Id: SpinChoiceFixedNodeAction.h 393 2011-08-04 14:39:41Z john.shumwayjr $
/*  Copyright (C) 2011 John B. Shumway, Jr. and Jianheng Liu

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
#ifndef __SpinChoiceFixedNodeAction_h_
#define __SpinChoiceFixedNodeAction_h_

#include "FixedNodeAction.h"
#include "ActionChoice.h"
#include "LinkSummable.h"
class SpinModelState;
class MPIManager;

/**
@version $Revision: 393 $
@author John Shumway and Jianheng Liu*/
class SpinChoiceFixedNodeAction : public FixedNodeAction,
                                  public ActionChoiceBase,
                                  public LinkSummable {
public:
    typedef blitz::Array<int,1> IArray;
    SpinChoiceFixedNodeAction(const SimulationInfo&, int initial,
	const Species&, NodeModel*,
        bool withNodalAction, bool useDistDerivative, int maxlevel,
        bool useManyBodyDistance,const MPIManager* mpi);

    virtual ~SpinChoiceFixedNodeAction();

    virtual double getActionDifference(const Paths &paths, int ipart);

//    virtual double getTotalAction(const Paths& paths) const;

    virtual void initCalc(const int nslice, const int firstSlice);
    
    virtual void handleLink(const LinkSummable::Vec &start, 
        const LinkSummable::Vec &end, int ipart, 
        int islice, const Paths &paths);

private:
    SpinModelState *spinModelState;
    double totalAction;
    const MPIManager* mpi;
};
#endif
