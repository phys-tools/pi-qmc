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
#ifndef __CompositeAction_h_
#define __CompositeAction_h_
class MultiLevelSampler;
class DisplaceMoveSampler;
class SectionChooser;
class Paths;
#include "Action.h"
#include <vector>

/** Action class for composing multiple Action classes.
  * @version $Revision$
  * @author John Shumway. */
class CompositeAction : public Action {
protected:
  typedef std::vector<Action*> ActionContainer;
  typedef ActionContainer::iterator ActionIter;
  typedef ActionContainer::const_iterator ConstActionIter;
public:
  /// Constructor, optionally providing the number of slots to reserve.
  CompositeAction(const int nreserve=0);
  /// Virtual destructor deletes all Action objects.
  virtual ~CompositeAction();
  /// Calculate the difference in action.
  virtual double getActionDifference(const MultiLevelSampler&,
                                     const int level);
  /// Calculate the difference in action.
  virtual double getActionDifference(const Paths&, const VArray &displacement,
    int nmoving, const IArray &movingIndex, int iFirstSlice, int iLastSlice);
  /// Calculate the total action.
  virtual double getTotalAction(const Paths&, const int level) const;
  /// Calculate the action and derivatives at a bead.
  virtual void getBeadAction(const Paths&, const int ipart, const int islice,
       double& u, double& utau, double& ulambda, Vec& fm, Vec& fp) const;
  /// Add an action object.
  void addAction(Action* a) {actions.push_back(a);}
  /// Initialize for a sampling section.
  virtual void initialize(const SectionChooser&);
  /// Accept last move.
  virtual void acceptLastMove();
  /// Returns pointer to an Action of type type, otherwise returns null pointer.
  virtual const Action* getAction(const std::type_info &type) const;
  /// Returns the number of action objects.
  int getCount() {return actions.size();}
protected:
  /// Pointers to the Action objects.
  ActionContainer actions;
};
#endif
