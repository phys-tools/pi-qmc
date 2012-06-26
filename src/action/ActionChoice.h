#ifndef __ActionChoice_h_
#define __ActionChoice_h_
class MultiLevelSampler;
class DisplaceMoveSampler;
class SectionChooser;
class Paths;
class ModelState;
class EnumeratedModelState;
#include "CompositeAction.h"
#include "LinkSummable.h"
#include <vector>

class ActionChoiceBase {
public:
    ActionChoiceBase() :
            modelState(0) {
    }
    virtual ~ActionChoiceBase() {}
    virtual double getActionChoiceDifference(const Paths&, int jmodel)=0;
    ModelState& getModelState() {
        return *modelState;
    }
    const ModelState& getModelState() const {
        return *modelState;
    }
protected:
    double actionDifference;
    ModelState *modelState;
};

/** Action class for representing a choice of action.
 * Used to compare free energies of two different path integrals.
 * @version $Revision$
 * @author John Shumway. */
class ActionChoice: public CompositeAction,
        public LinkSummable,
        public ActionChoiceBase {
public:
    /// Constructor, optionally providing the number of slots to reserve.
    ActionChoice(const int nreserve = 0);
    /// Virtual destructor deletes all Action objects.
    virtual ~ActionChoice();
    /// Calculate the difference in action.
    virtual double getActionDifference(const SectionSamplerInterface&,
            int level);
    /// Calculate the difference in action.
    virtual double getActionDifference(const Paths&, const VArray &displacement,
            int nmoving, const IArray &movingIndex, int iFirstSlice,
            int iLastSlice);
    /// Calculate the difference in action.
    virtual double getActionChoiceDifference(const Paths&, int jmodel);
    /// Calculate the total action.
    virtual double getTotalAction(const Paths&, const int level) const;
    /// Calculate the action and derivatives at a bead.
    virtual void getBeadAction(const Paths&, const int ipart, const int islice,
            double& u, double& utau, double& ulambda, Action::Vec& fm,
            Action::Vec& fp) const;
    /// Add an action object.
    void addAction(Action* a) {
        actions.push_back(a);
    }
    /// Initialize for a sampling section.
    virtual void initialize(const SectionChooser&);
    /// Accept last move.
    virtual void acceptLastMove();
    /// Get the number of model choices.
    virtual int getModelCount() const {
        return actions.size();
    }
    /// Initialize the calculation.
    virtual void initCalc(const int nslice, const int firstSlice);
    /// Add contribution from a link.
    virtual void handleLink(const LinkSummable::Vec& start,
            const LinkSummable::Vec& end, const int ipart, const int islice,
            const Paths&);
    /// Finalize the calculation.
    virtual void endCalc(const int nslice);
private:
    EnumeratedModelState *enumModelState;
    int jmodel;
};
#endif
