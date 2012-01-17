#ifndef __DoubleCollectiveSectionSampler_h_
#define __DoubleCollectiveSectionSampler_h_
template <int TDIM> class Beads;
class DoubleAction;
class Action;
class BeadFactory;
class Paths;
class DoubleSectionChooser;

class DoubleCollectiveSectionSampler : public DoubleMLSampler {
public:
  /// Constructor.
  DoubleCollectiveSectionSampler(const int nmoving, Paths&, DoubleSectionChooser&,
		  Action*, DoubleAction*, const bool sampleBoth, const int nrepeat,
		  const BeadFactory&);
/// Destructor.
  virtual ~DoubleCollectiveSectionSampler();
  /// Run the sampler.
  virtual void run();
private:
};
#endif
