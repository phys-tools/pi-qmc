#ifndef __Units_h_
#define __Units_h_
#include <string>
#include <map>
class Paths;
/// Class for unit conversions.
/// @author John Shumway
class Units {
public:
  /// Construct by giving the name.
  Units(const std::string& eunit, const std::string& lunit);
  /// Destructor.
  ~Units();
  /// Get conversion factor for length input with unit.
  double getLengthScaleIn(const std::string& unit, int factor=1) const;
  /// Get conversion factor for energy input with unit.
  double getEnergyScaleIn(const std::string& unit) const;
  /// Get conversion factor for mass input with unit.
  double getMassScaleIn(const std::string& unit) const;
protected:
  typedef std::map<std::string,double> unitMap;
  /// The name of the internal energy unit.
  const std::string internalEnergyUnit;
  /// The name of the internal length unit.
  const std::string internalLengthUnit;
  /// The conversion factors to different length units.
  unitMap lengthOut;
  /// The conversion factors to different energy units.
  unitMap energyOut;
  /// The conversion factors to different mass units.
  unitMap massOut;
};
#endif
