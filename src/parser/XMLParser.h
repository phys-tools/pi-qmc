#ifndef __XMLParser_h_
#define __XMLParser_h_
#include <string>
#include <libxml/tree.h>
#include <libxml/xpath.h>
#include <cstdlib>
#include <blitz/tinyvec.h>
#include <config.h>

/** Base class for classes that parse XML data.  
  * Includes helper functions.
  * @version $Revision$
  * @author John Shumway */
class XMLParser {
public:
  typedef blitz::TinyVector<double,NDIM> Vec; 
  typedef blitz::TinyVector<int,NDIM> IVec; 
  /// Destructor.
  virtual ~XMLParser() {};
  /// Parse some xml.
  virtual void parse(const xmlXPathContextPtr& ctxt) {};
public:
  /// Get the name of the node.
  static std::string getName(const xmlNodePtr& node);
  /// Get a double valued attribute.
  static double getDoubleAttribute(const xmlNodePtr& node,
                                   const std::string& attName);
  /// Get an integer valued attribute.
  static int getIntAttribute(const xmlNodePtr& node,
                             const std::string& attName);
  /// Get a string valued attribute.
  static std::string getStringAttribute(const xmlNodePtr& node,
                                        const std::string& attName);
  /// Get a bool valued attribute.
  static bool getBoolAttribute(const xmlNodePtr& node,
                               const std::string& attName);
  /// Get a vector valued attribute.
  static Vec getVecAttribute(const xmlNodePtr& node,
                             const std::string& attName);
  /// Get a vector valued attribute.
  static IVec getIVecAttribute(const xmlNodePtr& node,
                               const std::string& attName);
  /// Get the base of a link.
  static std::string getLinkBase(const xmlNodePtr& node,
                                 const std::string& attName="href");

  /// Get the path of a link.
  static std::string getLinkPath(const xmlNodePtr& node,
                                 const std::string& attName="href");
  /// Letters associated with directions in input file.
  static const std::string dimName;
};
#endif
