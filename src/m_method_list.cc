/*!
  \file   m_method_list.cc
  \author Stefan Buehler <sbuehler@uni-bremen.de>
  \date   Fri Nov 16 10:11:16 2001
  
  \brief  Workspace methods related to method lists.
  
*/

#include <map>
#include <stdexcept>
#include "arts.h"
#include "mystring.h"
#include "array.h"

//! Set up a method list.
/*!
  A method list just contains indices (in md_data) of methods
  intended for sequential execution. Only methods without keyword
  arguments are allowed. It is the task of this method to
  set this up. For example, it must be checked, whether the given
  names really correspond to methods.
  
  \param ml Output: The list of methods.
  \param ml_name The name of the variable for which the method has been called.
  \param methods An array of names of methods.
*/
void MethodListDefine(// WS Generic Output:
                      ArrayOfIndex& ml,
                      // WS Generic Output Names:
                      const String& ml_name,
                      // Control Parameters:
                      const ArrayOfString& methods)
{
  // The method data map. We can use this to find out the index of
  // each method.
  extern const std::map<String, Index> MdMap;

  // Make ml the right size:
  ml.resize(methods.nelem());

  // Loop through the given method names:
  for ( Index i=0; i<methods.nelem(); ++i )
    {
      // Find method id:
      const map<String, Index>::const_iterator im = MdMap.find(methods[i]);
      if ( im == MdMap.end() )
        {
          ostringstream os;
          os << "\"" << methods[i] << "\" is not a valid method. "
             << "Try \"arts -m all\" to\n"
             << "get a list of all ARTS methods.";
          throw runtime_error( os.str() );
        }

      // Assign to our method list:
      ml[i] = im->second;
    }
}
