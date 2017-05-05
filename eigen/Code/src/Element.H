#ifndef _ELEMENT_H_
#define _ELEMENT_H_
#define VERTICES 3
#include <array>
using namespace std;

class Element
{
public: 
  Element();
  /// Constructor. 
  Element(array<int, VERTICES> a_vertices);

  /// local indexing to get nodeNumber.
  const int& operator[](const int& a_localNodeNumber) const;
  const array<int, VERTICES>& vertices() const;

private:
  array<int, VERTICES> m_vertices;
};
#endif
