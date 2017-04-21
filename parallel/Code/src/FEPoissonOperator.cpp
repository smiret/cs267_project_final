#include "FEPoissonOperator.H"
FEPoissonOperator::FEPoissonOperator(){};
FEPoissonOperator::FEPoissonOperator(const FEGrid& a_grid)
{
  m_grid = a_grid;
  int nInt = a_grid.getNumInteriorNodes();
  SparseMatrix L(nInt,nInt);
  int nElt = m_grid.getNumElts();
  for(int i = 0; i < nElt; i++)
    {
      //const Element elem = m_grid.element(i);
      //const array<int, VERTICES> vert = elem.vertices();
      double area = m_grid.elementArea(i);
      for(int j = 0; j < VERTICES; j++)
        {
          const Node node1 = m_grid.getNode(i,j);
          if(node1.isInterior())
            {
              array<double, DIM> grad1 = m_grid.gradient(i,j);
              const int nodeID1 = node1.getInteriorNodeID();
              for(int k = 0; k < VERTICES; k++)
                {
                  const Node node2 = m_grid.getNode(i,k);
                  if(node2.isInterior())
                    {
                      array<double, DIM> grad2 = m_grid.gradient(i,k);
                      const int nodeID2 = node2.getInteriorNodeID();
                      L[{{nodeID1,nodeID2}}] += ((grad1[0] * grad2[0])+(grad1[1] * grad2[1])) * area;
                    }
                }
            }
        }
    }
  m_matrix = L;
};
void FEPoissonOperator::makeRHS(vector<double> & a_rhsAtNodes, const vector<double> & a_FCentroids) const
{
  int nInt = m_grid.getNumInteriorNodes();
  vector<double> b(nInt);
  int nElt = m_grid.getNumElts();
  double oneThird = 1.0/3.0;
  for(int i = 0; i < nElt; i++)
    {
      double area = m_grid.elementArea(i);
      for(int j = 0; j < VERTICES; j++)
        {
          const Node node = m_grid.getNode(i,j);
          const int nodeID = node.getInteriorNodeID();
          if(node.isInterior())
            {
              b[nodeID] += area * a_FCentroids[i] * oneThird;
            }
        }
    }
  a_rhsAtNodes = b;
};
const FEGrid& FEPoissonOperator::getFEGrid() const
{
  return m_grid;
};
const SparseMatrix& FEPoissonOperator::matrix() const
{
  return m_matrix;
};