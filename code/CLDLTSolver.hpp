#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include "CLinearBase.h"
template <typename T>
class CLDLTSolver
{
public:
    CLDLTSolver()
    {
		m_x = 0;
		m_b = 0;
		m_eigen_A = 0;
		m_decomposed = false;
    }
    ~CLDLTSolver()
    {
        if(m_dim != 0)
        {
            clear();
        }
    }
    void init(int dimension)
    {
        m_dim = dimension;
		m_b = new T[dimension];
		m_x = new T[dimension];
    }
    void add_data(int row, int col, T data)
    {
		m_vec_triplets.push_back(Eigen::Triplet<T>(row, col, data));
    }
	void decompose()
    {
		m_eigen_A = new Eigen::SparseMatrix<T>(m_dim, m_dim);
		m_eigen_A->setFromTriplets(m_vec_triplets.begin(), m_vec_triplets.end());
		m_vec_triplets.clear();
		m_eigen_solver = new Eigen::SimplicialLDLT<Eigen::SparseMatrix<T> >;
		m_eigen_solver->compute(*m_eigen_A);
		m_decomposed = true;
    }
    void solve(int correct_times = 0)
    {
		int i;
		Eigen::Matrix<T, -1, 1> eigen_b(m_dim), eigen_x(m_dim), eigen_e(m_dim), eigen_y(m_dim);
		for (i = 0; i < m_dim; ++i)
		{
			eigen_b(i) = m_b[i];
		}
		eigen_x = m_eigen_solver->solve(eigen_b);
		while (correct_times > 0)
		{
			eigen_e = eigen_b - (*m_eigen_A) * eigen_x;
			eigen_y = m_eigen_solver->solve(eigen_e);
			eigen_x = eigen_x + eigen_y;
			--correct_times;
		}
		for (i = 0; i < m_dim; ++i)
		{
			m_x[i] = eigen_x(i);
		}
    }
    void clear()
    {
		m_dim = 0;
		delete[]m_b;
		m_b = 0;
		delete[]m_x;
		m_x = 0;
		if (m_decomposed)
		{
			delete m_eigen_solver;
			m_eigen_solver = 0;
			delete m_eigen_A;
			m_eigen_A = 0;
		}
		else
		{
			m_vec_triplets.clear();
		}
		m_dim = 0;
    }
public:
	T *m_x;
	T *m_b;
	int m_dim;
protected:
	bool m_decomposed;
	std::vector<Eigen::Triplet<T> > m_vec_triplets;
	Eigen::SparseMatrix<T> *m_eigen_A;
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<T> > *m_eigen_solver;
private:
};
