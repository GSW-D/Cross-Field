#include "CMinimalNormSolver.hpp"
template <typename T>
struct CLeastSquareSolver
{
	CLeastSquareSolver()
	{
		m_b = 0;
		m_x = 0;
		m_dim = 0;
		m_decomposed = false;
		m_mat_size = 0;
	}
	~CLeastSquareSolver()
	{
		if (m_dim > 0)
		{
			clear();
		}
	}
	void init(int row, int col)
	{
		m_solver.init(col);
		m_x = m_solver.m_x;
		m_b = new T[row];
		m_dim = col;
	}
	void clear()
	{
		if (m_decomposed)
		{
			delete[]m_ordered_mat;
		}
		else
		{
			m_coef_mat.clear();
		}
		delete[]m_b;
		m_dim = 0;
		m_solver.clear();
	}
	void add_coef(int row, int col, T data)
	{
		MatElem<T> new_elem;
		new_elem.row = row;
		new_elem.col = col;
		new_elem.data = data;
		m_coef_mat.append(new_elem);
	}
	void decompose()
	{
		m_mat_size = sort_mat_row(m_coef_mat, &m_ordered_mat);
		fill_tensor();
		m_solver.decompose();
		m_decomposed = true;
	}
	void solve(int correct_times = 0)
	{
		MatElem<T>* mat_elem, * mat_end;
		int i;
		mat_end = m_ordered_mat + m_mat_size;
		for (i = 0; i < m_dim; ++i)
		{
			m_solver.m_b[i] = 0;
		}
		for (mat_elem = m_ordered_mat; mat_elem != mat_end; ++mat_elem)
		{
			m_solver.m_b[mat_elem->col] += mat_elem->data * m_b[mat_elem->row];
		}
		m_solver.solve(correct_times);
	}
private:
	void fill_tensor()
	{
		MatElem<T>* m0, * m1, * mat_end, * i, * j;
		m0 = m_ordered_mat;
		m1 = m_ordered_mat;
		mat_end = m_ordered_mat + m_mat_size;
		while (m0 < mat_end)
		{
			do
			{
				++m1;
			} while (m1 < mat_end && m1->row == m0->row);
			for (i = m0; i < m1; ++i)
			{
				for (j = m0; j < m1; ++j)
				{
					m_solver.add_data(i->col, j->col, i->data * j->data);
				}
			}
			m0 = m1;
		}
	}
public:
	T* m_b;
	T* m_x;
	int m_mat_size;
	int m_dim;
	MatElem <T>* m_ordered_mat;
private:
	bool m_decomposed;
	Pool<MatElem<T> > m_coef_mat;
	CLDLTSolver<T> m_solver;
};