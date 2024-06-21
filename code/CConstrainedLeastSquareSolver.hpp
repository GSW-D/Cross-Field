#include "CLeastSquareSolver.hpp"
template <typename T>
class CConstrainedLeastSquareSolver
{
public:
	CConstrainedLeastSquareSolver()
	{
		m_dim = 0;
		m_hard = 0;
		m_soft = 0;
		m_ordered_mat = 0;
		m_b = 0;
		m_x = 0;
		m_q = 0;
		m_e = 0;
	}
	~CConstrainedLeastSquareSolver()
	{
		if (m_dim != 0)
		{
			clear();
		}
	}
	void clear()
	{
		if (sorted)
		{
			delete[]m_ordered_mat;
		}
		else
		{
			m_coef_mat.clear();
		}
		delete[]m_x;
		delete[]m_b;
		delete[]m_q;
		m_dim = 0;
		m_solver.clear();
	}
	void init(int dim, int hard, int soft)
	{
		m_dim = dim;
		m_soft = soft;
		m_hard = hard;
		m_solver.init(hard, dim);
		m_q = new T[hard];//m_solver.m_b;
		m_b = new T[soft];
		m_x = new T[dim];
	}
	void add_hard_coef(int row, int col, T data)
	{
		m_solver.add_coef(row, col, data);
	}
	void add_soft_coef(int row, int col, T data)
	{
		MatElem<T> mat_elem;
		mat_elem.row = row;
		mat_elem.col = col;
		mat_elem.data = data;
		m_coef_mat.append(mat_elem);
	}
	void decompose()
	{
		m_solver.decompose();
	}
	void solve(T max_e = 1.0e-5, int max_iter = 0)
	{
		T d, t, s, e;
		T *r, *u, *v, *w, *x;
		int i, iter;
		if (max_iter == 0)
		{
			max_iter = int(sqrt(double(m_dim)) + 0.5);
		}
		r = new T[m_dim];
		u = new T[m_dim];
		v = new T[m_dim];
		w = new T[m_soft];
		x = new T[m_dim];
		m_mat_size = sort_mat(m_coef_mat, &m_ordered_mat);
		sorted = true;
		for (i = 0; i < m_hard; ++i)
		{
			m_solver.m_b[i] = m_q[i];
		}
		m_solver.solve();
		for (i = 0; i < m_dim; ++i)
		{
			x[i] = m_x[i] = m_solver.m_x[i];
		}
		for (i = 0; i < m_mat_size; ++i)
		{
			m_b[m_ordered_mat[i].row] -= m_ordered_mat[i].data * x[m_ordered_mat[i].col];
		}
		for (i = 0; i < m_dim; ++i)
		{
			r[i] = 0;
		}
		for (i = 0; i < m_mat_size; ++i)
		{
			r[m_ordered_mat[i].col] += m_ordered_mat[i].data * m_b[m_ordered_mat[i].row];
		}
		for (i = 0; i < m_hard; ++i)
		{
			m_solver.m_b[i] = 0;
		}
		for (i = 0; i < m_solver.m_mat_size; ++i)
		{
			m_solver.m_b[m_solver.m_ordered_mat[i].row] += m_solver.m_ordered_mat[i].data * r[m_solver.m_ordered_mat[i].col];
		}
		m_solver.solve();
		e = 0;
		for (i = 0; i < m_dim; ++i)
		{
			r[i] -= m_solver.m_x[i];
			v[i] = 0;
			e += r[i] * r[i];
		}
		m_e = e;
		s = 0;
		iter = 0;
		while (iter < max_iter || m_e > max_e)
		{
			for (i = 0; i < m_dim; ++i)
			{
				v[i] = r[i] - v[i] * s;
				u[i] = 0;
			}
			for (i = 0; i < m_soft; ++i)
			{
				w[i] = 0;
			}
			for (i = 0; i < m_mat_size; ++i)
			{
				w[m_ordered_mat[i].row] += m_ordered_mat[i].data * v[m_ordered_mat[i].col];
			}
			for (i = 0; i < m_mat_size; ++i)
			{
				u[m_ordered_mat[i].col] += m_ordered_mat[i].data * w[m_ordered_mat[i].row];
			}
			for (i = 0; i < m_hard; ++i)
			{
				m_solver.m_b[i] = 0;
			}
			for (i = 0; i < m_solver.m_mat_size; ++i)
			{
				m_solver.m_b[m_solver.m_ordered_mat[i].row] += m_solver.m_ordered_mat[i].data * u[m_solver.m_ordered_mat[i].col];
			}
			m_solver.solve();
			d = 0;
			for (i = 0; i < m_dim; ++i)
			{
				u[i] -= m_solver.m_x[i];
				d += u[i] * v[i];
			}
			t = e / d;
			e = 0;
			s = 0;
			for (i = 0; i < m_dim; ++i)
			{
				x[i] += v[i] * t;
				r[i] -= u[i] * t;
				e += r[i] * r[i];
				s += u[i] * r[i];
			}
			if (e < m_e)
			{
				m_e = e;
				for (i = 0; i < m_dim; ++i)
				{
					m_x[i] = x[i];
				}
				iter = 0;
			}
			else
			{
				++iter;
			}
			s /= d;
			//++iter;
		}
		delete[]x;
		delete[]w;
		delete[]v;
		delete[]u;
		delete[]r;
	}
	void reset_soft()
	{
		if (sorted)
		{
			delete[]m_ordered_mat;
		}
		else
		{
			m_coef_mat.clear();
		}
	}
public:
	T *m_b, *m_x, *m_q;
private:
	CMinimalNormSolver<T> m_solver;
	MatElem<T> *m_ordered_mat;
	Pool<MatElem <T> > m_coef_mat;
	T m_e;
	int m_dim, m_mat_size, m_soft, m_hard;
	bool sorted;
};