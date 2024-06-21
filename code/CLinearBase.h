#ifndef  CLINEARBASE_H
#define CLINEARBASE_H
template<typename T>
class Pool
{
public:
	Pool()
	{
		int i;
		m_size = 0;
		m_group = 0;
		m_pos = 0;
		m_capicity = 1;
		m_indices[0] = &m_first_element;
		for (i = 1; i < 32; ++i)
		{
			m_indices[i] = 0;
		}
	}
	T* append(T data)
	{
		if (m_pos == m_capicity)
		{
			++m_group;
			m_capicity <<= 1;
			m_indices[m_group] = new T[m_capicity];
			m_pos = 0;
		}
		m_indices[m_group][m_pos] = data;
		++m_size;
		++m_pos;
		return &(m_indices[m_group][m_pos - 1]);
	}
	void exact_data(T *output)
	{
		int i, group, pos, capicity;
		group = 0;
		pos = 0;
		capicity = 1;
		for (i = 0; i < m_size; ++i)
		{
			output[i] = m_indices[group][pos];
			++pos;
			if (pos == capicity)
			{
				capicity <<= 1;
				++group;
				pos = 0;
			}
		}
	}
	void clear()
	{
		int i;
		i = 1;
		while (m_indices[i] != 0 && i < 32)
		{
			delete[]m_indices[i];
			m_indices[i] = 0;
			++i;
		}
		m_size = 0;
		m_group = 0;
		m_pos = 0;
		m_capicity = 1;
	}
public:
	int m_size;
private:
	int m_group, m_pos, m_capicity;
	T m_first_element;
	T *m_indices[32];
};
template <typename T>
struct MatElem
{
	int row;
	int col;
	T data;
};
template <typename T>
void merge_sort(T *data, int *order, int length)
{
	int *buffer[2];
	int source, dest, step;
	int *m0, *m1, *m2, *i, *j, *k, *end;
	buffer[0] = new int[length];
	buffer[1] = new int[length];
	end = buffer[0] + length;
	k = order;
	for (i = buffer[0]; i  != end; ++i)
	{
		*i = *k;
		++k;
	}
	source = 0;
	dest = (source + 1) % 2;
	step = 1;
	while (step < length)
	{
		m0 = buffer[source];
		end = buffer[source] + length;
		m1 = m0 + step;
		m2 = m1 + step;
		if (step * 2 < length)
		{
			k = buffer[dest];
		}
		else
		{
			k = order;
		}
		while (m0 < end)
		{
			if (m1 > end)
			{
				m1 = end;
			}
			if (m2 > end)
			{
				m2 = end;
			}
			i = m0;
			j = m1;
			while (i < m1 || j < m2)
			{
				if (j == m2)
				{
					*k = *i;
					++i;
				}
				else if (i == m1)
				{
					*k = *j;
					++j;
				}
				else
				{
					if (data[*j] < data[*i])
					{
						*k = *j;
						++j;
					}
					else
					{
						*k = *i;
						++i;
					}
				}
				++k;
			}
			m0 = m2;
			m1 = m0 + step;
			m2 = m1 + step;
		}
		step *= 2;
		source = dest;
		dest = (source + 1) % 2;
	}
	delete[]buffer[1];
	delete[]buffer[0];
}
template <typename T>
int sort_mat(Pool<MatElem<T> > &coef_mat, MatElem<T> **ordered_mat)
{
	MatElem<T> *buffer[2];
	int source, dest, step, len, mat_size;
	MatElem<T> *m0, *m1, *m2, *i, *j, *k, *mat_end;
	buffer[0] = new MatElem<T>[coef_mat.m_size];
	buffer[1] = new MatElem<T>[coef_mat.m_size];
	len = coef_mat.m_size;
	coef_mat.exact_data(buffer[0]);
	coef_mat.clear();
	source = 0;
	dest = (source + 1) % 2;
	step = 1;
	while (step < len)
	{
		m0 = buffer[source];
		m1 = m0 + step;
		m2 = m1 + step;
		mat_end = buffer[source] + len;
		k = buffer[dest];
		while (m0 < mat_end)
		{
			if (m1 > mat_end)
			{
				m1 = mat_end;
			}
			if (m2 > mat_end)
			{
				m2 = mat_end;
			}
			i = m0;
			j = m1;
			while ((i < m1) || (j < m2))
			{
				if (j == m2)
				{
					*k = *i;
					++i;
				}
				else if (i == m1)
				{
					*k = *j;
					++j;
				}
				else
				{
					if (i->col < j->col)
					{
						*k = *i;
						++i;
					}
					else if (j->col < i->col)
					{
						*k = *j;
						++j;
					}
					else
					{
						if (i->row < j->row)
						{
							*k = *i;
							++i;
						}
						else
						{
							*k = *j;
							++j;
						}
					}
				}
				++k;
			}
			m0 = m2;
			m1 = m0 + step;
			m2 = m1 + step;
		}
		step *= 2;
		source = dest;
		dest = (source + 1) % 2;
	}
	buffer[dest][0] = buffer[source][0];
	i = buffer[dest];
	j = buffer[source] + 1;
	mat_end = buffer[source] + len;
	while (j < mat_end)
	{
		if (j->col == i->col && j->row == i->row)
		{
			i->data += j->data;
		}
		else
		{
			++i;
			*i = *j;
		}
		++j;
	}
	mat_size = i - buffer[dest] + 1;
	*ordered_mat = new MatElem<T>[mat_size];
	mat_end = *ordered_mat + mat_size;
	i = buffer[dest];
	for (j = *ordered_mat; j < mat_end; ++j)
	{
		*j = *i;
		++i;
	}
	delete[]buffer[1];
	delete[]buffer[0];
	return mat_size;
}

template <typename T>
int sort_mat_row(Pool<MatElem<T> >& coef_mat, MatElem<T>** ordered_mat)
{
	MatElem<T>* buffer[2];
	int source, dest, step, len, mat_size;
	MatElem<T>* m0, * m1, * m2, * i, * j, * k, * mat_end;
	buffer[0] = new MatElem<T>[coef_mat.m_size];
	buffer[1] = new MatElem<T>[coef_mat.m_size];
	len = coef_mat.m_size;
	coef_mat.exact_data(buffer[0]);
	coef_mat.clear();
	source = 0;
	dest = (source + 1) % 2;
	step = 1;
	while (step < len)
	{
		m0 = buffer[source];
		m1 = m0 + step;
		m2 = m1 + step;
		mat_end = buffer[source] + len;
		k = buffer[dest];
		while (m0 < mat_end)
		{
			if (m1 > mat_end)
			{
				m1 = mat_end;
			}
			if (m2 > mat_end)
			{
				m2 = mat_end;
			}
			i = m0;
			j = m1;
			while ((i < m1) || (j < m2))
			{
				if (j == m2)
				{
					*k = *i;
					++i;
				}
				else if (i == m1)
				{
					*k = *j;
					++j;
				}
				else
				{
					if (i->row < j->row)
					{
						*k = *i;
						++i;
					}
					else if (j->row < i->row)
					{
						*k = *j;
						++j;
					}
					else
					{
						if (i->col < j->col)
						{
							*k = *i;
							++i;
						}
						else
						{
							*k = *j;
							++j;
						}
					}
				}
				++k;
			}
			m0 = m2;
			m1 = m0 + step;
			m2 = m1 + step;
		}
		step *= 2;
		source = dest;
		dest = (source + 1) % 2;
	}
	buffer[dest][0] = buffer[source][0];
	i = buffer[dest];
	j = buffer[source] + 1;
	mat_end = buffer[source] + len;
	while (j < mat_end)
	{
		if (j->col == i->col && j->row == i->row)
		{
			i->data += j->data;
		}
		else
		{
			++i;
			*i = *j;
		}
		++j;
	}
	mat_size = i - buffer[dest] + 1;
	*ordered_mat = new MatElem<T>[mat_size];
	mat_end = *ordered_mat + mat_size;
	i = buffer[dest];
	for (j = *ordered_mat; j < mat_end; ++j)
	{
		*j = *i;
		++i;
	}
	delete[]buffer[1];
	delete[]buffer[0];
	return mat_size;
}
#endif