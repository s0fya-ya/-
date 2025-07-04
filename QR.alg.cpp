// QR.alg.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

class Exception : public exception
{
protected:
	char* str;
public:
	Exception(const char* s)
	{
		str = new char[strlen(s) + 1];
		strcpy_s(str, strlen(s) + 1, s);
	}
	Exception(const Exception& e)
	{
		str = new char[strlen(e.str) + 1];
		strcpy_s(str, strlen(e.str) + 1, e.str);
	}
	~Exception()
	{
		delete[] str;
	}

	virtual void print()
	{
		cout << "Exception: " << str << "; " << what();
	}
};

class InvalidOperandException : public Exception
{
protected:
	int width;
	int height;
public:
	InvalidOperandException(const char* s, int Height, int Width) :
		Exception(s) {
		width = Width; height = Height;
	}
	virtual void print()
	{
		cout << "InvalidOperandException: " << str << "; width: " << width <<
			", height: " << height << "; " << what();
	}
};

class IndexOutOfBoundsException : public Exception
{
protected:
	int row;
	int column;
public:
	IndexOutOfBoundsException(const char* s, int Row, int Column) : Exception(s)
	{
		row = Row; column = Column;
	}
	virtual void print()
	{
		cout << "IndexOutOfBoundsException: " << str << "; row: " << row <<
			", column: " << column << "; " << what();
	}
};

class WrongSizeException : public Exception
{
protected:
	int width;
	int height;
public:
	WrongSizeException(const char* s, int Height, int Width) :Exception(s)
	{
		width = Width; height = Height;
	}
	virtual void print()
	{
		cout << "WrongSizeException: " << str << "; width: " << width <<
			", height: " << height << "; " << what();
	}
};

class NonPositiveException : public WrongSizeException
{
public:
	NonPositiveException(const char* s, int Height, int Width) :WrongSizeException(s, Height, Width) {}
	virtual void print()
	{
		cout << "NonPositiveException: " << str << "; width: " << width <<
			", height: " << height << "; " << what();
	}
};

class TooLargeSizeException : public WrongSizeException
{
public:
	TooLargeSizeException(const char* s, int Height, int Width) :WrongSizeException(s, Height, Width) {}
	virtual void print()
	{
		cout << "TooLargeSizeException: " << str << "; width: " << width <<
			", height: " << height << "; " << what();
	}
};

template <class T>
class Vector
{
protected:
	T* ptr;
	int length;

public:

	Vector<T>(int len = 2)
	{
		Vector<T>::length = len;
		Vector<T>::ptr = new T[len];
	}

	~Vector()
	{
		delete[] Vector<T>::ptr;
	}

	T getPtr(int index)
	{
		return ptr[index];
	}
	int getlength()
	{
		return length;
	}
	int setLength(int newlength)
	{
		length = newlength;
	}
	void setPtr(int index, T value)
	{
		ptr[index] = value;
	}

	Vector<T>(const Vector<T>& v)
	{
		Vector<T>::length = v.length;
		Vector<T>::ptr = new T[length];
		for (int i = 0; i < Vector<T>::length; i++)
			Vector<T>::ptr[i] = v.ptr[i];
	}

	T& operator[](int i) {
		return ptr[i];
	}
	T& operator[](int i) const {
		return ptr[i];
	}

	virtual void print()
	{
		for (int i = 0; i < length; i++)
			cout << ptr[i] << " ";
		cout << "\n";
	}
};

template<class T>
class Matrix : Vector<T>
{
protected:
	T** ptr;
	int height;
	int width;
public:

	Matrix <T>(int Height = 2, int Width = 2)
	{
		if (Height <= 0 || Width <= 0)
			throw NonPositiveException("Attempt to create matrix of negative or zero size", Height, Width);
		if (Height >= 20 || Width >= 20)
			throw TooLargeSizeException("Matrix is too large", Height, Width);
		Matrix<T>::height = Height;
		Matrix<T>::width = Width;
		Matrix<T>::ptr = new T * [Matrix<T>::height];
		for (int i = 0; i < Matrix<T>::height; i++)
			Matrix<T>::ptr[i] = new T[Matrix<T>::width];
	}

	Matrix <T>(const Matrix<T>& M)
	{
		Matrix<T>::width = M.width;
		Matrix<T>::height = M.height;
		Matrix<T>::ptr = new T * [Matrix<T>::height];
		for (int i = 0; i < Matrix<T>::height; i++)
		{
			Matrix<T>::ptr[i] = new T[Matrix<T>::width];
			for (int j = 0; j < Matrix<T>::width; j++)
				Matrix<T>::ptr[i][j] = M.ptr[i][j];
		}
	}

	~Matrix()
	{
		if (Matrix<T>::ptr != NULL)
		{
			for (int i = 0; i < Matrix<T>::height; i++)
				delete[] Matrix<T>::ptr[i];
			delete[] Matrix<T>::ptr;
			Matrix<T>::ptr = NULL;
		}
	}

	T getPtr(int Height, int Width)
	{
		return ptr[Height][Width];
	}
	int getHeight()
	{
		return height;
	}
	int getWidth()
	{
		return width;
	}
	int setHeight(int newheight)
	{
		height = newheight;
	}
	int setWidth(int newwidth)
	{
		width = newwidth;
	}
	void setPtr(int Height, int Width, T value)
	{
		ptr[Height][Width] = value;
	}

	Matrix<T>& operator=(const Matrix<T>& m)
	{
		Matrix<T>::height = m.height;
		Matrix<T>::width = m.width;
		delete[] Matrix<T>::ptr;
		Matrix<T>::ptr = new T * [Matrix<T>::height];
		for (int i = 0; i < Matrix<T>::height; i++)
			Matrix<T>::ptr[i] = new T[Matrix<T>::width];

		for (int i = 0; i < Matrix<T>::height; i++)
			for (int j = 0; j < Matrix<T>::width; j++)
				Matrix<T>::ptr[i][j] = m.ptr[i][j];
		return*this;
	}

	Matrix<T> operator+(const Matrix<T>& m)
	{
		if (Matrix<T>::height != m.height || Matrix<T>::width != m.width)
			throw InvalidOperandException("Width and height aren't equal in operator*", m.height, this->width);

		Matrix<T>::height = m.height;
		Matrix<T>::width = m.width;
		Matrix res(Matrix<T>::height, Matrix<T>::width);

		for (int i = 0; i < Matrix<T>::height; i++)
			for (int j = 0; j < Matrix<T>::width; j++)
				res[i][j] = Matrix<T>::ptr[i][j] + m.ptr[i][j];

		return res;
	}

	Matrix<T> operator-(const Matrix<T>& m)
	{
		if (Matrix<T>::height != m.height || Matrix<T>::width != m.width)
			throw InvalidOperandException("Width and height aren't equal in operator*", m.height, this->width);

		Matrix<T>::height = m.height;
		Matrix<T>::width = m.width;
		Matrix res(Matrix<T>::height, Matrix<T>::width);

		for (int i = 0; i < Matrix<T>::height; i++)
			for (int j = 0; j < Matrix<T>::width; j++)
				res[i][j] = Matrix<T>::ptr[i][j] - m.ptr[i][j];

		return res;
	}

	Matrix<T> operator*(const Matrix<T>& m)
	{
		if (Matrix<T>::width != m.height)
			throw InvalidOperandException("Width and height aren't equal in operator*", m.height, this->width);

		Matrix<T> res(Matrix<T>::height, m.width);

		for (int i = 0; i < Matrix<T>::height; i++)
			for (int j = 0; j < m.width; j++)
			{
				double sum = 0;
				for (int k = 0; k < Matrix<T>::width; k++)
					sum += Matrix<T>::ptr[i][k] * m.ptr[k][j];
				if (abs(sum) < 0.0000001)
					sum = 0;
				res.ptr[i][j] = sum;
			}
		return res;
	}

	virtual void print()
	{
		for (int i = 0; i < Matrix<T>::height; i++)
		{
			for (int j = 0; j < Matrix<T>::width; j++)
				cout << Matrix<T>::ptr[i][j] << " ";
			cout << "\n";
		}
	}

	T* operator[](int index)
	{
		if (index < 0 || index >= Matrix<T>::height)
			throw IndexOutOfBoundsException("Wrong index in operator[]", index, -1);
		return Matrix<T>::ptr[index];
	}
	T& operator()(int row, int column) const
	{
		if (row < 0 || row >= height || column < 0 || column >= width)
			throw IndexOutOfBoundsException("Wrong index in operator()", row, column);
		return Matrix<T>::ptr[row][column];
	}

	template<class T1>
	friend bool operator==(Matrix<T1>& m, Matrix<T1>& n);
	template<class T1>
	friend bool operator!=(Matrix<T1>& m, Matrix<T1>& n);

	Matrix<T> transpose()
	{
		Matrix n(Matrix<T>::width, Matrix<T>::height);

		for (int i = 0; i < Matrix<T>::height; i++)
			for (int j = 0; j < Matrix<T>::width; j++)
				n[j][i] = ptr[i][j];

		return n;
	}

	Matrix<T> QRDecomposition()
	{
		Matrix<T> Q(Matrix<T>::height, Matrix<T>::width);
		for (int j = 0; j < Matrix<T>::width; j++)
		{
			Matrix<T> v(Matrix<T>::height, 1);
			for (int i = 0; i < Matrix<T>::height; i++)
			{
				v(i, 0) = Matrix<T>::ptr[i][j];
			}
			for (int i = 0; i < j; i++)
			{
				T dot_product = 0;
				for (int k = 0; k < Matrix<T>::height; k++)
				{
					dot_product += Q(k, i) * Matrix<T>::ptr[k][j];
				}
				for (int k = 0; k < Matrix<T>::height; k++)
				{
					v(k, 0) -= dot_product * Q(k, i);
				}
			}
			T norm = 0;
			for (int i = 0; i < Matrix<T>::height; i++)
			{
				norm += v(i, 0) * v(i, 0);
			}
			norm = sqrt(norm);
			if (norm == 0)
				for (int i = 0; i < Matrix<T>::height; i++)
					Q(i, j) = 0;
			else
				for (int i = 0; i < Matrix<T>::height; i++)
					Q(i, j) = v(i, 0) / norm;
		}
		return Q;
	}

	Matrix<T> Identity(int n)
	{
		Matrix<T> result(n, n);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				result(i, j) = (i == j) ? 1 : 0;
			}
		}
		return result;
	}

	Matrix<T> eigenvalues(const double eps = 0.000001)
	{
		int n = Matrix<T>::height;
		//Matrix<T> Q = Identity(n);
		Matrix<T> R = *this;
		Matrix<T> Q1, R1;
		while (true)
		{
			// QR-разложение матрицы R
			Q1 = R.QRDecomposition();
			R1 = Q1.transpose() * R;
			//Q = Q * Q1;
			R = R1 * Q1;
			// Проверка на сходимость
			T sum = 0;
			for (int j = 0; j < n - 1; j++)
			{
				for (int i = j + 1; i < n; i++)
					sum += abs(R(i, j));
			}
			if (sum < eps)
				break;
		}
		Matrix<T> eigenvals(n, 1);
		for (int i = 0; i < n; i++)
			eigenvals(i, 0) = R(i, i);

		return eigenvals;
	}

	template <class T1>
	friend ostream& operator<<(ostream& s, Matrix<T1> M);
	template <class T1>
	friend istream& operator>>(istream& s, Matrix<T1>& M);

	Matrix<T> eigenvectors(Matrix<T> eigenvalues, const double eps = 1e-10)
	{
		int n = Matrix<T>::height;
		Matrix<T> eigenvectors = Identity(n);

		while (true)
		{
			double max_offdiag = 0.0;
			int p = 0, q = 0;

			for (int i = 0; i < n; i++)
				for (int j = i + 1; j < n; j++)
				{
					double offdiag = abs(Matrix<T>::ptr[i][j]);
					if (offdiag > max_offdiag)
					{
						max_offdiag = offdiag;
						p = i;
						q = j;
					}
				}

			if (max_offdiag < eps)
				break;

			double t = 0.5 * atan2(2 * Matrix<T>::ptr[p][q], eigenvalues.ptr[q][0] - eigenvalues.ptr[p][0]);
			double c = cos(t);
			double s = sin(t);

			double Apq = Matrix<T>::ptr[p][q];
			Matrix<T>::ptr[p][q] = 0.0;
			Matrix<T>::ptr[q][p] = 0.0;
			for (int k = 0; k < n; k++)
				if (k != p && k != q)
				{
					double Akp = Matrix<T>::ptr[k][p];
					double Akq = Matrix<T>::ptr[k][q];
					Matrix<T>::ptr[k][p] = c * Akp - s * Akq;
					Matrix<T>::ptr[p][k] = Matrix<T>::ptr[k][p];
					Matrix<T>::ptr[k][q] = c * Akq + s * Akp;
					Matrix<T>::ptr[q][k] = Matrix<T>::ptr[k][q];
				}

			double dp = eigenvalues.ptr[p][0];
			double dq = eigenvalues.ptr[q][0];
			eigenvalues.ptr[p][0] = c * c * dp - 2 * s * c * Apq + s * s * dq;
			eigenvalues.ptr[q][0] = s * s * dp + 2 * s * c * Apq + c * c * dq;

			for (int i = 0; i < n; i++)
			{
				double Vip = eigenvectors[i][p];
				double Viq = eigenvectors[i][q];
				eigenvectors[i][p] = c * Vip - s * Viq;
				eigenvectors[i][q] = s * Vip + c * Viq;
			}
			return eigenvectors;
		}
	}
	template <class T>
	friend void printEigenvaluesAndEigenvectors(const Matrix<T>& eigenvalues, const Matrix<T>& eigenvectors);
};

template<class T>
bool operator==(Matrix<T>& m, Matrix<T>& n)
{
	if (m.height != n.height || m.width != m.width)
		return false;
	for (int i = 0; i < m.height; i++)
		for (int j = 0; j < m.width; j++)
			if (m.ptr[i][j] != n.ptr[i][j])
				return false;

	return true;
}

template<class T>
bool operator!=(Matrix<T>& m, Matrix<T>& n)
{
	return !(m == n);
}

template <class T>
ostream& operator<< (ostream& s, Matrix<T> M)
{
	if (typeid(s) == typeid(ofstream))
	{
		//вывод в файл
		s << M.height << " " << M.width << '\n';
		for (int i = 0; i < M.height; i++)
		{
			for (int j = 0; j < M.width; j++)
				s << M.ptr[i][j] << " ";
			s << '\n';
		}
	}
	else
	{
		//вывод в консоль
		for (int i = 0; i < M.height; i++)
		{
			for (int j = 0; j < M.width; j++)
				s << M(i, j) << " ";
			s << "\n";
		}
	}
	return s;
}

template <class T>
istream& operator>>(istream& s, Matrix<T>& M)
{
	if (typeid(s) == typeid(ifstream))
	{
		int rows, columns;
		s >> rows >> columns;
		if (rows != M.height || columns != M.width)
		{
			for (int i = 0; i < M.height; i++)
				delete[] M.ptr[i];
			delete[] M.ptr;
			M.ptr = NULL;
			M.height = rows;
			M.width = columns;
			M.ptr = new T * [M.height];
			for (int i = 0; i < M.height; i++)
			{
				M.ptr[i] = new T[M.width];
				for (int j = 0; j < M.width; j++)
					s >> M.ptr[i][j];
			}
		}
		else
		{
			for (int i = 0; i < rows; i++)
				for (int j = 0; j < columns; j++)
					s >> M.ptr[i][j];
		}
	}
	else
	{
		for (int i = 0; i < M.height; i++)
			for (int j = 0; j < M.width; j++)
				s >> M.ptr[i][j];
	}
	return s;
}

template<class T>
Vector<T> operator*(Matrix<T>& m, Vector<T>& v)
{
	Vector<T> result(m.getHeight());

	for (int i = 0; i < m.getHeight(); i++)
	{
		T sum = 0;
		for (int j = 0; j < m.getWidth(); j++)
			sum += m.getPtr(i, j) * v.getPtr(j);

		result.setPtr(i, sum);
	}

	return result;
}

template<class T>
void printEigenvaluesAndEigenvectors(Matrix<T>& eigenvalues, Matrix<T>& eigenvectors)
{
	int numEigenvalues = eigenvalues.getHeight();

	int* indices = new int[numEigenvalues];
	for (int i = 0; i < numEigenvalues; i++)
		indices[i] = i;

	for (int i = 0; i < numEigenvalues - 1; i++)
		for (int j = 0; j < numEigenvalues - i - 1; j++)
			if (eigenvalues(indices[j], 0) > eigenvalues(indices[j + 1], 0))
			{
				int temp = indices[j];
				indices[j] = indices[j + 1];
				indices[j + 1] = temp;
			}

	for (int i = 0; i < numEigenvalues; i++)
	{
		int index = indices[i];
		if (i != numEigenvalues - 1 && (eigenvalues(indices[i], 0) == eigenvalues(indices[i + 1], 0)))
		{
			int k = 0;
			for (int j = 0; j < numEigenvalues; j++)
				if (eigenvectors(j, indices[i]) == eigenvectors(j, indices[i + 1]))
					k++;
			if (k != numEigenvalues)
			{
				cout << "Eigenvalue: " << eigenvalues(index, 0) << endl;
				cout << "Eigenvector: ";
				for (int j = 0; j < eigenvectors.getHeight(); j++)
					cout << eigenvectors(j, index) << " ";
				cout << endl << endl;
			}
		}
		else
		{
			cout << "Eigenvalue: " << eigenvalues(index, 0) << endl;
			cout << "Eigenvector: ";
			for (int j = 0; j < eigenvectors.getHeight(); j++)
				cout << eigenvectors(j, index) << " ";
			cout << endl << endl;
		}
	}
	delete[] indices;
}

int main()
{
	try
	{
		Matrix <double> T(3, 3);
		T(0, 0) = 1;
		T(0, 1) = 0;
		T(1, 0) = 0;
		T(1, 1) = 0;
		T(2, 0) = 0;
		T(2, 1) = 0;
		T(2, 2) = 0;
		T(1, 2) = 0;
		T(0, 2) = 0;

		cout << "\nMatrix: " << endl;
		T.print();
		Matrix<double> Q = T.QRDecomposition(); //матрица Q
		Matrix<double> R = Q.transpose() * T; //матрица R
		Matrix<double> v = T.eigenvalues(); //собственные числа
		Matrix<double> V1 = T.eigenvectors(v); //собственные вектора
		cout << "\n";
		printEigenvaluesAndEigenvectors(v, V1);
	}
	catch (InvalidOperandException e)
	{
		e.print();
	}
	catch (IndexOutOfBoundsException e)
	{
		e.print();
	}
	catch (WrongSizeException e)
	{
		e.print();
	}
	catch (NonPositiveException e)
	{
		e.print();
	}
	catch (TooLargeSizeException e)
	{
		e.print();
	}
	catch (Exception e)
	{
		e.print();
	}

	Matrix<double> A(2, 2);
	cout << "\nEnter elements of the matrix: " << endl;
	cin >> A;
	cout << "\nYour matrix: " << endl << A;

	Vector <double> vec1;
	vec1[0] = 3;
	vec1[1] = 8;
	cout << "\nVector: " << endl;
	vec1.print();

	Vector<double> vec = A * vec1;
	cout << "\nMatrix * Vector" << endl;
	vec.print();

	char c; cin >> c;
	return 0;
}
