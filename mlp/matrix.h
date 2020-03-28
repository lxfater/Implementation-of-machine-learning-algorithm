#ifndef _MATRIX_H
#define _MATRIX_H
#include <ostream>
#include <vector>
#include <istream>
#include<iostream>

using namespace std;

template <class T>
class Matrix
{

public:
    Matrix();
    Matrix(int rows, int columns);
    Matrix(int row, int col, T (*generator)());
    Matrix(int row, int col, T *array);
    Matrix<T> operator+(const Matrix<T> &matrix);
    Matrix<T> operator-(const Matrix<T> &matrix);
    Matrix<T> operator*(Matrix<T> &matrix);
    Matrix<T> multiply(Matrix<T> &matrix);
    Matrix<T> multiply(T num);
    Matrix<T> Hadamard(const Matrix<T> &matrix);
    Matrix<T> transpose();
    Matrix<T> apply(T (*func)(T));

    vector<T> getCol(int colIndex);
    vector<T> getRow(int colIndex);
    int getRowSize();
    int getColSize();
    void matrixInfo();
    vector<vector<T>> getMatrix();

private:
    int rowSize;
    int colSize;
    vector<vector<T>> matrix;
};

// 重载输出符号
template <class T>
inline ostream& operator<<(ostream &os, Matrix<T> &matrix){

    vector<vector<T>> temp = matrix.getMatrix();
    for (class vector<vector<T>>::iterator ite = temp.begin(); ite != temp.end(); ite++)
    {

        vector<T> temp_vect = *ite;
        for(class vector<T>::iterator itee = temp_vect.begin(); itee != temp_vect.end(); itee++)
             os << *itee << " ";
         os << endl;
    }
    return os;
}

template <class T>
Matrix<T>::Matrix()
{
}
template <class T>
Matrix<T>::Matrix(int row, int col)
{
    rowSize = row;
    colSize = col;

    for (int i = 0; i < row; i++)
    {
        vector<T> tmp(col);
        matrix.push_back(tmp);
    }
}

template <class T>
Matrix<T>::Matrix(int row, int col, T (*generator)())
{
    rowSize = row;
    colSize = col;

    for (int i = 0; i < row; i++)
    {
        vector<T> tmp;
        for (int i = 0; i < col; i++)
        {
            T num = generator();
            tmp.push_back(num);
        }
        
        matrix.push_back(tmp);
    }
}

template <class T>
Matrix<T>::Matrix(int row, int col, T *array)
{
    rowSize = row;
    colSize = col;

    for (int i = 0; i < row; i++)
    {
        vector<T> tmp;
        for (int j = 0; j < col; j++)
        {
            T num = *(array+i*col+j);
            tmp.push_back(num);
        }
        
        matrix.push_back(tmp);
    }
}

// 辅助函数
template <class T> void 
Matrix<T>::matrixInfo()
{
    cout << "rowSize: " << rowSize << " colSize: " << colSize << endl;
    for (class vector<vector<T>>::iterator ite = matrix.begin(); ite != matrix.end(); ite++)
    {
        vector<T> temp_vect = *ite;
        for(class vector<T>::iterator itee = temp_vect.begin(); itee != temp_vect.end(); itee++)
             cout << *itee << " ";
         cout << endl;
    }
}

template <class T> int
Matrix<T>::getColSize(){
    return colSize;
}
template <class T> int 
Matrix<T>::getRowSize(){
    return rowSize;
}
template <class T> vector<vector<T>> 
Matrix<T>::getMatrix(){
    return matrix;
}

template <class T> vector<T>
Matrix<T>::getCol(int colIndex){
    if(colIndex<0 || colIndex > colSize){
        throw "colIndex out of range";
    }
    vector<T> temp;
    for (int i = 0; i < rowSize; i++)
    {
        temp.push_back(matrix[i][colIndex]);
    }
    return temp;

}

template <class T> vector<T>
Matrix<T>::getRow(int rowIndex){
    if(rowIndex<0 || rowIndex > rowSize){
        throw "rowIndex out of range";
    }
    vector<T> temp;
    for (int i = 0; i < colSize; i++)
    {
        temp.push_back(matrix[rowIndex][i]);
    }
    return temp;

}

// 计算
template <class T> Matrix<T> 
Matrix<T>::operator+(const Matrix<T> &matrix){
    if((this->colSize != matrix.colSize) || (this->rowSize != matrix.rowSize)){
        throw "col or row do not match +";
    }
    
    Matrix newMatrix(this->rowSize,this->colSize);

    for (int i = 0; i < this->rowSize; i++)
    {
        for (int j = 0; j < this->colSize; j++)
        {
            
            newMatrix.matrix[i][j] = this->matrix[i][j] + matrix.matrix[i][j];
        }
        
    }

    return newMatrix;
}

template <class T> Matrix<T> 
Matrix<T>::operator*(Matrix &matrix){
    if(this->colSize != matrix.rowSize){
        throw "do not match *";
    }
    
    Matrix newMatrix(this->rowSize,matrix.colSize);

    for (int i = 0; i < this->rowSize; i++)
    {
        vector<T> row = this->getRow(i);
        for (int j = 0; j < matrix.colSize; j++)
        {
            
            vector<T> col = matrix.getCol(j);
            T sum = 0;
            for (int k = 0; k < this->colSize; k++)
            {
                sum += row[k]*col[k];
            }
            
            newMatrix.matrix[i][j] = sum;
        }
        
    }

    return newMatrix;
}

template <class T> Matrix<T> 
Matrix<T>::multiply(Matrix<T> &matrix){
 
    if(this->colSize != matrix.rowSize){
        throw "do not match *matrix";
    }
    
    Matrix<T> newMatrix(this->rowSize,matrix.colSize);

    for (int i = 0; i < this->rowSize; i++)
    {
        vector<T> row = this->getRow(i);
        for (int j = 0; j < matrix.colSize; j++)
        {
            
            vector<T> col = matrix.getCol(j);
            T sum = 0;
            for (int k = 0; k < this->colSize; k++)
            {
                sum += row[k]*col[k];
            }
            
            newMatrix.matrix[i][j] = sum;
        }
        
    }

    return newMatrix;
}

template <class T> Matrix<T> 
Matrix<T>::multiply(T num){
 
    Matrix<T> newMatrix(this->rowSize,this->colSize);

    for (int i = 0; i < this->rowSize; i++)
    {
        for (int j = 0; j < this->colSize; j++)
        {
            
            newMatrix.matrix[i][j] = this->matrix[i][j]*num;
        }
        
    }

    return newMatrix;
}

template <class T> Matrix<T> 
Matrix<T>::apply(T (*func)(T)){
 
    Matrix<T> newMatrix(this->rowSize,this->colSize);

    for (int i = 0; i < this->rowSize; i++)
    {
        for (int j = 0; j < this->colSize; j++)
        {
            
            newMatrix.matrix[i][j] = func(this->matrix[i][j]);
        }
        
    }

    return newMatrix;
}

template <class T> Matrix<T> 
Matrix<T>::operator-(const Matrix<T> &matrix){
    if((this->colSize != matrix.colSize) || (this->rowSize != matrix.rowSize)){
        throw "col or row do not match -";
    }
    
    Matrix<T> newMatrix(this->rowSize,this->colSize);

    for (int i = 0; i < this->rowSize; i++)
    {

        for (int j = 0; j < this->colSize; j++)
        {

            newMatrix.matrix[i][j] = this->matrix[i][j] - matrix.matrix[i][j];
        }
    }

    return newMatrix;
}

template <class T> Matrix<T> 
Matrix<T>::Hadamard(const Matrix<T> &matrix){

    if((this->colSize != matrix.colSize) || (this->rowSize != matrix.rowSize)){
        throw "col or row do not match h";
    }
    
    Matrix<T> newMatrix(this->rowSize,this->colSize);

    for (int i = 0; i < this->rowSize; i++)
    {

        for (int j = 0; j < this->colSize; j++)
        {

            newMatrix.matrix[i][j] = this->matrix[i][j] * matrix.matrix[i][j];
        }
    }

    return newMatrix;
}

template <class T> Matrix<T> 
Matrix<T>::transpose(){

    Matrix<T> newMatrix(this->colSize,this->rowSize);

    for (int i = 0; i < this->rowSize; i++)
    {

        for (int j = 0; j < this->colSize; j++)
        {

            newMatrix.matrix[j][i] = this->matrix[i][j];
        }
    }

    return newMatrix;
}
#endif