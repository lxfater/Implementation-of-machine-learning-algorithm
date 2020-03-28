#include "matrix.h"
#include <algorithm>
#include <math.h>
#include <string>
#include <fstream>
#include <sstream>
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock

using namespace std;

template <class T>
class Network
{
private:
    int layerNum;
    int *size;
    vector<Matrix<T>> biases;
    vector<Matrix<T>> weight;

public:
    Network(int layerNumber, int *neuronSize, T (*generator)());
    Matrix<T> feedforward(const Matrix<T> &matrix);
    void SGD(vector<vector<Matrix<T>>> trainingData, int epochs, T eta);
    void SGD(vector<vector<Matrix<T>>> trainingData, int epochs, T eta, vector<vector<Matrix<T>>> testData );
    void BGD(vector<vector<Matrix<T>>> trainingData, int epochs, T eta);
    void BGD(vector<vector<Matrix<T>>> trainingData, int epochs, T eta, vector<vector<Matrix<T>>> testData );
    void MBGD(vector<vector<Matrix<T>>> trainingData, int miniBatchSize, int epochs, T eta);
    void MBGD(vector<vector<Matrix<T>>> trainingData, int miniBatchSize, int epochs, T eta, vector<vector<Matrix<T>>> testData);
    void update(vector<vector<Matrix<T>>> trainingData, T eta);
    vector<T> colVectorMatrixToVector(Matrix<T> &x);
    int test(vector<vector<Matrix<T>>> testData);
    vector<vector<Matrix<T>>> backprop(Matrix<T> &x, Matrix<T> &y);
    static T sigmoidFuction(T x);
    static T sigmoidDerivativeFuction(T x);
    vector<Matrix<T>> getWeight();
    vector<Matrix<T>> getBiases();
    void print();
};
template <class T>
void Network<T>::print()
{
    vector<Matrix<double>> weight = this->getWeight();
    cout << "weight:" << endl;
    for (int i = 0; i < weight.size(); i++)
    {
        cout << "weight"<< i <<":" << endl;
        cout << weight[i];
    }
    cout << endl;

    vector<Matrix<double>> biases = this->getBiases();
    cout << "biases:" << endl;
    for (int i = 0; i < biases.size(); i++)
    {
        cout << "biases"<< i <<":" << endl;
        cout << biases[i];
    }
    cout << endl;
}
template <class T>
vector<Matrix<T>> Network<T>::getWeight()
{
    return weight;
}
template <class T>
vector<Matrix<T>> Network<T>::getBiases()
{
    return biases;
}
template <class T>
Network<T>::Network(int layerNumber, int *neuronSize, T (*generator)())
{
    layerNum = layerNumber;
    size = neuronSize;

    vector<int> dimensionCol(neuronSize, neuronSize + (layerNumber - 1));
    vector<int> dimensionRow(neuronSize + 1, neuronSize + layerNumber);

    for (int i = 0; i < layerNumber - 1; i++)
    {
        Matrix<T> tempWeight(dimensionRow[i], dimensionCol[i], generator);
        Matrix<T> tempBiases(dimensionRow[i], 1, generator);
        weight.push_back(tempWeight);
        biases.push_back(tempBiases);
    }
}

template <class T>
T Network<T>::sigmoidFuction(T x)
{
    return (1 / (1 + exp(-x)));
}

template <class T>
Matrix<T>
Network<T>::feedforward(const Matrix<T> &matrix)
{
    Matrix<T> temp = matrix;
    for (int i = 0; i < layerNum - 1; i++)
    {
        Matrix<T> r = weight[i] * temp + biases[i];
        temp = r.apply(sigmoidFuction);
    }
    return temp;
}
template <class T>
void Network<T>::SGD(vector<vector<Matrix<T>>> trainingData, int epochs, T eta)
{
    
    for (int i = 1; i < epochs+1; i++)
    {
        // 每回合前处理数据
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        shuffle(trainingData.begin(),trainingData.end(),default_random_engine(seed));
        for (int j = 0; j < trainingData.size(); j++)
        {
            vector<vector<Matrix<T>>> temp;
            temp.push_back(trainingData[j]);
            this->update(temp,eta);
        }
        cout << "Complete epoch:" << i << endl;
    }    
}
template <class T>
void Network<T>::SGD(vector<vector<Matrix<T>>> trainingData, int epochs, T eta, vector<vector<Matrix<T>>> testData)
{
    int biggest = 0;
    int epochIndex = -1;
    for (int i = 1; i < epochs+1; i++)
    {
        // 每回合前处理数据
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        shuffle(trainingData.begin(),trainingData.end(),default_random_engine(seed));
        int size = trainingData.size();
        for (int j = 0; j < size; j++)
        {
            vector<vector<Matrix<T>>> temp;
            temp.push_back(trainingData[j]);
            this->update(temp,eta);
        }
        int correct = test(testData);
        cout << "SGD epoch-" << i << " "<< correct <<"/"<<testData.size() << endl;
        if(biggest<correct){
            biggest = correct;
            epochIndex = i;
        }
    }
    cout << "SGD Best epoch-" << epochIndex << " "<< biggest <<"/"<<testData.size() << endl; 
}

template <class T>
void Network<T>::BGD(vector<vector<Matrix<T>>> trainingData, int epochs, T eta)
{

}
template <class T>
void Network<T>::BGD(vector<vector<Matrix<T>>> trainingData, int epochs, T eta, vector<vector<Matrix<T>>> testData)
{
    int biggest = 0;
    int epochIndex = -1;
    for (int i = 1; i < epochs+1; i++)
    {
        // 每回合前处理数据
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        shuffle(trainingData.begin(),trainingData.end(),default_random_engine(seed));
        this->update(trainingData,eta);
        int correct = test(testData);
        cout << "BGD epoch-" << i << " "<< correct <<"/"<<testData.size() << endl;
        // this->print();
        if(biggest<correct){
            biggest = correct;
            epochIndex = i;
        }
    } 
    cout << endl;
    cout << "BGD Best epoch-" << epochIndex << " "<< biggest <<"/"<<testData.size() << endl; 
}
template <class T>
void Network<T>::MBGD(vector<vector<Matrix<T>>> trainingData, int miniBatchSize, int epochs, T eta)
{
}
template <class T>
void Network<T>::MBGD(vector<vector<Matrix<T>>> trainingData, int miniBatchSize, int epochs, T eta, vector<vector<Matrix<T>>> testData)
{
}

template <class T>
void Network<T>::update(vector<vector<Matrix<T>>> trainingData, T eta)
{
    // 同意的维度
    int deminsion = layerNum - 1;
    // 创建梯度向量的数组,B为偏置向量，W为权重向量
    vector<Matrix<T>> gradientBVector;
    vector<Matrix<T>> gradientWVector;

    vector<int> dimensionCol(size, size + (layerNum - 1));
    vector<int> dimensionRow(size + 1, size + layerNum);

    for (int i = 0; i < layerNum - 1; i++)
    {
        Matrix<T> tempWeight(dimensionRow[i], dimensionCol[i]);
        Matrix<T> tempBiases(dimensionRow[i], 1);
        gradientWVector.push_back(tempWeight);
        gradientBVector.push_back(tempBiases);
    }

    //新权，新偏置,不停累加梯度
    for (int i = 0; i < trainingData.size(); i++)
    {
        vector<vector<Matrix<T>>> gradientVector = backprop(trainingData[i][0], trainingData[i][1]);
        for (int j = 0; j < deminsion; j++)
        {
            int a;
            gradientBVector[j] = gradientBVector[j] + gradientVector[0][j];
            gradientWVector[j] = gradientWVector[j] + gradientVector[1][j];
        }
    }
    for (int j = 0; j < deminsion; j++)
    {
        double p = eta / trainingData.size();
        biases[j] = biases[j] - (gradientBVector[j].multiply(p));
        weight[j] = weight[j] - (gradientWVector[j].multiply(p));
    }
}

template <class T>
vector<vector<Matrix<T>>> Network<T>::backprop(Matrix<T> &x, Matrix<T> &y)
{
    // 同意的维度
    int deminsion = layerNum - 1;
    // 创建梯度向量的数组,B为偏置向量，W为权重向量
    vector<Matrix<T>> gradientBVector(deminsion);
    vector<Matrix<T>> gradientWVector(deminsion);

    // 前向传播，记录中间值与激活值方便后面根据公式求导
    // 中间值
    vector<Matrix<T>> zVector;
    // 激活值
    vector<Matrix<T>> activationVector;
    // 输入为第一个激活值
    activationVector.push_back(x);

    Matrix<T> z;
    Matrix<T> activation = x;
    for (int i = 0; i < deminsion; i++)
    {
        z = weight[i] * activation + biases[i];
        zVector.push_back(z);
        activation = z.apply(sigmoidFuction);
        activationVector.push_back(activation);
    }

    // 反向传播

    // 激活值比中间值出一个，因为input也算
    Matrix<T> costDerivative = activationVector[deminsion] - y;
    Matrix<T> sigmoidDerivative = zVector[deminsion - 1].apply(sigmoidDerivativeFuction);
    Matrix<T> delta = costDerivative.Hadamard(sigmoidDerivative);
    gradientBVector[deminsion - 1] = delta;
    Matrix<T> transpose = activationVector[deminsion - 1].transpose();
    gradientWVector[deminsion - 1] = delta.multiply(transpose);

    for (int i = 2; i < layerNum; i++)
    {
        Matrix<T> z = zVector[deminsion - i];
        Matrix<T> sigmoidDerivative = z.apply(sigmoidDerivativeFuction);
        Matrix<T> weightTranspose = weight[deminsion - i + 1].transpose();
        delta = (weightTranspose * delta).Hadamard(sigmoidDerivative);

        gradientBVector[deminsion - i] = delta;
        Matrix<T> transpose = activationVector[deminsion - i].transpose();
        gradientWVector[deminsion - i] = delta * (transpose);
    }

    vector<vector<Matrix<T>>> gradientVector(2);

    gradientVector[0] = gradientBVector;
    gradientVector[1] = gradientWVector;

    return gradientVector;
}

// 测试多输出
template <class T> int
Network<T>::test(vector<vector<Matrix<T>>> testData)
{
        int sum =0;
        for (int j = 0; j < testData.size(); j++)
        {
            Matrix<T> result = this->feedforward(testData[j][0]);
            // result.matrix.
            vector<T> x = colVectorMatrixToVector(result);
            typename vector<T>::iterator xBiggest = max_element(begin(x),end(x));
            int xIndex = distance(begin(x), xBiggest);

            vector<T> y = colVectorMatrixToVector(testData[j][1]);
            typename vector<T>::iterator yBiggest = max_element(begin(y),end(y));
            int yIndex = distance(begin(y), yBiggest); 
            if(xIndex == yIndex){
                sum += 1;
            } 
        }
        return sum;        
}
// 测试当输出
// template <class T> int
// Network<T>::test(vector<vector<Matrix<T>>> testData)
// {
//         int sum =0;
//         for (int j = 0; j < testData.size(); j++)
//         {
//             Matrix<T> result = this->feedforward(testData[j][0]);
//             vector<T> x = colVectorMatrixToVector(result);
//             vector<T> y = colVectorMatrixToVector(testData[j][1]);
//             if(abs(y[0]-x[0])<=0.5){
//                 sum += 1;
//             } 
//         }
//         return sum;        
// }

template <class T> vector<T>
Network<T>::colVectorMatrixToVector(Matrix<T> &x)
{
    if(x.getColSize() !=1){
        throw "It is not colVectorMatrix";
    }
    vector<T> temp;
    for (int i = 0; i < x.getRowSize(); i++)
    {
        temp.push_back(x.getMatrix()[i][0]);
    }
    return temp;
    
}

template <class T>
T Network<T>::sigmoidDerivativeFuction(T x)
{
    return (sigmoidFuction(x) * (1 - sigmoidFuction(x)));
}

vector<double> rowGenerator(int col)
{
    vector<double> tmp(col, 1);
    return tmp;
}

double generator()
{
    // srand((unsigned)time(NULL));
    double a = rand()%10*0.1;
    return a;
}
// 加载iris
vector<vector<Matrix<double>>>
irisLoadData(string path)
{
    fstream f(path);
    string line;
    vector<vector<Matrix<double>>> tempVectorVectorMatrix;
    while (getline(f, line, '\n'))
    {
        istringstream iss(line);
        string temp;

        vector<Matrix<double>> tempVectorMatrix;

        Matrix<double> y;
        double arrayX[4][1];

        int index = 0;
        while (getline(iss, temp, ','))
        {
            if (temp == "Iris-setosa")
            {
                double arrayY[3][1] = {{1}, {0}, {0}};
                Matrix<double> tempMatrix(3, 1, (double *)arrayY);
                y = tempMatrix;
            }
            else if (temp == "Iris-versicolor")
            {
                double arrayY[3][1] = {{0}, {1}, {0}};
                Matrix<double> tempMatrix(3, 1, (double *)arrayY);
                y = tempMatrix;
            }
            else if (temp == "Iris-virginica")
            {
                double arrayY[3][1] = {{0}, {0}, {1}};
                Matrix<double> tempMatrix(3, 1, (double *)arrayY);
                y = tempMatrix;
            }
            else
            {
                arrayX[index][0] = atof(temp.c_str());
            }
            index++;
        }
        Matrix<double> x(4, 1, (double *)arrayX);

        tempVectorMatrix.push_back(x);
        tempVectorMatrix.push_back(y);

        tempVectorVectorMatrix.push_back(tempVectorMatrix);
    }
    f.close();
    return tempVectorVectorMatrix;
}
vector<vector<Matrix<double>>>
hiLoadData(string path)
{
    fstream f(path);
    string line;
    vector<vector<Matrix<double>>> tempVectorVectorMatrix;
    while (getline(f, line, '\n'))
    {
        istringstream iss(line);
        string temp;

        vector<Matrix<double>> tempVectorMatrix;

        Matrix<double> y;
        double arrayX[2][1];

        int index = 0;
        while (getline(iss, temp, '\t'))
        {
            if (index == 2){

                if( stoi(temp) == 0){
                  double arrayY[1][1] = {{0}};
                  Matrix<double> tempMatrix(1, 1, (double *)arrayY);
                  y = tempMatrix;
                }else{
                  double arrayY[1][1] = {{1}};
                  Matrix<double> tempMatrix(1, 1, (double *)arrayY);
                  y = tempMatrix;
                }

            }
            else
            {
                    arrayX[index][0] = atof(temp.c_str());    
            }
            index++;
            
        }
        Matrix<double> x(2, 1, (double *)arrayX);

        tempVectorMatrix.push_back(x);
        tempVectorMatrix.push_back(y);

        tempVectorVectorMatrix.push_back(tempVectorMatrix);
    }
    f.close();
    return tempVectorVectorMatrix;
}
int main()
{
    try
    {
        string path("iris.txt");
        vector<vector<Matrix<double>>> data = irisLoadData(path);
        unsigned seed = chrono::system_clock::now().time_since_epoch().count();
        shuffle(data.begin(),data.end(),default_random_engine(seed));
        const int cut = 100;
        vector<vector<Matrix<double>>> trainingData;
        vector<vector<Matrix<double>>> testData;
        for (size_t i = 0; i < cut; i++)
        {
            trainingData.push_back(data[i]);
        }
        for (size_t i = cut; i < data.size(); i++)
        {
            testData.push_back(data[i]);
        }
        
        int c[3] = {4, 20,3};
        Network<double> a(3, c, generator);

        a.SGD(trainingData,400,1.5,testData);
        // a.print();


        // double arrayY[2][1] = {{1}, {0}};
        // Matrix<double> y(2, 1, (double *)arrayY);
        // double arrayX[3][1] = {{1},{1},{3}};
        // Matrix<double> x(3, 1, (double *)arrayX);

        // int c[3] = {3, 2, 2};
        // Network<double> a(3, c, generator);
        // vector<vector<Matrix<double>>> gradientVector = a.backprop(x,y);
        // Matrix<double> sss = a.feedforward(x);
        return 0;
    }
    catch (char const *e)
    {
        cout << e << endl;
    }
}
