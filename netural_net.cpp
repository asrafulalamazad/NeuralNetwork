#include <vector>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <fstream>
#include <sstream>

using namespace std;



class TrainingData{
  public:
    TrainingData(const string filename);
    bool isEof(void){ return m_trainingDataFile.eof(); }
    void getTopology(vector<unsigned> &topology);
    //return the number of input values read from the file;
    unsigned getNextInputs(vector<double> &inputVals);
    unsigned getTargetOutPuts(vector<double> &targetOutputVals);
  private:
    ifstream m_trainingDataFile;
};





struct Connection{
  double weight;
  double deltaWeight;
};



class Neuron;
typedef vector<Neuron> Layer;
// ----------------------------------------------
//              Neurom Class start
//-----------------------------------------------
class Neuron{
  public:
    Neuron(unsigned numOutputs,unsigned m_myindex );
    void setOutputVal(double val) { m_outputVal = val; }
    double getOutputVal(void) const { return m_outputVal ;}

    void feedForward(const Layer &prevLayer);

  private:
    static double treansferFunction(double x);
    static double treansferFunctionDerivative(double x);
    static double rendomWeight(void) { return  rand() / double(RAND_MAX); }
    double m_outputVal;
    vector<Connection> m_outputWeights;
    unsigned m_myindex;

};

double Neuron::treansferFunction(double x){
  return tanh(x);

}
double Neuron::treansferFunctionDerivative(double x){
  return 1.0 - x * x; 

}


void Neuron::feedForward(const Layer &prevLayer){
  double sum = 0.0;

  for(unsigned n=0;n < prevLayer.size() ; ++n){
    sum += prevLayer[n].getOutputVal() * prevLayer[n].m_outputWeights[m_myindex].weight;
  }
  m_outputVal = treansferFunction(sum);

}

Neuron::Neuron(unsigned numOutputs,unsigned m_myindex){
  for (unsigned c = 0;c < numOutputs; ++c){
    m_outputWeights.push_back(Connection());
    m_outputWeights.back().weight = rendomWeight();
  }

  m_myindex = m_myindex;

}

// ----------------------------------------------
//              Netueal_net Class start
//-----------------------------------------------
class Netural_net{
  private:
    vector<Layer> m_layer;
    double m_error;
    double m_recentAverageError;
    static double m_reacentAverageSmoothingFacetor;


  public:
    Netural_net(const vector<unsigned> &topology);

    void feedForward(const vector<double> &inputVals); 
    void backProp(const vector<double> &targetVals) {};
    void getResults(vector<double> &resultVals) const {};

};

double Netural_net::m_reacentAverageSmoothingFacetor = 100.0;

void Netural_net::getResults(vector<double> &resultVals) const {
  resultVals.clear();
  for (unsigned n=0;n<m_layer.back().size() -1 ;++n){
    resultVals.push_back(m_layer.back()[n].getOutputVal());
  }
}

void Netural_net::backProp(const vector<double> &targetVals){
  Layer &outputLayer = m_layer.back();
  m_error = 0.0;

  for (unsigned n=0;n < outputLayer.size() - 1; ++n){
    double delta = targetVals[n] - outputLayer[n].getOutputVal();
    m_error += delta *delta;
  }
  m_error /= outputLayer.size() - 1;
  m_error += sqrt(m_error);
  
  //Implement a recent avarage measurment
  m_recentAverageError = (m_recentAverageError * m_reacentAverageSmoothingFacetor + m_error) / (m_reacentAverageSmoothingFacetor + 1.0);
  for(unsigned n =0;n<outputLayer.size() - 1 ; ++n){
    outputLayer[n].clacOutputGradients(targetVals[n]);
  }

  for(unsigned layerNum = m_layer.size() - 2; layerNum > 0;--layerNum){
    Layer &hiddenLayer = m_layer[layerNum];
    Layer &nextLayer = m_layer[layerNum + 1];
    for(unsigned n=0;n< hiddenLayer.size();++n){
      hiddenLayer[n].calcHiddenGradients(nextLayer);
    }
  }

  //update connection weight
  for(unsigned layerNum = m_layer.size() - 1; layerNum > 0;--layerNum){
    Layer &layer = m_layer[layerNum];
    Layer &prevLayer = m_layer[layer - 1];

    for(unsigned n=0;n<layer.size() -1;++n){
      layer[n].updateInputWeights(prevLayer);
    }
  }


}

void Netural_net::feedForward(const vector<double> &inputVals){
  assert(inputVals.size() == m_layer[0].size() - 1 );
  //assign (latch) the input value into the input neuron
  for (unsigned i=0;i< inputVals.size() ; ++i){
    m_layer[0][i].setOutputVal(inputVals[i]);
  }

  //forward propagate
  for(unsigned layerNum = 1; m_layer.size(); ++layerNum){
    Layer &prevLayer = m_layer[layerNum - 1];
    for(unsigned n=0;n<m_layer[layerNum].size() - 1; ++n){
      m_layer[layerNum][n].feedForward(prevLayer);
    }
  }

}

Netural_net::Netural_net(const vector<unsigned> &topology){
  unsigned numLayers = topology.size();
  for (unsigned layerNum=0;layerNum < numLayers;++layerNum){
    m_layer.push_back(Layer());
    unsigned numOutputs = layerNum == topology.size() - 1 ? 0 : topology[layerNum + 1 ] ;
        //add new layer

    for (unsigned neuronNum = 0; neuronNum <= topology[layerNum]; ++neuronNum){
      m_layer.back().push_back(Neuron(numOutputs,neuronNum));
      cout << "Made a neuron !!" << endl;
    }
    m_layer.back().back().setOutputVal(1.0);
  }
}

void showVectorVals(string lable,vector<double> &v){
  cout << lable << " ";
  for(unsigned i=0;i<v.size();++i){
    cout << v[i] << " ";
  }
  cout << endl;
}


int main(int argc, char const *argv[]){

    TrainingData trainData("trainingData.txt");
    vector<unsigned> topology;
    // topology.push_back(3);
    // topology.push_back(2);
    // topology.push_back(1);

    

    Netural_net myNet(topology);

    vector<double> inputVals;
    myNet.feedForward(inputVals);

    vector<double> targetVals;
    myNet.backProp(targetVals);
    
    vector<double> resultVals;
    myNet.getResults(resultVals);

}
