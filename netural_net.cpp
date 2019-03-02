#include <vector>
#include <iostream>
using namespace std;

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
    Neuron();
  
  private:
    double m_outputVal;
    vector<Connection> m_outputWeights;

};



// ----------------------------------------------
//              Netueal_net Class start
//-----------------------------------------------
class Netural_net{
  private:
    vector<Layer> m_layer;

  public:
    Netural_net(const vector<unsigned> &topology);

    void feedForward(const vector<double> &inputVals) {};
    void backProp(const vector<double> &targetVals) {};
    void getResults(vector<double> &resultVals) const {};

};

Netural_net::Netural_net(const vector<unsigned> &topology){
  unsigned numLayers = topology.size();
  for (unsigned layerNum=0;layerNum < numLayers;++layerNum){
    m_layer.push_back(Layer());
    //add new layer

    for (unsigned neuronNum = 0;neuronNum <= topology[layerNum]; ++neuronNum){
      m_layer.back().push_back(Neuron());
      cout << "Made a neuron !!" << endl;
    }
  }
}




int main(int argc, char const *argv[]){
    vector<unsigned> topology;
    topology.push_back(3);
    topology.push_back(2);
    topology.push_back(1);

    Netural_net myNet(topology);

    vector<double> inputVals;
    myNet.feedForward(inputVals);

    vector<double> targetVals;
    myNet.backProp(targetVals);
    
    vector<double> resultVals;
    myNet.getResults(resultVals);

}
