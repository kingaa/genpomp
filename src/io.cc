#include "io.h"
#include "omp.h"

using namespace std;

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

double StringToFloat_Type ( const string &Text )
{                               
	stringstream ss(Text);
	double result;
	return ss >> result ? result : 0;
}

void readParams( map<string,double> &paramMap, char *fileName){

  string line, parName, parVal;
  vector<string> vessel;  
  map<string,double> params;

  ifstream myfile (fileName);
  if (myfile.is_open())
    {
      while ( getline (myfile,line) )
      {
	vessel = split(line,' ');
	paramMap[vessel.at(0)] = StringToFloat_Type(vessel.at(1));
      }
      myfile.close();
    } else {
    cout << "error: Failed to read in parameters" << endl;
    cout << fileName << endl;
  }
}

// Function to read in random walk parameters, transformations for each parameter, and how to perturb each parameter
void read_random_walk_params( map<string,double> &paramMap, vector<string> &transforms, vector<string> &perturbations, char *fileName){

  string line,parName,parVal;
  vector<string> vessel;  
  map<string,double> params;

  ifstream myfile (fileName);
  if (myfile.is_open())
    {
      while ( getline (myfile,line) )
      {
	vessel = split(line,' ');
	paramMap[vessel.at(0)] = StringToFloat_Type(vessel.at(1));
	transforms.push_back(vessel.at(2));
	perturbations.push_back(vessel.at(3));
      }
      myfile.close();
    } else {
    cout << "error: Failed to read in parameters" << endl;
    cout << fileName << endl;
  }
}


//Function to read in sequences from a text file
void readSeqs( vector<double> &times, vector<string> &seqs, const char *seqFile) {
  string line,time,seq;
  vector<string> vessel;
  ifstream myfile (seqFile);
  int i = 0;
  if (myfile.is_open())
  {
    while ( getline (myfile,line) )
    {
      vessel = split(line,' ');
      if(vessel.at(0) != "times"){ // burn header if it exists
	times.push_back(StringToFloat_Type(vessel.at(0)));
	seqs.push_back(vessel.at(1));
      } else {
	cout << "Burning header: " << vessel.at(0) << " ... " << endl;
      }
    }
    myfile.close();
  } else {
    cout << "error: Failure to open in sequence datafile" << endl;
  }
}

//Function to write execution start and end times to file
void printStatus(ofstream &outfile, bool start, string functionName, int iter){
  time_t now = time(0);
  outfile << functionName << ' ' << iter << ' ' << start << ' ' << now << ' ' <<  omp_get_thread_num() << endl;
}
