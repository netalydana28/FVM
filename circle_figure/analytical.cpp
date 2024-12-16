#include <bits/stdc++.h>  
using namespace std;
const int node_num=313;
map<int, vector<double> > nodes;
vector<double> split_double(string line, char delimeter){
    string my_str = "";
    vector<double> v;
    for(int i; i<line.size(); i++){
        if(((line[i] == delimeter) or (i==line.size()-1))){
            if(i==line.size()-1){
                my_str+=line[i];  
            }
            char* temp = new char[my_str.length()+1];
            strcpy(temp, my_str.c_str());
            double to_insert = atof(temp);
            v.push_back(to_insert);
            my_str = "";
        }else{
            my_str+=line[i];  
        }
    }
    return v;
}
void nodes_coor(){
    string myText;
    ifstream NodesFile("coordinates.txt");
    vector<double> v; 
    int cnt = 1;
    while (getline(NodesFile, myText)) {
        v = split_double(myText, ' ');
        // cout<<v[0]<<" "<<v[1]<<endl;
        v.erase(v.end()-1);
        nodes[cnt] = v;
        v.clear();
        cnt++;
    }
    NodesFile.close();
}

int main(){
    
    ofstream MyFile("pressure_analitycal.dat");
    MyFile << "VARIABLES= \"X\", \"Y\", \"P\"" <<endl;
    for(int i=1; i<node_num+1; i++){
        MyFile << nodes[i][0] << " " << nodes[i][1] <<" "<<sin(nodes[i][0] + 2 * nodes[i][1]) + pow(M_E, (2 * nodes[i][0] + 3 * nodes[i][1]))<<endl;
        // cout<<i+1<< " "<<P_faces[i]<<endl;
    }
    // for(int i=0; i<face_num; i++){
    //     MyFile <<P_faces[i]<<endl;
    //     // cout<<i+1<< " "<<P_faces[i]<<endl;
    // }
    MyFile.close();
    return 0;
}