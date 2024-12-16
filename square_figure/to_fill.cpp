#include <bits/stdc++.h>   
using namespace std; 
 
 
vector<int> split_int(string line, char delimeter){ 
    string my_str = ""; 
    vector<int> v; 
    for(int i; i<line.size(); i++){ 
        if((line[i] == delimeter) or (i==line.size()-1)){ 
            if(i==line.size()-1){ 
                my_str+=line[i];   
            } 
            char* temp = new char[my_str.length()+1]; 
            strcpy(temp, my_str.c_str()); 
            int to_insert = atof(temp); 
            v.push_back(to_insert); 
            my_str = ""; 
        }else{ 
            my_str+=line[i];   
        } 
    } 
    return v; 
} 
 
int main(){ 
    string myText; 
 
        //Nodes filling 
    ifstream NodesFile("node_to_cell.txt"); 
 
    vector<int> v;  
 
    while (getline(NodesFile, myText)) { 
        v = split_int(myText, ','); 
        for(int i=1; i<5;i++){ 
            cout<<v[0];
            cout<<",";
            cout<<i;
            cout<<",";
            cout<<v[i]<<endl; 
            // cout<<i<<endl; 
        } 
        v.clear(); 
    } 
 
    NodesFile.close(); 
    return 0; 
}