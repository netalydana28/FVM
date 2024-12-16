//Navier Stokes on unstructured mesh 
#include <bits/stdc++.h>  
using namespace std;
const int node_num=4015, center_num=7828;
map<int, vector<double> > nodes, centers;
int face_num;
map<vector<int> , int> l_cell_to_face, l_cell_to_node, l_face_to_cell, l_face_to_node;
set<int> l_node_to_nodes[node_num], l_center_to_centers[center_num];
map<int, vector<int> > node_to_cell;
map<vector<int> , double> sn;
vector<int> boundary_face;
vector<int> nodes_boundary;
double P_nodes[node_num], U_nodes[node_num], V_nodes[node_num];
vector<double> P_faces, U_faces, V_faces;
double P_centers_old[center_num], P_centers_new[center_num];
double U_centers_old[center_num], U_centers_new[center_num];
double V_centers_old[center_num], V_centers_new[center_num];
double f[center_num];
double eps=0.0001, d = 1.0/200.0, dt = 0.000001, ro = 1.0, max_diff;
vector<double> split_double(string line, char delimeter){
    string my_str = "";
    vector<double> v;
    for(int i; i<line.size(); i++){
        if(((line[i] == delimeter) or (i==line.size()-1))){
            if(line[i+1] == delimeter){
                continue;
            }
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
vector<int> split_int(string line, char delimeter){
    string my_str = "";
    vector<int> v;
    for(int i; i<line.size(); i++){
        if((line[i] == delimeter) or (i==line.size()-1)){
            if(line[i+1] == delimeter){
                continue;
            }
            if(i==line.size()-1){
                my_str+=line[i];  
            }
            if(my_str[0]== 't'){
                my_str = "";
                continue;
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
void nodes_coor(){
    string myText;
    ifstream NodesFile("sq1_1_coor.txt");
    vector<double> v; 
    int cnt = 1;
    while (getline(NodesFile, myText)) {
        v = split_double(myText, ' ');
        // cnt = int(v[v.begin()]);
        v.erase(v.begin());
        v.erase(v.end()-1);
        nodes[cnt] = v;
        v.clear();
        cnt++;
    }
    NodesFile.close();
}
void cell_to_node(){
    string myText;
    ifstream NodesFile("sq1_1_links.txt");
    vector<int> v; 
    vector<int> v1; 
    while (getline(NodesFile, myText)) {
        v = split_int(myText, ' ');
        int cell = v[0];
        v.erase(v.begin());
        v.erase(v.begin());
        for(int i=1; i<4; i++){
            v1.push_back(cell);
            v1.push_back(i);
            l_cell_to_node[v1] = v[i-1];
            if(node_to_cell.find(v[i-1]) != node_to_cell.end()){
               node_to_cell[v[i-1]].push_back(cell);
            }else{
                vector<int> temp;
                temp.push_back(cell);
                node_to_cell[v[i-1]] = temp;
                temp.clear();
            }
            v1.clear();
        }
        v.clear();
    }
    NodesFile.close();
}
void center_coordinates_generator(){
    vector<int> v;
    vector<double> v1;
    for(int i = 1; i < center_num + 1; i++){
        v.push_back(i);
        v.push_back(1);
        register int node1 = l_cell_to_node[v];
        v.erase(v.end()-1);
        v.push_back(2);
        register int node2 = l_cell_to_node[v];
        v.erase(v.end()-1);
        v.push_back(3);
        register int node3 = l_cell_to_node[v];
        v1.push_back((nodes[node1][0] + nodes[node2][0] + nodes[node3][0]) / 3);
        v1.push_back((nodes[node1][1] + nodes[node2][1] + nodes[node3][1]) / 3);
        centers[i]=v1;
        v.clear();
        v1.clear();
    }
}
void face_links(){
    map< set<int>, vector<int> > nodes_to_cell;
    vector<int> key;
    vector<int> key1;
    register int cnt = 1;
    vector<int> trg_nodes;
    set<int> nodes_ind;
    vector< vector<int> > c_t_f(center_num);
    for(int i=1; i<center_num+1; i++){
        key.push_back(i);
        for(int j=1; j<4; j++){
            key.push_back(j);
            register int trg_node = l_cell_to_node[key];
            trg_nodes.push_back(trg_node);
            key.erase(key.end()-1);
        }
        key.clear();
        register int m;
        for(int j=1; j<4; j++){
            if(j==3){
                m = 0;
            }else{
                m = j;
            }
            nodes_ind.insert(trg_nodes[j-1]); 
            nodes_ind.insert(trg_nodes[m]);
            if(nodes_to_cell.find(nodes_ind) == nodes_to_cell.end()){
                key1.push_back(i);
                nodes_to_cell[nodes_ind] = key1; 
                key1.clear(); 
            }else{
                nodes_to_cell[nodes_ind].push_back(i);
            }
            nodes_ind.clear();           
        }
        trg_nodes.clear();
    }
    face_num = nodes_to_cell.size();
    map< set<int>, vector<int> >::iterator m_it=nodes_to_cell.begin();
    register int f_cnt = 1;
    while(m_it!=nodes_to_cell.end()){
        key1 = m_it->second; 
        set <int>::iterator s_it = m_it->first.begin();
        register int s_cnt = 1;
        int node1;
        int node2;
        while(s_it != m_it->first.end()){
            key.push_back(f_cnt); 
            key.push_back(s_cnt);
            l_face_to_node[key] = *s_it;
            if(s_cnt == 1){
                node1 = *s_it;
            }else if(s_cnt == 2){
                node2 = *s_it;
            }
            s_cnt++;
            s_it++;
            key.clear();
            
        }
        l_node_to_nodes[node1-1].insert(node2);
        l_node_to_nodes[node2-1].insert(node1);
        key.push_back(f_cnt);
        for(int i=1; i<key1.size()+1; i++){
            c_t_f[key1[i-1]-1].push_back(f_cnt);
            key.push_back(i);
            l_face_to_cell[key] = key1[i-1]; 
            key.erase(key.end()-1);
            if(key1.size() == 1){
                boundary_face.push_back(f_cnt);
                key.push_back(i+1);
                key.erase(key.end()-1);
                l_face_to_cell[key] = key1[i-1]; 
            }   
            l_center_to_centers[key1[0]].insert(key1[1]);
            l_center_to_centers[key1[1]].insert(key1[0]);
            
        }
        f_cnt++;
        m_it++;
        key1.clear();
        key.clear();
    }
    key.push_back(boundary_face[1]);
    key.push_back(1);
    key1.push_back(boundary_face[2]);
    key1.push_back(0);
    if(l_face_to_node[key] != l_face_to_node[key1]){
        int temp = boundary_face[0];
        boundary_face[0] = boundary_face[1];
        boundary_face[1] = temp;
    }
    key.clear();
    key1.clear();
    for(int i=1; i<center_num+1; i++){
        key.push_back(i);
        for(int j=1; j<c_t_f[i-1].size()+1; j++){
            key.push_back(j);
            l_cell_to_face[key] = c_t_f[i-1][j-1];
            key.erase(key.end()-1);
        }
        key.clear();
    }
    for(int i=0; i<boundary_face.size(); i++){
        int b_f = boundary_face[i];
        key.push_back(b_f); 
        key.push_back(1);
        int n1 = l_face_to_node[key];
        double x1 = nodes[l_face_to_node[key]][0];
        double y1 = nodes[l_face_to_node[key]][1]; 
        key.erase(key.end()-1);
        key.push_back(2);
        int n2 = l_face_to_node[key];
        double x2 = nodes[l_face_to_node[key]][0];
        double y2 = nodes[l_face_to_node[key]][1]; 
        nodes_boundary.push_back(n1);
        nodes_boundary.push_back(n2);
        key.clear();
    }
    for(int i=0; i<center_num; i++){
        f[i] = 0;
    }
}
vector<double> face_center(int face_ind){
    vector <int> key; 
    key.push_back(face_ind); 
    key.push_back(1);
    vector<double> n1 = nodes[l_face_to_node[key]];
    key.erase(key.end()-1);
    key.push_back(2);
    vector<double> n2 = nodes[l_face_to_node[key]];
    vector<double> ans; 
    ans.push_back((n1[0] + n2[0])/2);
    ans.push_back((n1[1] + n2[1])/2);
    return ans;
}
vector<double> normal(int face_ind){
    vector <double> ans; 
    vector <int> key; 
    key.push_back(face_ind); 
    key.push_back(1);
    vector<double> n1 = nodes[l_face_to_node[key]];
    key.erase(key.end()-1);
    key.push_back(2);
    vector<double> n2 = nodes[l_face_to_node[key]];
    ans.push_back(n1[1] - n2[1]); 
    ans.push_back(n2[0] - n1[0]);  
    double len = pow(pow(ans[0], 2) + pow(ans[1], 2), 0.5); 
    for(int i=0; i<ans.size(); i++){
        ans[i] /= len;
    }
    return ans; 
}
double area(int face_ind){
    vector <int> key; 
    key.push_back(face_ind); 
    key.push_back(1);
    vector<double> n1 = nodes[l_face_to_node[key]];
    key.erase(key.end()-1);
    key.push_back(2);
    vector<double> n2 = nodes[l_face_to_node[key]];
    return pow((pow((n1[0] - n2[0]), 2) + pow((n1[1] - n2[1]), 2)) , 0.5);
}
void cell_to_vertex_interpolation_old(){ 
    for(int node_ind=1; node_ind<node_num+1; node_ind++){
        vector<double> node_coor = nodes[node_ind]; 
        vector<double> center_coor;  
        int cell_ind; 
        double nom_U = 0, denom = 0; 
        double nom_V = 0; 
        double nom_P = 0; 
        for(int i=0; i<node_to_cell[node_ind].size(); i++){ 
            cell_ind = node_to_cell[node_ind][i]; 
            center_coor = centers[cell_ind]; 
            double d = pow(pow((center_coor[0] - node_coor[0]), 2) + pow((center_coor[1] - node_coor[1]), 2), 0.5); 
            double P = P_centers_old[cell_ind - 1];
            double U = U_centers_old[cell_ind - 1]; 
            double V = V_centers_old[cell_ind - 1]; 
            nom_P += P/d;
            nom_U += U/d; 
            nom_V += V/d; 
            denom += 1/d; 
        } 
        U_nodes[node_ind-1] = nom_U/denom;
        V_nodes[node_ind-1] = nom_V/denom;
        P_nodes[node_ind-1] = nom_P/denom;
    }
}
void cell_to_vertex_interpolation_new(){ 
    for(int node_ind=1; node_ind<node_num+1; node_ind++){
        vector<double> node_coor = nodes[node_ind]; 
        vector<double> center_coor;  
        int cell_ind; 
        double nom_U = 0, denom = 0; 
        double nom_V = 0; 
        double nom_P = 0; 
        for(int i=0; i<node_to_cell[node_ind].size(); i++){ 
            cell_ind = node_to_cell[node_ind][i]; 
            center_coor = centers[cell_ind]; 
            double d = pow(pow((center_coor[0] - node_coor[0]), 2) + pow((center_coor[1] - node_coor[1]), 2), 0.5); 
            double P = P_centers_new[cell_ind - 1];
            double U = U_centers_new[cell_ind - 1]; 
            double V = V_centers_new[cell_ind - 1]; 
            nom_P += P/d;
            nom_U += U/d; 
            nom_V += V/d; 
            denom += 1/d; 
        } 
        U_nodes[node_ind-1] = nom_U/denom;
        V_nodes[node_ind-1] = nom_V/denom;
        P_nodes[node_ind-1] = nom_P/denom;
    }
}
void set_boundary_conditions_U(){// add boundary derivative
    vector<int> v;
    vector<int> key; 
    for(int l=0; l<boundary_face.size(); l++){
        int b_f = boundary_face[l];
        v.push_back(b_f); 
        v.push_back(1);
        int n[2];
        n[0] = l_face_to_node[v];
        double x1 = nodes[n[0]][0];
        double y1 = nodes[n[0]][1]; 
        v.erase(v.end()-1);
        v.push_back(2);
        n[1] = l_face_to_node[v];
        double x2 = nodes[n[1]][0];
        double y2 = nodes[n[1]][1]; 
        if(((x1==0.0) and (x2==0.0)) or ((y1==0.0) and (y2==0.0)) or ((y1==1.0) and (y2==1.0))){
            U_nodes[n[0]-1] = 0;
            U_nodes[n[1]-1] = 0;
        }else if((x1==1.0) and (x2==1.0)){
            // cout<<"nodes "<<n[0]<<" "<<n[1]<<"\n";
            v.clear();
            for(int k = 0; k<2; k++){
                for(int i=0; i<node_to_cell[n[k]].size(); i++){
                    v.push_back(node_to_cell[n[k]][i]);
                    // cout<<"node# "<<n[k]<<" center: "<<node_to_cell[n[k]][i]<<"\n";
                    double center_coor[2]; 
                    center_coor[0] = centers[node_to_cell[n[k]][i]][0];
                    center_coor[1] = centers[node_to_cell[n[k]][i]][1];
                    double nom = 0, denom = 0, valid_denom;
                    for(int j=1; j<4; j++){
                        v.push_back(j);
                        // cout<<"key "<<v[0]<<" "<<v[1]<<"\n";
                        int face_ind = l_cell_to_face[v];
                        v.erase(v.end()-1);
                        // cout<<node_to_cell[n[k]][i]<<" "<<j<<" "<<face_ind<<" index"<<"\n";
                        key.push_back(face_ind);
                        key.push_back(1);
                        int fn1 = l_face_to_node[key];
                        key.erase(key.end()-1);
                        key.push_back(2);
                        int fn2 = l_face_to_node[key];
                        key.clear();
                        // cout<<"face nodes "<<fn1<<" "<<fn2<<"\n";
                        if(fn2 == n[k]){
                            fn2 = fn1; 
                        }else if ((fn2 != n[k]) and (fn1 != n[k])){
                            // cout<<"Faceno "<<face_ind<<"\n";
                            continue;
                        }
                        // cout<<"Faceyes "<<face_ind<<"\n";
                        vector<double> face_coor = face_center(face_ind);
                        double nx = face_coor[1] - center_coor[1];
                        double ny = center_coor[0] - face_coor[0];
                        nx /= pow(pow(nx, 2) + pow(ny, 2), 0.5);
                        ny /= pow(pow(nx, 2) + pow(ny, 2), 0.5);
                        double area = pow(pow(face_coor[0] - center_coor[0], 2) + pow(face_coor[1] - center_coor[1], 2), 0.5);
                        double dx = nodes[fn2][0] - nodes[n[k]][0];
                        double dy = nodes[fn2][1] - nodes[n[k]][1];
                        if(((dx>0) and (nx<0)) or ((dx<0) and (nx>0))){
                            nx*=(-1);
                        }
                        if(((dy>0) and (ny<0)) or ((dy<0) and (ny>0))){
                            ny*=(-1);
                        }
                        valid_denom=0;
                        if(dx != 0){
                            valid_denom += (nx/dx);
                        }
                        if(dy != 0){
                            valid_denom += (ny/dy);
                        }
                        nom += U_nodes[fn2-1] * area * valid_denom;
                        denom += area * valid_denom;
                        // cout<<" end "<<v.size()<<"\n";
                        // cout<<"Area "<<area<<" valid denom "<<valid_denom<<"\n";
                        // cout<<"FACE " <<j<<" "<<nom<<" "<<denom<<"\n";
                    }
                    U_nodes[n[k]-1] = nom/denom;
                    // cout<<" final "<<nom<<" "<<denom<<"\n";
                    v.clear();
                }
            }
        }
        v.clear();
        if((x1==0.0) and (x2==0.0) and (y1>=0.4) and (y1<=0.6) and (y2>=0.4) and (y2<=0.6)){
            U_nodes[n[0]-1] = 1;
            U_nodes[n[1]-1] = 1;
        }
        U_faces[b_f-1] = (U_nodes[n[0]-1] + U_nodes[n[1]-1])/2;
        // cout<< U_faces[b_f-1]<<"\n";
    }
}
void set_boundary_conditions_V(){// add boundary derivative
    vector<int> v;
    vector<int> key; 
    for(int l=0; l<boundary_face.size(); l++){
        int b_f = boundary_face[l];
        v.push_back(b_f); 
        v.push_back(1);
        int n[2];
        n[0] = l_face_to_node[v];
        double x1 = nodes[n[0]][0];
        double y1 = nodes[n[0]][1]; 
        v.erase(v.end()-1);
        v.push_back(2);
        n[1] = l_face_to_node[v];
        double x2 = nodes[n[1]][0];
        double y2 = nodes[n[1]][1]; 
        if(((x1==0.0) and (x2==0.0)) or ((y1==0.0) and (y2==0.0)) or ((y1==1.0) and (y2==1.0))){
            V_nodes[n[0]-1] = 0;
            V_nodes[n[1]-1] = 0;
        }else if((x1==1.0) and (x2==1.0)){
            // cout<<"nodes "<<n[0]<<" "<<n[1]<<"\n";
            v.clear();
            for(int k = 0; k<2; k++){
                for(int i=0; i<node_to_cell[n[k]].size(); i++){
                    v.push_back(node_to_cell[n[k]][i]);
                    // cout<<"node# "<<n[k]<<" center: "<<node_to_cell[n[k]][i]<<"\n";
                    double center_coor[2]; 
                    center_coor[0] = centers[node_to_cell[n[k]][i]][0];
                    center_coor[1] = centers[node_to_cell[n[k]][i]][1];
                    double nom = 0, denom = 0, valid_denom;
                    for(int j=1; j<4; j++){
                        v.push_back(j);
                        // cout<<"key "<<v[0]<<" "<<v[1]<<"\n";
                        int face_ind = l_cell_to_face[v];
                        v.erase(v.end()-1);
                        // cout<<node_to_cell[n[k]][i]<<" "<<j<<" "<<face_ind<<" index"<<"\n";
                        key.push_back(face_ind);
                        key.push_back(1);
                        int fn1 = l_face_to_node[key];
                        key.erase(key.end()-1);
                        key.push_back(2);
                        int fn2 = l_face_to_node[key];
                        key.clear();
                        // cout<<"face nodes "<<fn1<<" "<<fn2<<"\n";
                        if(fn2 == n[k]){
                            fn2 = fn1; 
                        }else if ((fn2 != n[k]) and (fn1 != n[k])){
                            // cout<<"Faceno "<<face_ind<<"\n";
                            continue;
                        }
                        // cout<<"Faceyes "<<face_ind<<"\n";
                        vector<double> face_coor = face_center(face_ind);
                        double nx = face_coor[1] - center_coor[1];
                        double ny = center_coor[0] - face_coor[0];
                        nx /= pow(pow(nx, 2) + pow(ny, 2), 0.5);
                        ny /= pow(pow(nx, 2) + pow(ny, 2), 0.5);
                        double area = pow(pow(face_coor[0] - center_coor[0], 2) + pow(face_coor[1] - center_coor[1], 2), 0.5);
                        double dx = nodes[fn2][0] - nodes[n[k]][0];
                        double dy = nodes[fn2][1] - nodes[n[k]][1];
                        if(((dx>0) and (nx<0)) or ((dx<0) and (nx>0))){
                            nx*=(-1);
                        }
                        if(((dy>0) and (ny<0)) or ((dy<0) and (ny>0))){
                            ny*=(-1);
                        }
                        valid_denom=0;
                        if(dx != 0){
                            valid_denom += (nx/dx);
                        }
                        if(dy != 0){
                            valid_denom += (ny/dy);
                        }
                        nom += V_nodes[fn2-1] * area * valid_denom;
                        denom += area * valid_denom;
                        // cout<<" end "<<v.size()<<"\n";
                        // cout<<"Area "<<area<<" valid denom "<<valid_denom<<"\n";
                        // cout<<"FACE " <<j<<" "<<nom<<" "<<denom<<"\n";
                    }
                    V_nodes[n[k]-1] = nom/denom;
                    // cout<<" final "<<nom<<" "<<denom<<"\n";
                    v.clear();
                }
            }
        }
        v.clear();
        V_faces[b_f-1] = (V_nodes[n[0]-1] + V_nodes[n[1]-1])/2;
        // cout<< V_faces[b_f-1]<<"\n";
    }
}
void set_boundary_conditions_P(){// add boundary derivative
    vector<int> v;
    vector<int> key; 
    for(int l=0; l<boundary_face.size(); l++){
        int b_f = boundary_face[l];
        v.push_back(b_f); 
        v.push_back(1);
        int n[2];
        n[0] = l_face_to_node[v];
        double x1 = nodes[n[0]][0];
        double y1 = nodes[n[0]][1]; 
        v.erase(v.end()-1);
        v.push_back(2);
        n[1] = l_face_to_node[v];
        double x2 = nodes[n[1]][0];
        double y2 = nodes[n[1]][1]; 
        if((x1==0.0) and (x2==0.0)){
            P_nodes[n[0]-1] = 0;
            P_nodes[n[1]-1] = 0;
        }else if(((x1==1.0) and (x2==1.0)) or ((y1==1.0) and (y2==1.0)) or ((y1==0.0) and (y2==0.0))){
            // cout<<"nodes "<<n[0]<<" "<<n[1]<<"\n";
            v.clear();
            for(int k = 0; k<2; k++){
                for(int i=0; i<node_to_cell[n[k]].size(); i++){
                    v.push_back(node_to_cell[n[k]][i]);
                    // cout<<"node# "<<n[k]<<" center: "<<node_to_cell[n[k]][i]<<"\n";
                    double center_coor[2]; 
                    center_coor[0] = centers[node_to_cell[n[k]][i]][0];
                    center_coor[1] = centers[node_to_cell[n[k]][i]][1];
                    double nom = 0, denom = 0, valid_denom;
                    for(int j=1; j<4; j++){
                        v.push_back(j);
                        // cout<<"key "<<v[0]<<" "<<v[1]<<"\n";
                        int face_ind = l_cell_to_face[v];
                        v.erase(v.end()-1);
                        // cout<<node_to_cell[n[k]][i]<<" "<<j<<" "<<face_ind<<" index"<<"\n";
                        key.push_back(face_ind);
                        key.push_back(1);
                        int fn1 = l_face_to_node[key];
                        key.erase(key.end()-1);
                        key.push_back(2);
                        int fn2 = l_face_to_node[key];
                        key.clear();
                        // cout<<"face nodes "<<fn1<<" "<<fn2<<"\n";
                        if(fn2 == n[k]){
                            fn2 = fn1; 
                        }else if ((fn2 != n[k]) and (fn1 != n[k])){
                            // cout<<"Faceno "<<face_ind<<"\n";
                            continue;
                        }
                        // cout<<"Faceyes "<<face_ind<<"\n";

                        vector<double> face_coor = face_center(face_ind);
                        double nx = face_coor[1] - center_coor[1];
                        double ny = center_coor[0] - face_coor[0];
                        nx /= pow(pow(nx, 2) + pow(ny, 2), 0.5);
                        ny /= pow(pow(nx, 2) + pow(ny, 2), 0.5);
                        double area = pow(pow(face_coor[0] - center_coor[0], 2) + pow(face_coor[1] - center_coor[1], 2), 0.5);
                        double dx = nodes[fn2][0] - nodes[n[k]][0];
                        double dy = nodes[fn2][1] - nodes[n[k]][1];
                        if(((dx>0) and (nx<0)) or ((dx<0) and (nx>0))){
                            nx*=(-1);
                        }
                        if(((dy>0) and (ny<0)) or ((dy<0) and (ny>0))){
                            ny*=(-1);
                        }
                        valid_denom=0;
                        if(dx != 0){
                            valid_denom += (nx/dx);
                        }
                        if(dy != 0){
                            valid_denom += (ny/dy);
                        }
                        nom += P_nodes[fn2-1] * area * valid_denom;
                        denom += area * valid_denom;
                        // cout<<" end "<<v.size()<<"\n";
                        // cout<<"Area "<<area<<" valid denom "<<valid_denom<<"\n";
                        // cout<<"FACE " <<j<<" "<<nom<<" "<<denom<<"\n";
                    }
                    P_nodes[n[k]-1] = nom/denom;
                    // cout<<" final "<<nom<<" "<<denom<<"\n";
                    v.clear();
                }
            }
        }
        v.clear();
        P_faces[b_f-1] = (P_nodes[n[0]-1] + P_nodes[n[1]-1])/2;
        // cout<< P_faces[b_f-1]<<"\n";
    }
}
void fvm(){ //add eps and check cell to vertex interpolation
    for(int main=0; main < 1000; main++){
        cout<<main<<"\n";
        // cout<<"8.1"<<"\n";
        for(int i=1; i<center_num+1; i++){
            vector<int> key;
            key.push_back(i);
            double au=0;
            double av=0;
            vector<double> center_coor = centers[i];
            for(int j=1; j<4; j++){
                key.push_back(j);
                register int face_ind = l_cell_to_face[key];
                vector<int> key1;
                key1.push_back(face_ind);
                double deltx, delty, nx, ny, valid_denom;
                if(find(boundary_face.begin(), boundary_face.end(), face_ind) != boundary_face.end()){   
                    vector<double> face_coor = face_center(face_ind);
                    deltx = face_coor[0] - center_coor[0];
                    delty = face_coor[1] - center_coor[1];
                    key1.push_back(1);
                    nx = sn[key1];
                    key1.erase(key1.end()-1);
                    key1.push_back(2);
                    ny = sn[key1];
                    if(((deltx>0) and (nx<0)) or ((deltx<0) and (nx>0))){
                        nx*=(-1);
                    }
                    if(((delty>0) and (ny<0)) or ((delty<0) and (ny>0))){
                        ny*=(-1);
                    }
                    valid_denom=0;
                    if(deltx != 0){
                        valid_denom += (nx/deltx);
                    }
                    if(delty != 0){
                        valid_denom += (ny/delty);
                    }
                    // cout<<U_faces[face_ind-1]<<" "<<V_faces[face_ind-1]<<"\n";
                    au += (area(face_ind))*(U_faces[face_ind-1] - U_centers_old[i-1])*((d*valid_denom) - U_centers_old[i-1]*nx - V_centers_old[i-1]*ny);
                    av += (area(face_ind))*(V_faces[face_ind-1] - V_centers_old[i-1])*((d*valid_denom) - U_centers_old[i-1]*nx - V_centers_old[i-1]*ny);
                    key1.clear();
                }else{
                    key1.push_back(1);
                    nx = sn[key1];
                    register int out_center_ind = l_face_to_cell[key1]; 
                    key1.erase(key1.end()-1);
                    key1.push_back(2);
                    vector<double> c2 = centers[out_center_ind];
                    ny = sn[key1];
                    if(out_center_ind == i){
                        out_center_ind = l_face_to_cell[key1];
                        c2 = centers[out_center_ind];
                    }
                    deltx = c2[0] - center_coor[0];
                    delty = c2[1] - center_coor[1];
                    if(((deltx>0) and (nx<0)) or ((deltx<0) and (nx>0))){
                        nx*=(-1);
                    }
                    if(((delty>0) and (ny<0)) or ((delty<0) and (ny>0))){
                        ny*=(-1);
                    }
                    valid_denom=0;
                    if(deltx != 0){
                        valid_denom += (nx/deltx);
                    }
                    if(delty != 0){
                        valid_denom += (ny/delty);
                    }
                    au += ((U_centers_old[out_center_ind-1] - U_centers_old[i-1]) * area(face_ind))*((d*valid_denom) - U_centers_old[i-1]*nx - V_centers_old[i-1]*ny);
                    av += ((V_centers_old[out_center_ind-1] - V_centers_old[i-1]) * area(face_ind))*((d*valid_denom) - U_centers_old[i-1]*nx - V_centers_old[i-1]*ny);
                    key1.clear();
                }
                key.erase(key.end()-1);
            }
            U_centers_new[i-1] = U_centers_old[i-1] + dt * au;
            V_centers_new[i-1] = V_centers_old[i-1] + dt * av; 
            key.clear();
        }
        // cout<<"8.2"<<"\n";
        ofstream yFile1("u_centers_new.dat");
        yFile1 << "VARIABLES= \"X\", \"Y\", \"P\", \"U\", \"L\"" <<"\n";
        for(int i = 0; i<center_num; i++){
            yFile1 << centers[i+1][0] << " " << centers[i+1][1] <<" "<<U_centers_new[i]<<" "<<V_centers_new[i]<<"\n";
        }
        yFile1.close();
        cell_to_vertex_interpolation_new();
        set_boundary_conditions_U();
        set_boundary_conditions_V();
        // cout<<"8.3"<<"\n";
        for(int i=1; i<center_num+1; i++){
            vector<int> key;
            key.push_back(i);
            double au=0;
            double av=0;
            vector<double> center_coor = centers[i];
            for(int j=1; j<4; j++){
                key.push_back(j);
                register int face_ind = l_cell_to_face[key];
                vector<int> key1;
                key1.push_back(face_ind);
                double deltx, delty, nx, ny;
                if(find(boundary_face.begin(), boundary_face.end(), face_ind) != boundary_face.end()){   
                    vector<double> face_coor = face_center(face_ind);
                    deltx = face_coor[0] - center_coor[0];
                    delty = face_coor[1] - center_coor[1];
                    key1.push_back(1);
                    nx = sn[key1];
                    key1.erase(key1.end()-1);
                    key1.push_back(2);
                    ny = sn[key1];
                    if(((deltx>0) and (nx<0)) or ((deltx<0) and (nx>0))){
                        nx*=(-1);
                    }
                    if(((delty>0) and (ny<0)) or ((delty<0) and (ny>0))){
                        ny*=(-1);
                    }
                    au += (area(face_ind))*(U_faces[face_ind-1] - U_centers_new[i-1])*nx;
                    av += (area(face_ind))*(V_faces[face_ind-1] - V_centers_new[i-1])*ny;
                    key1.clear();
                }else{
                    key1.push_back(1);
                    nx = sn[key1];
                    register int out_center_ind = l_face_to_cell[key1]; 
                    key1.erase(key1.end()-1);
                    key1.push_back(2);
                    vector<double> c2 = centers[out_center_ind];
                    ny = sn[key1];
                    if(out_center_ind == i){
                        out_center_ind = l_face_to_cell[key1];
                        c2 = centers[out_center_ind];
                    }
                    deltx = c2[0] - center_coor[0];
                    delty = c2[1] - center_coor[1];
                    if(((deltx>0) and (nx<0)) or ((deltx<0) and (nx>0))){
                        nx*=(-1);
                    }
                    if(((delty>0) and (ny<0)) or ((delty<0) and (ny>0))){
                        ny*=(-1);
                    }
                    au += ((U_centers_new[out_center_ind-1] - U_centers_new[i-1]) * area(face_ind))*nx;
                    av += (V_centers_new[out_center_ind-1] - V_centers_new[i-1]) * area(face_ind)*ny;
                    key1.clear();
                }
                key.erase(key.end()-1);
            }
            f[i-1] = (ro/dt) * (au + av);
            au = 0;
            av = 0; 
            key.clear();
        }
        // cout<<"8.4"<<"\n";
        do{
            cell_to_vertex_interpolation_old();
            set_boundary_conditions_P();
            for(int i=1; i<center_num+1; i++){
                vector<int> key;
                key.push_back(i);
                double nom=0, denom=0;
                vector<double> center_coor = centers[i];
                for(int j=1; j<4; j++){
                    key.push_back(j);
                    register int face_ind = l_cell_to_face[key];
                    vector<int> key1;
                    key1.push_back(face_ind);
                    double deltx, delty, nx, ny, valid_denom;
                    if(find(boundary_face.begin(), boundary_face.end(), face_ind) != boundary_face.end()){   
                        vector<double> face_coor = face_center(face_ind);
                        deltx = face_coor[0] - center_coor[0];
                        delty = face_coor[1] - center_coor[1];
                        key1.push_back(1);
                        nx = sn[key1];
                        key1.erase(key1.end()-1);
                        key1.push_back(2);
                        ny = sn[key1];
                        if(((deltx>0) and (nx<0)) or ((deltx<0) and (nx>0))){
                            nx*=(-1);
                        }
                        if(((delty>0) and (ny<0)) or ((delty<0) and (ny>0))){
                            ny*=(-1);
                        }
                        valid_denom=0;
                        if(deltx != 0){
                            valid_denom += (nx/deltx);
                        }
                        if(delty != 0){
                            valid_denom += (ny/delty);
                        }
                        nom += (area(face_ind))*P_faces[face_ind-1]*valid_denom;
                        denom += (area(face_ind))*valid_denom;
                        key1.clear();
                    }else{
                        key1.push_back(1);
                        nx = sn[key1];
                        register int out_center_ind = l_face_to_cell[key1]; 
                        key1.erase(key1.end()-1);
                        key1.push_back(2);
                        vector<double> c2 = centers[out_center_ind];
                        ny = sn[key1];
                        if(out_center_ind == i){
                            out_center_ind = l_face_to_cell[key1];
                            c2 = centers[out_center_ind];
                        }
                        deltx = c2[0] - center_coor[0];
                        delty = c2[1] - center_coor[1];
                        if(((deltx>0) and (nx<0)) or ((deltx<0) and (nx>0))){
                            nx*=(-1);
                        }
                        if(((delty>0) and (ny<0)) or ((delty<0) and (ny>0))){
                            ny*=(-1);
                        }
                        valid_denom=0;
                        if(deltx != 0){
                            valid_denom += (nx/deltx);
                        }
                        if(delty != 0){
                            valid_denom += (ny/delty);
                        }
                        nom = (area(face_ind))*valid_denom*P_centers_old[out_center_ind-1];
                        denom += (area(face_ind))*valid_denom;
                        key1.clear();
                    }
                    key.erase(key.end()-1);
                }
                nom -= f[i-1];
                P_centers_new[i-1] = nom/denom;
                key.clear();
            }
            max_diff = 0.0; 
            for(int i = 0; i< center_num; i++){
                if(max_diff<fabs(P_centers_new[i] - P_centers_old[i])){
                    max_diff = fabs(P_centers_new[i] - P_centers_old[i]);
                }
                P_centers_old[i] = P_centers_new[i];
            }
            cell_to_vertex_interpolation_old();
            set_boundary_conditions_P();
        }while(max_diff>eps);
        // cout<<"8.5"<<"\n";
        for(int i=1; i<center_num+1; i++){
            vector<int> key;
            key.push_back(i);
            double au=0, av=0;
            vector<double> center_coor = centers[i];
            for(int j=1; j<4; j++){
                key.push_back(j);
                register int face_ind = l_cell_to_face[key];
                vector<int> key1;
                key1.push_back(face_ind);
                double deltx, delty, nx, ny, valid_denom;
                if(find(boundary_face.begin(), boundary_face.end(), face_ind) != boundary_face.end()){   
                    vector<double> face_coor = face_center(face_ind);
                    deltx = face_coor[0] - center_coor[0];
                    delty = face_coor[1] - center_coor[1];
                    key1.push_back(1);
                    nx = sn[key1];
                    key1.erase(key1.end()-1);
                    key1.push_back(2);
                    ny = sn[key1];
                    if(((deltx>0) and (nx<0)) or ((deltx<0) and (nx>0))){
                        nx*=(-1);
                    }
                    if(((delty>0) and (ny<0)) or ((delty<0) and (ny>0))){
                        ny*=(-1);
                    }
                    valid_denom=0;
                    if(deltx != 0){
                        valid_denom += (nx/deltx);
                    }
                    if(delty != 0){
                        valid_denom += (ny/delty);
                    }
                    au += (area(face_ind))*(P_faces[face_ind-1] - P_centers_new[i-1])*nx;
                    av += (area(face_ind))*(P_faces[face_ind-1] - P_centers_new[i-1])*ny;
                    key1.clear();
                }else{
                    key1.push_back(1);
                    nx = sn[key1];
                    register int out_center_ind = l_face_to_cell[key1]; 
                    key1.erase(key1.end()-1);
                    key1.push_back(2);
                    vector<double> c2 = centers[out_center_ind];
                    ny = sn[key1];
                    if(out_center_ind == i){
                        out_center_ind = l_face_to_cell[key1];
                        c2 = centers[out_center_ind];
                    }
                    deltx = c2[0] - center_coor[0];
                    delty = c2[1] - center_coor[1];
                    if(((deltx>0) and (nx<0)) or ((deltx<0) and (nx>0))){
                        nx*=(-1);
                    }
                    if(((delty>0) and (ny<0)) or ((delty<0) and (ny>0))){
                        ny*=(-1);
                    }
                    valid_denom=0;
                    if(deltx != 0){
                        valid_denom += (nx/deltx);
                    }
                    if(delty != 0){
                        valid_denom += (ny/delty);
                    }
                    au += (area(face_ind))*(P_centers_new[out_center_ind-1] - P_centers_new[i-1])*nx;
                    av += (area(face_ind))*(P_centers_new[out_center_ind-1] - P_centers_new[i-1])*ny;
                    key1.clear();
                }
                key.erase(key.end()-1);
            }
            U_centers_old[i-1] = U_centers_new[i-1] - (dt/ro) * au;
            V_centers_old[i-1] = V_centers_new[i-1] - (dt/ro) * av;
            key.clear();
        }
        // cout<<"8.6"<<"\n";
        for(int i = 0; i< center_num; i++){
            P_centers_old[i] = P_centers_new[i];
        }
        cell_to_vertex_interpolation_old();
        set_boundary_conditions_P();
        set_boundary_conditions_U();
        set_boundary_conditions_V();
    }
    cell_to_vertex_interpolation_old();
    set_boundary_conditions_P();
    set_boundary_conditions_U();
    set_boundary_conditions_V();
    // cout<<"8.7"<<"\n";
    ofstream yFile("out6.dat");
    yFile << "VARIABLES= \"X\", \"Y\", \"P\", \"U\", \"V\"" <<"\n";
    for(int i = 0; i<node_num; i++){
        yFile << nodes[i+1][0] << " " << nodes[i+1][1] <<" "<<P_nodes[i]<<" "<<U_nodes[i]<<" "<<V_nodes[i]<<"\n";
    }
    yFile.close();
    ofstream yFilec("out_centers.dat");
    yFilec << "VARIABLES= \"X\", \"Y\", \"P\", \"U\", \"V\"" <<"\n";
    for(int i = 0; i<center_num; i++){
        yFilec << centers[i+1][0] << " " << centers[i+1][1] <<" "<<P_centers_old[i]<<" "<<U_centers_old[i]<<" "<<V_centers_old[i]<<"\n";
    }
    yFilec.close();
}
int main(){// correct function calling order
    nodes_coor();//collecting coordinates and links
    cout<<"1"<<"\n";
    cell_to_node();
    cout<<"2"<<"\n";
    center_coordinates_generator();
    cout<<"3"<<"\n";
    face_links();
    cout<<"4"<<"\n";
    vector<double> v1;
    vector<int> v1_;
    cout<<"Boundary face num "<<boundary_face.size()<<"\n";
    for(int i=1; i<face_num+1; i++){
        v1 = normal(i);
        v1_.push_back(i);
        v1_.push_back(1);
        sn[v1_] = v1[0];
        v1_.erase(v1_.end()-1);
        v1_.push_back(2);
        sn[v1_] = v1[1];
        v1_.clear();
        v1.clear();
    }
    for(int i=0; i<center_num; i++){
        U_centers_new[i]=0;
        U_centers_old[i]=0;
    }
    for(int i=0; i<node_num; i++){
        U_nodes[i] = 0;
    }
    for(int i=0; i<face_num; i++){
        U_faces.push_back(0);
    } 
    for(int i=0; i<center_num; i++){
        P_centers_new[i]=0;
        P_centers_old[i]=0;
    } 
    for(int i=0; i<node_num; i++){
        P_nodes[i] = 0;
    }
    for(int i=0; i<face_num; i++){
        P_faces.push_back(0);
    } 
    for(int i=0; i<center_num; i++){
        V_centers_new[i]=0;
        V_centers_old[i]=0;
    } 
    for(int i=0; i<node_num; i++){
        V_nodes[i] = 0;
    }
    for(int i=0; i<face_num; i++){
        V_faces.push_back(0);
    } 
    set_boundary_conditions_U();
    cout<<"5"<<"\n";
    set_boundary_conditions_V();
    cout<<"6"<<"\n";
    set_boundary_conditions_P();
    cout<<"7"<<"\n";
    fvm();
    cout<<"8"<<"\n";
    return 0;
}