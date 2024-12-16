//Navier Stokes on unstructured mesh optimization 
//решение в ущерб памяти, в пользу времени
#include <bits/stdc++.h>  
using namespace std;
const int node_num=4015, center_num=7828;
double nodes[node_num][2], centers[center_num][2];
int face_num;
int l_cell_to_face[center_num][3];
int l_cell_to_node[center_num][3];
map<vector<int>, int> l_face_to_cell, l_face_to_node;
set<int> l_node_to_nodes[node_num], l_center_to_centers[center_num];
vector<int> node_to_cell[node_num];
vector<vector<double> > sn;
vector<vector<double> > f_centers;
vector<double> face_areas;
vector<int> boundary_face;
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
    for(int i=0; i<line.size(); i++){
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
    for(int i=0; i<line.size(); i++){
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
    int cnt = 0;
    while (getline(NodesFile, myText)) {
        v = split_double(myText, ' ');
        // cnt = int(v[v.begin()]);
        v.erase(v.begin());
        v.erase(v.end()-1);
        nodes[cnt][0] = v[0];
        nodes[cnt][1] = v[1];
        v.clear();
        cnt++;
    }
    NodesFile.close();
}
void cell_to_node(){
    string myText;
    ifstream NodesFile("sq1_1_links.txt");
    vector<int> v; 
    // vector<int> v1; 
    while (getline(NodesFile, myText)) {
        v = split_int(myText, ' ');
        int cell = v[0];
        v.erase(v.begin());
        v.erase(v.begin());
        for(int i=0; i<3; i++){
            // v1.push_back(cell);
            // v1.push_back(i);
            l_cell_to_node[cell-1][i] = v[i];
            node_to_cell[v[i]-1].push_back(cell);
            // if(node_to_cell.find(v[i]) != node_to_cell.end()){
            //    node_to_cell[v[i]].push_back(cell);
            // }else{
            //     vector<int> temp;
            //     temp.push_back(cell);
            //     node_to_cell[v[i]] = temp;
            //     temp.clear();
            // }
            // v1.clear();
        }
        v.clear();
    }
    NodesFile.close();
}
void center_coordinates_generator(){
    for(int i = 0; i < center_num; i++){
        register int node1 = l_cell_to_node[i][0];
        register int node2 = l_cell_to_node[i][1];
        register int node3 = l_cell_to_node[i][2];
        centers[i][0] = ((nodes[node1-1][0] + nodes[node2-1][0] + nodes[node3-1][0]) / 3);
        centers[i][1] = ((nodes[node1-1][1] + nodes[node2-1][1] + nodes[node3-1][1]) / 3);
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
        // key.push_back(i);
        for(int j=1; j<4; j++){
            // key.push_back(j);
            register int trg_node = l_cell_to_node[i-1][j-1];
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
        register int s_cnt = 1, node1, node2;
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
        // cout<<f_cnt<<endl;
        for(int i=1; i<key1.size()+1; i++){
            c_t_f[key1[i-1]-1].push_back(f_cnt);
            key.push_back(i);
            l_face_to_cell[key] = key1[i-1]; 
            key.erase(key.end()-1);
            if(key1.size() == 1){
                boundary_face.push_back(f_cnt);
                key.push_back(i+1);
                l_face_to_cell[key] = key1[i-1]; 
                key.erase(key.end()-1);
                continue;
                // cout<<"SJHB "<<key1[i-1]<<"\n";
            }   
            // l_center_to_centers[key1[0]].insert(key1[1]);
            // // cout<<f_cnt<<endl;
            // l_center_to_centers[key1[1]].insert(key1[0]); 
            // cout<<key1[1]<<endl;
            // cout<<f_cnt<<endl;
        }
        f_cnt++;
        m_it++;
        key1.clear();
        key.clear();
    }
    // key.push_back(boundary_face[1]);
    // key.push_back(1);
    // key1.push_back(boundary_face[2]);
    // key1.push_back(0);
    // if(l_face_to_node[key] != l_face_to_node[key1]){
    //     int temp = boundary_face[0];
    //     boundary_face[0] = boundary_face[1];
    //     boundary_face[1] = temp;
    // }
    // key.clear();
    // key1.clear();
    for(int i=0; i<center_num; i++){
        // key.push_back(i);
        for(int j=0; j<c_t_f[i].size(); j++){
            // key.push_back(j);
            l_cell_to_face[i][j] = c_t_f[i][j];
            key.erase(key.end()-1);
        }
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
    double n1[2];
    n1[0] = nodes[l_face_to_node[key]-1][0];
    n1[1] = nodes[l_face_to_node[key]-1][1];
    key.erase(key.end()-1);
    key.push_back(2);
    double n2[2];
    n2[0] = nodes[l_face_to_node[key]-1][0];
    n2[1] = nodes[l_face_to_node[key]-1][1];    vector<double> ans; 
    ans.push_back((n1[0] + n2[0])/2);
    ans.push_back((n1[1] + n2[1])/2);
    return ans;
}
vector<double> normal(int face_ind){
    vector <double> ans; 
    vector <int> key; 
    key.push_back(face_ind); 
    key.push_back(1);
    double n1[2];
    n1[0] = nodes[l_face_to_node[key]-1][0];
    n1[1] = nodes[l_face_to_node[key]-1][1];
    key.erase(key.end()-1);
    key.push_back(2);
    double n2[2];
    n2[0] = nodes[l_face_to_node[key]-1][0];
    n2[1] = nodes[l_face_to_node[key]-1][1];
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
    double n1[2];
    n1[0] = nodes[l_face_to_node[key]-1][0];
    n1[1] = nodes[l_face_to_node[key]-1][1];
    key.erase(key.end()-1);
    key.push_back(2);
    double n2[2];
    n2[0] = nodes[l_face_to_node[key]-1][0];
    n2[1] = nodes[l_face_to_node[key]-1][1];
    return pow((pow((n1[0] - n2[0]), 2) + pow((n1[1] - n2[1]), 2)) , 0.5);
}
double volume(int cen_ind){
    int n1 = l_cell_to_node[cen_ind-1][0];
    int n2 = l_cell_to_node[cen_ind-1][1];
    int n3 = l_cell_to_node[cen_ind-1][2];
    double x1 = nodes[n1-1][0];
    double x2 = nodes[n2-1][0];
    double x3 = nodes[n3-1][0];
    double y1 = nodes[n1-1][1];
    double y2 = nodes[n2-1][1];
    double y3 = nodes[n3-1][1];
    return (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) / 2;
}
void cell_to_vertex_interpolation_old(){ 
    for(int node_ind=1; node_ind<node_num+1; node_ind++){
        double node_coor[2];
        node_coor[0]= nodes[node_ind-1][0]; 
        node_coor[1]= nodes[node_ind-1][1];  
        double center_coor[2];  
        int cell_ind; 
        double nom_U = 0, denom = 0; 
        double nom_V = 0; 
        double nom_P = 0; 
        for(int i=0; i<node_to_cell[node_ind-1].size(); i++){ 
            cell_ind = node_to_cell[node_ind-1][i]; 
            center_coor[0] = centers[cell_ind-1][0]; 
            center_coor[1] = centers[cell_ind-1][1]; 
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
        double node_coor[2];
        node_coor[0]= nodes[node_ind-1][0]; 
        node_coor[1]= nodes[node_ind-1][1]; 
        double center_coor[2];  
        int cell_ind; 
        double nom_U = 0, denom = 0; 
        double nom_V = 0; 
        double nom_P = 0; 
        for(int i=0; i<node_to_cell[node_ind-1].size(); i++){ 
            cell_ind = node_to_cell[node_ind-1][i]; 
            center_coor[0] = centers[cell_ind-1][0]; 
            center_coor[1] = centers[cell_ind-1][1]; 
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
        double x1 = nodes[n[0]-1][0];
        double y1 = nodes[n[0]-1][1]; 
        v.erase(v.end()-1);
        v.push_back(2);
        n[1] = l_face_to_node[v];
        double x2 = nodes[n[1]-1][0];
        double y2 = nodes[n[1]-1][1]; 
        if(((x1==0.0) and (x2==0.0)) or ((y1==0.0) and (y2==0.0)) or ((y1==1.0) and (y2==1.0))){
            U_nodes[n[0]-1] = 0;
            U_nodes[n[1]-1] = 0;
        }else if((x1==1.0) and (x2==1.0)){
            // cout<<"nodes "<<n[0]<<" "<<n[1]<<"\n";
            v.clear();
            for(int k = 0; k<2; k++){
                for(int i=0; i<node_to_cell[n[k]-1].size(); i++){
                    v.push_back(node_to_cell[n[k]-1][i]);
                    // cout<<"node# "<<n[k]<<" center: "<<node_to_cell[n[k]][i]<<"\n";
                    double center_coor[2]; 
                    center_coor[0] = centers[node_to_cell[n[k]-1][i]-1][0];
                    center_coor[1] = centers[node_to_cell[n[k]-1][i]-1][1];
                    double nom = 0, denom = 0, valid_denom;
                    for(int j=1; j<4; j++){
                        v.push_back(j);
                        // cout<<"key "<<v[0]<<" "<<v[1]<<"\n";
                        int face_ind = l_cell_to_face[node_to_cell[n[k]-1][i]-1][j-1];
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
                        vector<double> face_coor = f_centers[face_ind-1];
                        double nx = face_coor[1] - center_coor[1];
                        double ny = center_coor[0] - face_coor[0];
                        nx /= pow(pow(nx, 2) + pow(ny, 2), 0.5);
                        ny /= pow(pow(nx, 2) + pow(ny, 2), 0.5);
                        double area = pow(pow(face_coor[0] - center_coor[0], 2) + pow(face_coor[1] - center_coor[1], 2), 0.5);
                        double dx = nodes[fn2-1][0] - nodes[n[k]-1][0];
                        double dy = nodes[fn2-1][1] - nodes[n[k]-1][1];
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
        double x1 = nodes[n[0]-1][0];
        double y1 = nodes[n[0]-1][1]; 
        v.erase(v.end()-1);
        v.push_back(2);
        n[1] = l_face_to_node[v];
        double x2 = nodes[n[1]-1][0];
        double y2 = nodes[n[1]-1][1]; 
        if(((x1==0.0) and (x2==0.0)) or ((y1==0.0) and (y2==0.0)) or ((y1==1.0) and (y2==1.0))){
            V_nodes[n[0]-1] = 0;
            V_nodes[n[1]-1] = 0;
        }else if((x1==1.0) and (x2==1.0)){
            // cout<<"nodes "<<n[0]<<" "<<n[1]<<"\n";
            v.clear();
            for(int k = 0; k<2; k++){
                for(int i=0; i<node_to_cell[n[k]-1].size(); i++){
                    v.push_back(node_to_cell[n[k]-1][i]);
                    // cout<<"node# "<<n[k]<<" center: "<<node_to_cell[n[k]][i]<<"\n";
                    double center_coor[2]; 
                    center_coor[0] = centers[node_to_cell[n[k]-1][i]-1][0];
                    center_coor[1] = centers[node_to_cell[n[k]-1][i]-1][1];
                    double nom = 0, denom = 0, valid_denom;
                    for(int j=1; j<4; j++){
                        v.push_back(j);
                        // cout<<"key "<<v[0]<<" "<<v[1]<<"\n";
                        int face_ind = l_cell_to_face[node_to_cell[n[k]-1][i]-1][j-1];
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
                        vector<double> face_coor = f_centers[face_ind-1];
                        double nx = face_coor[1] - center_coor[1];
                        double ny = center_coor[0] - face_coor[0];
                        nx /= pow(pow(nx, 2) + pow(ny, 2), 0.5);
                        ny /= pow(pow(nx, 2) + pow(ny, 2), 0.5);
                        double area = pow(pow(face_coor[0] - center_coor[0], 2) + pow(face_coor[1] - center_coor[1], 2), 0.5);
                        double dx = nodes[fn2-1][0] - nodes[n[k]-1][0];
                        double dy = nodes[fn2-1][1] - nodes[n[k]-1][1];
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
        double x1 = nodes[n[0]-1][0];
        double y1 = nodes[n[0]-1][1]; 
        v.erase(v.end()-1);
        v.push_back(2);
        n[1] = l_face_to_node[v];
        double x2 = nodes[n[1]-1][0];
        double y2 = nodes[n[1]-1][1]; 
        if((x1==0.0) and (x2==0.0)){
            P_nodes[n[0]-1] = 0;
            P_nodes[n[1]-1] = 0;
        }else if(((x1==1.0) and (x2==1.0)) or ((y1==1.0) and (y2==1.0)) or ((y1==0.0) and (y2==0.0))){
            // cout<<"nodes "<<n[0]<<" "<<n[1]<<"\n";
            v.clear();
            for(int k = 0; k<2; k++){
                for(int i=0; i<node_to_cell[n[k]-1].size(); i++){
                    v.push_back(node_to_cell[n[k]-1][i]);
                    // cout<<"node# "<<n[k]<<" center: "<<node_to_cell[n[k]][i]<<"\n";
                    double center_coor[2]; 
                    center_coor[0] = centers[node_to_cell[n[k]-1][i]-1][0];
                    center_coor[1] = centers[node_to_cell[n[k]-1][i]-1][1];
                    double nom = 0, denom = 0, valid_denom;
                    for(int j=1; j<4; j++){
                        v.push_back(j);
                        // cout<<"key "<<v[0]<<" "<<v[1]<<"\n";
                        int face_ind = l_cell_to_face[node_to_cell[n[k]-1][i]-1][j-1];
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

                        vector<double> face_coor = f_centers[face_ind-1];
                        double nx = face_coor[1] - center_coor[1];
                        double ny = center_coor[0] - face_coor[0];
                        nx /= pow(pow(nx, 2) + pow(ny, 2), 0.5);
                        ny /= pow(pow(nx, 2) + pow(ny, 2), 0.5);
                        double area = pow(pow(face_coor[0] - center_coor[0], 2) + pow(face_coor[1] - center_coor[1], 2), 0.5);
                        double dx = nodes[fn2-1][0] - nodes[n[k]-1][0];
                        double dy = nodes[fn2-1][1] - nodes[n[k]-1][1];
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
bool is_boundary_face(int face_ind){
    // cout<<"HERE"<<"\n";
    vector<int> key1;
    key1.push_back(face_ind);
    key1.push_back(1);
    vector<int> key2; 
    key2.push_back(face_ind);
    key2.push_back(2);
    // cout<<key1[0]<<" "<<key2[0]<<" "<<key1[1]<<" "<<key2[1]<<"\n";
    // cout<<"ans "<<l_face_to_cell[key1]<<" "<<l_face_to_cell[key2]<<"\n";
    if(l_face_to_cell[key1] == l_face_to_cell[key2]){
        return true;
    }
    return false;
}
void fvm(){
    for(int main=0; main < 100000; main++){
        cout<<main<<"\n";
        // cout<<"8.1"<<"\n";
        for(int i=1; i<center_num+1; i++){
            // vector<int> key;
            // key.push_back(i);
            double au=0, av=0, center_coor[2];
            center_coor[0] = centers[i-1][0];
            center_coor[1] = centers[i-1][1];
            for(int j=1; j<4; j++){
                // key.push_back(j);
                register int face_ind = l_cell_to_face[i-1][j-1];
                vector<int> key1;
                key1.push_back(face_ind);
                double deltx, delty, nx, ny, valid_denom;
                if(is_boundary_face(face_ind)){   
                    vector<double> face_coor = f_centers[face_ind-1];
                    deltx = face_coor[0] - center_coor[0];
                    delty = face_coor[1] - center_coor[1];
                    // key1.push_back(1);
                    nx = sn[face_ind-1][0];
                    // key1.erase(key1.end()-1);
                    // key1.push_back(2);
                    ny = sn[face_ind-1][1];
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
                    au += face_areas[face_ind-1]*(U_faces[face_ind-1] - U_centers_old[i-1])*((d*valid_denom) - U_centers_old[i-1]*nx - V_centers_old[i-1]*ny);
                    av += face_areas[face_ind-1]*(V_faces[face_ind-1] - V_centers_old[i-1])*((d*valid_denom) - U_centers_old[i-1]*nx - V_centers_old[i-1]*ny);
                    // key1.clear();
                }else{
                    key1.push_back(1);
                    nx = sn[face_ind-1][0];
                    register int out_center_ind = l_face_to_cell[key1]; 
                    key1.erase(key1.end()-1);
                    key1.push_back(2);
                    double c2[2];
                    c2[0] = centers[out_center_ind-1][0];
                    c2[1] = centers[out_center_ind-1][1];
                    ny = sn[face_ind-1][1];
                    if(out_center_ind == i){
                        out_center_ind = l_face_to_cell[key1];
                        c2[0] = centers[out_center_ind-1][0];
                        c2[1] = centers[out_center_ind-1][1];
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
                    au += ((U_centers_old[out_center_ind-1] - U_centers_old[i-1]) * face_areas[face_ind-1])*((d*valid_denom) - U_centers_old[i-1]*nx - V_centers_old[i-1]*ny);
                    av += ((V_centers_old[out_center_ind-1] - V_centers_old[i-1]) * face_areas[face_ind-1])*((d*valid_denom) - U_centers_old[i-1]*nx - V_centers_old[i-1]*ny);
                    key1.clear();
                }
                // key.erase(key.end()-1);
            }
            U_centers_new[i-1] = U_centers_old[i-1] + (dt/volume(i)) * au;
            V_centers_new[i-1] = V_centers_old[i-1] + (dt/volume(i)) * av; 
            // key.clear();
        }
        // cout<<"8.2"<<"\n";
        ofstream yFile1("u_centers_new.dat");
        yFile1 << "VARIABLES= \"X\", \"Y\", \"P\", \"U\", \"W\"" <<"\n";
        for(int i = 0; i<center_num; i++){
            yFile1 << centers[i][0] << " " << centers[i][1] <<" "<<U_centers_new[i]<<" "<<V_centers_new[i]<<"\n";
        }
        yFile1.close();
        cell_to_vertex_interpolation_new();
        set_boundary_conditions_U();
        set_boundary_conditions_V();
        // cout<<"8.3"<<"\n";
        for(int i=1; i<center_num+1; i++){
            // vector<int> key;
            // key.push_back(i);
            double au=0, av=0, center_coor[2];
            center_coor[0] = centers[i-1][0];
            center_coor[1] = centers[i-1][1];
            for(int j=1; j<4; j++){
                // key.push_back(j);
                register int face_ind = l_cell_to_face[i-1][j-1];
                vector<int> key1;
                key1.push_back(face_ind);
                double deltx, delty, nx, ny;
                if(is_boundary_face(face_ind)){   
                    vector<double> face_coor = f_centers[face_ind-1];
                    deltx = face_coor[0] - center_coor[0];
                    delty = face_coor[1] - center_coor[1];
                    // key1.push_back(1);
                    nx = sn[face_ind-1][0];
                    // key1.erase(key1.end()-1);
                    // key1.push_back(2);
                    ny = sn[face_ind-1][1];
                    if(((deltx>0) and (nx<0)) or ((deltx<0) and (nx>0))){
                        nx*=(-1);
                    }
                    if(((delty>0) and (ny<0)) or ((delty<0) and (ny>0))){
                        ny*=(-1);
                    }
                    au += (face_areas[face_ind-1])*(U_faces[face_ind-1] - U_centers_new[i-1])*nx;
                    av += (face_areas[face_ind-1])*(V_faces[face_ind-1] - V_centers_new[i-1])*ny;
                    // key1.clear();
                }else{
                    key1.push_back(1);
                    nx = sn[face_ind-1][0];
                    register int out_center_ind = l_face_to_cell[key1]; 
                    key1.erase(key1.end()-1);
                    key1.push_back(2);
                    double c2[2];
                    c2[0] = centers[out_center_ind-1][0];
                    c2[1] = centers[out_center_ind-1][1];
                    ny = sn[face_ind-1][1];
                    if(out_center_ind == i){
                        out_center_ind = l_face_to_cell[key1];
                        c2[0] = centers[out_center_ind-1][0];
                        c2[1] = centers[out_center_ind-1][1];
                    }
                    deltx = c2[0] - center_coor[0];
                    delty = c2[1] - center_coor[1];
                    if(((deltx>0) and (nx<0)) or ((deltx<0) and (nx>0))){
                        nx*=(-1);
                    }
                    if(((delty>0) and (ny<0)) or ((delty<0) and (ny>0))){
                        ny*=(-1);
                    }
                    au += (U_centers_new[out_center_ind-1] - U_centers_new[i-1]) * face_areas[face_ind-1]*nx;
                    av += (V_centers_new[out_center_ind-1] - V_centers_new[i-1]) * face_areas[face_ind-1]*ny;
                    key1.clear();
                }
                // key.erase(key.end()-1);
            }
            f[i-1] = (ro/dt) * (au + av);
            // au = 0;
            // av = 0; 
            // key.clear();
        }
        // cout<<"8.4"<<"\n";
        do{
            cell_to_vertex_interpolation_old();
            set_boundary_conditions_P();
            for(int i=1; i<center_num+1; i++){
                // vector<int> key;
                // key.push_back(i);
                double nom=0, denom=0, center_coor[2];
                center_coor[0] = centers[i-1][0];
                center_coor[1] = centers[i-1][1];
                for(int j=1; j<4; j++){
                    // key.push_back(j);
                    register int face_ind = l_cell_to_face[i-1][j-1];
                    vector<int> key1;
                    key1.push_back(face_ind);
                    double deltx, delty, nx, ny, valid_denom;
                    if(is_boundary_face(face_ind)){   
                        vector<double> face_coor = f_centers[face_ind-1];
                        deltx = face_coor[0] - center_coor[0];
                        delty = face_coor[1] - center_coor[1];
                        // key1.push_back(1);
                        nx = sn[face_ind-1][0];
                        // key1.erase(key1.end()-1);
                        // key1.push_back(2);
                        ny = sn[face_ind-1][1];
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
                        nom += face_areas[face_ind-1]*P_faces[face_ind-1]*valid_denom;
                        denom += face_areas[face_ind-1]*valid_denom;
                        // key1.clear();
                    }else{
                        key1.push_back(1);
                        nx = sn[face_ind-1][0];
                        register int out_center_ind = l_face_to_cell[key1]; 
                        key1.erase(key1.end()-1);
                        key1.push_back(2);
                        double c2[2];
                        c2[0] = centers[out_center_ind-1][0];
                        c2[1] = centers[out_center_ind-1][1];
                        ny = sn[face_ind-1][1];
                        if(out_center_ind == i){
                            out_center_ind = l_face_to_cell[key1];
                            c2[0] = centers[out_center_ind-1][0];
                            c2[1] = centers[out_center_ind-1][1];
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
                        nom += face_areas[face_ind-1]*valid_denom*P_centers_old[out_center_ind-1];
                        denom += face_areas[face_ind-1]*valid_denom;
                        key1.clear();
                    }
                    // key.erase(key.end()-1);
                }
                nom -= f[i-1]*volume(i);
                P_centers_new[i-1] = nom/denom;
                // key.clear();
            }
            max_diff = 0.0; 
            for(int i = 0; i< center_num; i++){
                if(max_diff<fabs(P_centers_new[i] - P_centers_old[i])){
                    max_diff = fabs(P_centers_new[i] - P_centers_old[i]);
                }
                P_centers_old[i] = P_centers_new[i];
            }
            // cout<<max_diff<<"\n";
            // printf("Error: %.6f\n", max_diff);
            // cell_to_vertex_interpolation_old();
            // set_boundary_conditions_P();
        }while(max_diff>eps);
        // cout<<"8.5"<<"\n";
        cell_to_vertex_interpolation_new();
        set_boundary_conditions_P();
        for(int i=1; i<center_num+1; i++){
            // vector<int> key;
            // key.push_back(i);
            double au=0, av=0, center_coor[2];
            center_coor[0] = centers[i-1][0];
            center_coor[1] = centers[i-1][1];
            for(int j=1; j<4; j++){
                // key.push_back(j);
                register int face_ind = l_cell_to_face[i-1][j-1];
                vector<int> key1;
                key1.push_back(face_ind);
                double deltx, delty, nx, ny, valid_denom;
                if(is_boundary_face(face_ind)){   
                    vector<double> face_coor = f_centers[face_ind-1];
                    deltx = face_coor[0] - center_coor[0];
                    delty = face_coor[1] - center_coor[1];
                    // key1.push_back(1);
                    nx = sn[face_ind-1][0];
                    // key1.erase(key1.end()-1);
                    // key1.push_back(2);
                    ny = sn[face_ind-1][1];
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
                    au += face_areas[face_ind-1]*(P_faces[face_ind-1] - P_centers_new[i-1])*nx;
                    av += face_areas[face_ind-1]*(P_faces[face_ind-1] - P_centers_new[i-1])*ny;
                    // key1.clear();
                }else{
                    key1.push_back(1);
                    nx = sn[face_ind-1][0];
                    register int out_center_ind = l_face_to_cell[key1]; 
                    key1.erase(key1.end()-1);
                    key1.push_back(2);
                    double c2[2];
                    c2[0] = centers[out_center_ind-1][0];
                    c2[1] = centers[out_center_ind-1][1];
                    ny = sn[face_ind-1][1];
                    if(out_center_ind == i){
                        out_center_ind = l_face_to_cell[key1];
                        c2[0] = centers[out_center_ind-1][0];
                        c2[1] = centers[out_center_ind-1][1];
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
                    au += face_areas[face_ind-1]*(P_centers_new[out_center_ind-1] - P_centers_new[i-1])*nx;
                    av += face_areas[face_ind-1]*(P_centers_new[out_center_ind-1] - P_centers_new[i-1])*ny;
                    key1.clear();
                }
                // key.erase(key.end()-1);
            }
            U_centers_old[i-1] = U_centers_new[i-1] - (dt/(ro * volume(i))) * au;
            V_centers_old[i-1] = V_centers_new[i-1] - (dt/(ro * volume(i))) * av;
            // key.clear();
        }
        // cout<<"8.6"<<"\n";
        for(int i = 0; i< center_num; i++){
            P_centers_old[i] = P_centers_new[i];
        }
        cell_to_vertex_interpolation_old();
        set_boundary_conditions_P();
        set_boundary_conditions_U();
        set_boundary_conditions_V();
        if(main == 10){
            ofstream yFile10("out10.dat");
            yFile10 << "VARIABLES= \"X\", \"Y\", \"P\", \"U\", \"V\"" <<"\n";
            for(int i = 0; i<node_num; i++){
                yFile10 << nodes[i][0] << " " << nodes[i][1] <<" "<<P_nodes[i]<<" "<<U_nodes[i]<<" "<<V_nodes[i]<<"\n";
            }
            yFile10.close();
        }
        if(main == 100){
            ofstream yFile100("out100.dat");
            yFile100 << "VARIABLES= \"X\", \"Y\", \"P\", \"U\", \"V\"" <<"\n";
            for(int i = 0; i<node_num; i++){
                yFile100 << nodes[i][0] << " " << nodes[i][1] <<" "<<P_nodes[i]<<" "<<U_nodes[i]<<" "<<V_nodes[i]<<"\n";
            }
            yFile100.close();
        }
        if(main == 500){
            ofstream yFile500("out500.dat");
            yFile500 << "VARIABLES= \"X\", \"Y\", \"P\", \"U\", \"V\"" <<"\n";
            for(int i = 0; i<node_num; i++){
                yFile500 << nodes[i][0] << " " << nodes[i][1] <<" "<<P_nodes[i]<<" "<<U_nodes[i]<<" "<<V_nodes[i]<<"\n";
            }
            yFile500.close();
        }
        if(main == 1000){
            ofstream yFile1000("out1000.dat");
            yFile1000 << "VARIABLES= \"X\", \"Y\", \"P\", \"U\", \"V\"" <<"\n";
            for(int i = 0; i<node_num; i++){
                yFile1000 << nodes[i][0] << " " << nodes[i][1] <<" "<<P_nodes[i]<<" "<<U_nodes[i]<<" "<<V_nodes[i]<<"\n";
            }
            yFile1000.close();
        }
        if(main == 50000){
            ofstream yFile50000("out50000.dat");
            yFile50000 << "VARIABLES= \"X\", \"Y\", \"P\", \"U\", \"V\"" <<"\n";
            for(int i = 0; i<node_num; i++){
                yFile50000 << nodes[i][0] << " " << nodes[i][1] <<" "<<P_nodes[i]<<" "<<U_nodes[i]<<" "<<V_nodes[i]<<"\n";
            }
            yFile50000.close();
        }
    }
    // cell_to_vertex_interpolation_old();
    // set_boundary_conditions_P();
    // set_boundary_conditions_U();
    // set_boundary_conditions_V();
    // cout<<"8.7"<<"\n";
    ofstream yFile("out6.dat");
    yFile << "VARIABLES= \"X\", \"Y\", \"P\", \"U\", \"V\"" <<"\n";
    for(int i = 0; i<node_num; i++){
        yFile << nodes[i][0] << " " << nodes[i][1] <<" "<<P_nodes[i]<<" "<<U_nodes[i]<<" "<<V_nodes[i]<<"\n";
    }
    yFile.close();
    // ofstream yFilec("out_centers.dat");
    // yFilec << "VARIABLES= \"X\", \"Y\", \"P\", \"U\", \"V\"" <<"\n";
    // for(int i = 0; i<center_num; i++){
    //     yFilec << centers[i][0] << " " << centers[i][1] <<" "<<P_centers_old[i]<<" "<<U_centers_old[i]<<" "<<V_centers_old[i]<<"\n";
    // }
    // yFilec.close();
}
int main(){
    nodes_coor();
    cout<<"1"<<"\n";
    cell_to_node();
    cout<<"2"<<"\n";
    center_coordinates_generator();
    cout<<"3"<<"\n";
    face_links();
    cout<<"4"<<"\n";
    vector<double> v1;
    vector<double> v2;
    cout<<"Boundary face num "<<boundary_face.size()<<"\n";
    for(int i=0; i<face_num; i++){
        v1 = normal(i+1);
        v2 = face_center(i+1);
        f_centers.push_back(v2);
        double A = area(i+1);
        face_areas.push_back(A);
        sn.push_back(v1);
        v1.clear();
        v2.clear();
        // if(is_boundary_face(i+1)){
        //     cout<<"tRUE"<<is_boundary_face(i+1)<<"\n";
        // }
    }
    // cout<<"BOUN "<<is_boundary_face(636)<<endl;
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