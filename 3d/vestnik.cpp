//Unstructured 3d-mesh
//Solving equation from vestnik
#include <bits/stdc++.h>  
using namespace std;
#define e 2.71828
const int node_num=4750, center_num=23298;
double nodes[node_num][3], centers[center_num][3], U_nodes[node_num],  
    U_centers_old[center_num], U_centers_new[center_num], 
    eps=0.0001, max_diff, f_volume[center_num], lim_x = 0.25, lim_y = 0.5, lim_z = 0.5;
int face_num, l_cell_to_face[center_num][4], l_cell_to_node[center_num][4], source_node, i=0, j=0;
vector<vector<int> > l_face_to_cell, l_face_to_node;
vector<int> node_to_cell[node_num], boundary_face;
vector<vector<double> > sn, f_centers;
vector<double> face_areas, U_faces;
vector<double> split_double(string line, char delimeter){
    string my_str = "";
    vector<double> v;
    for(i=0; i<line.size(); i++){
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
    for(i=0; i<line.size(); i++){
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
    ifstream NodesFile("vestnik_coor.txt");
    vector<double> v; 
    int cnt = 0;
    while (getline(NodesFile, myText)) {
        v = split_double(myText, ' ');
        v.erase(v.begin());
        nodes[cnt][0] = v[0];
        nodes[cnt][1] = v[1];
        nodes[cnt][2] = v[2];
        v.clear();
        cnt++;
    }
    NodesFile.close();
}
void cell_to_node(){
    string myText;
    ifstream NodesFile("vestnik_links.txt");
    vector<int> v; 
    while (getline(NodesFile, myText)) {
        v = split_int(myText, ' ');
        int cell = v[0];
        v.erase(v.begin());
        v.erase(v.begin());
        for(i=0; i<4; i++){
            l_cell_to_node[cell-1][i] = v[i]-1;
            node_to_cell[v[i]-1].push_back(cell-1);
        }
        v.clear();
    }
    NodesFile.close();
}
void center_coordinates_generator(){
    for(i = 0; i < center_num; i++){
        int node1 = l_cell_to_node[i][0], node2 = l_cell_to_node[i][1], 
            node3 = l_cell_to_node[i][2], node4 = l_cell_to_node[i][3];
        centers[i][0] = ((nodes[node1][0] + nodes[node2][0] + nodes[node3][0] + nodes[node4][0]) / 4);
        centers[i][1] = ((nodes[node1][1] + nodes[node2][1] + nodes[node3][1] + nodes[node4][1]) / 4);
        centers[i][2] = ((nodes[node1][2] + nodes[node2][2] + nodes[node3][2] + nodes[node4][2]) / 4);
    }
}
void face_links(){
    map< set<int>, vector<int> > nodes_to_cell;
    vector<int> key1, trg_nodes;
    int cnt = 1;
    set<int> nodes_ind;
    vector< vector<int> > c_t_f(center_num);
    int arr[4][3] = {{0,1,2}, {0,1,3}, {0,2,3}, {1,2,3}}, m;
    for(i=0; i<center_num; i++){
        for(j=0; j<4; j++){
            register int trg_node = l_cell_to_node[i][j];
            trg_nodes.push_back(trg_node);
        }
        for(j=0; j<4; j++){
            nodes_ind.insert(trg_nodes[arr[j][0]]);
            nodes_ind.insert(trg_nodes[arr[j][1]]); 
            nodes_ind.insert(trg_nodes[arr[j][2]]);
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
    register int f_cnt = 0;
    while(m_it!=nodes_to_cell.end()){
        key1 = m_it->second; 
        set <int>::iterator s_it = m_it->first.begin();
        vector<int> temp;
        while(s_it != m_it->first.end()){
            temp.push_back(*s_it);
            s_it++;
        }
        l_face_to_node.push_back(temp);
        temp.clear();
        for(i=0; i<key1.size(); i++){
            c_t_f[key1[i]].push_back(f_cnt);
            temp.push_back(key1[i]); 
            if(key1.size() == 1){
                boundary_face.push_back(f_cnt);
                temp.push_back(key1[i]); 
                continue;
            }   
        }
        l_face_to_cell.push_back(temp);
        temp.clear();
        f_cnt++;
        m_it++;
        key1.clear();
    }
    for(i=0; i<center_num; i++){
        for(j=0; j<c_t_f[i].size(); j++){
            l_cell_to_face[i][j] = c_t_f[i][j];
        }
    }
    cout<<l_face_to_node.size()<<endl;
}
vector<double> face_center(int face_ind){
    double n1[3], n2[3], n3[3];
    for(i=0; i<3; i++){
        n1[i] = nodes[l_face_to_node[face_ind][0]][i];
        n2[i] = nodes[l_face_to_node[face_ind][1]][i];
        n3[i] = nodes[l_face_to_node[face_ind][2]][i]; 
    }    
    vector<double> ans; 
    ans.push_back((n1[0] + n2[0] + n3[0])/3);
    ans.push_back((n1[1] + n2[1] + n3[1])/3);
    ans.push_back((n1[2] + n2[2] + n3[2])/3);
    return ans;
}
vector<double> normal(int face_ind){
    vector <double> ans; 
    double n1[3], n2[3], n3[3];
    for(i=0; i<3; i++){
        n1[i] = nodes[l_face_to_node[face_ind][0]][i];
        n2[i] = nodes[l_face_to_node[face_ind][1]][i];
        n3[i] = nodes[l_face_to_node[face_ind][2]][i]; 
    }  
    ans.push_back((n2[1] - n1[1])*(n3[2]-n1[2]) - (n2[2] - n1[2])*(n3[1] - n1[1])); 
    ans.push_back((n2[2] - n1[2])*(n3[0]-n1[0]) - (n2[0] - n1[0])*(n3[2] - n1[2])); 
    ans.push_back((n2[0] - n1[0])*(n3[1]-n1[1]) - (n2[1] - n1[1])*(n3[0] - n1[0])); 
    double len = pow(pow(ans[0], 2) + pow(ans[1], 2) + pow(ans[2], 2), 0.5); 
    for(i=0; i<ans.size(); i++){
        ans[i] /= len;
    }
    return ans; 
}
double area(int face_ind){
    vector <double> ans; 
    double n1[3], n2[3], n3[3];
    for(i=0; i<3; i++){
        n1[i] = nodes[l_face_to_node[face_ind][0]][i];
        n2[i] = nodes[l_face_to_node[face_ind][1]][i];
        n3[i] = nodes[l_face_to_node[face_ind][2]][i]; 
    }    
    ans.push_back((n2[1] - n1[1])*(n3[2]-n1[2]) - (n2[2] - n1[2])*(n3[1] - n1[1])); 
    ans.push_back((n2[2] - n1[2])*(n3[0]-n1[0]) - (n2[0] - n1[0])*(n3[2] - n1[2])); 
    ans.push_back((n2[0] - n1[0])*(n3[1]-n1[1]) - (n2[1] - n1[1])*(n3[0] - n1[0])); 
    return 0.5 * pow((pow(ans[0], 2) + pow(ans[1], 2) + pow(ans[2], 2)) , 0.5);
}
double volume(int cell_ind){
    double x[3], y[3], z[3];
    for(i=0; i<3; i++){
        x[i] = nodes[l_cell_to_node[cell_ind][i+1]][0] - nodes[l_cell_to_node[cell_ind][0]][0];
        y[i] = nodes[l_cell_to_node[cell_ind][i+1]][1] - nodes[l_cell_to_node[cell_ind][0]][1];
        z[i] = nodes[l_cell_to_node[cell_ind][i+1]][2] - nodes[l_cell_to_node[cell_ind][0]][2];
    }
    return (x[0]*y[1]*z[2] + x[1]*y[2]*z[0] + y[0]*z[1]*x[2] - x[2]*y[1]*z[0] - x[1]*y[0]*z[2] - x[0]*y[2]*z[1]) / 6.0;
}
void cell_to_vertex_interpolation_new(){ 
    for(int node_ind=0; node_ind<node_num; node_ind++){
        double node_coor[3], center_coor[3], nom_U = 0, denom = 0, U, d; 
        for(i=0; i<3; i++){
            node_coor[i]= nodes[node_ind][i]; 
        }
        int cell_ind; 
        for(i=0; i<node_to_cell[node_ind].size(); i++){ 
            cell_ind = node_to_cell[node_ind][i]; 
            for(j=0; j<3; j++){
                center_coor[j] = centers[cell_ind][j]; 
            }
            d = pow(pow((center_coor[0] - node_coor[0]), 2) + pow((center_coor[1] - node_coor[1]), 2)  + pow((center_coor[2] - node_coor[2]), 2), 0.5); 
            U = U_centers_new[cell_ind]; 
            nom_U += U/d; 
            denom += 1/d; 
        } 
        U_nodes[node_ind] = nom_U/denom;
    }
}
void set_boundary_conditions_U(){
    double x[3], y[3], z[3];
    int b_f, n[3];
    vector <double> face_coor;
    for(i=0; i<boundary_face.size(); i++){
        b_f = boundary_face[i];
        for(j = 0; j<3; j++){
            n[j] = l_face_to_node[b_f][j];
            x[j] = nodes[n[j]][0], y[j] = nodes[n[j]][1], z[j] = nodes[n[j]][2];
            if(!((z[j]==0.5) and (z[j]==0.5) and (z[j]==0.5))){
                U_nodes[n[j]] = cos(3*x[j] + y[j] - 2*z[j]) + exp(x[j]-z[j]) + 1;
            }   
        }
        U_faces[b_f] = (U_nodes[n[0]] + U_nodes[n[1]] + U_nodes[n[2]])/3;
    }
}
bool is_boundary_face(int face_ind){
    if(l_face_to_cell[face_ind][0] == l_face_to_cell[face_ind][1]){
        return true;
    }
    return false;
}
void fvm(){
    double center_coor[3], delt[3], n[3], valid_denom, nom, denom, c2[3];
    register int out_center_ind, face_ind;
    vector<double> face_coor;
    do{
        for(i=0; i<center_num; i++){
            nom = 0, denom = 0;
            for(int l=0; l<3; l++){
                center_coor[l] = centers[i][l];
            }
            for(j=0; j<4; j++){
                face_ind = l_cell_to_face[i][j];
                if(is_boundary_face(face_ind)){   
                    face_coor = f_centers[face_ind];
                    for(int l=0; l<3; l++){
                        delt[l] = face_coor[l] - center_coor[l];
                        n[0] = sn[face_ind][0];
                        if(((delt[l]>0) and (n[l]<0)) or ((delt[l]<0) and (n[l]>0))){
                            n[l]*=(-1);
                        }
                    }
                    valid_denom=0.0;
                    if(delt[0] != 0){
                        valid_denom += (n[0]/delt[0]);
                    }
                    if(delt[1] != 0){
                        valid_denom += (n[1]/delt[1]);
                    }
                    if(face_coor[2] == 0.5){
                        nom += (U_faces[face_ind]*valid_denom + (2*sin( 3 * face_coor[0] + face_coor[1] - 1) - exp(face_coor[0] - 0.5)))*face_areas[face_ind];
                    }else{
                        if(delt[2] != 0){
                            valid_denom += (n[2]/delt[2]);
                        }
                        nom += face_areas[face_ind]*U_faces[face_ind]*valid_denom;
                    }  
                }else{
                    out_center_ind = l_face_to_cell[face_ind][0]; 
                    if(out_center_ind == i){
                        out_center_ind = l_face_to_cell[face_ind][1];
                    }
                    valid_denom=0.0;
                    for(int l =0; l<3; l++){
                        n[l] = sn[face_ind][l];
                        c2[l] = centers[out_center_ind][l];
                        delt[l] = c2[l] - center_coor[l];
                        if(((delt[l]>0) and (n[l]<0)) or ((delt[l]<0) and (n[l]>0))){
                            n[l]*=(-1);
                        }
                        if(delt[l] != 0){
                            valid_denom += (n[l]/delt[l]);
                        }
                    }
                    nom += U_centers_old[out_center_ind]* face_areas[face_ind]*valid_denom;
                }
                denom += face_areas[face_ind]*valid_denom;
            }
            nom += f_volume[i] * (4*cos(3*center_coor[0] + center_coor[1] - 2*center_coor[2]) - pow(12, (e - center_coor[2])) - 10);
            denom -= 10*f_volume[i];
            U_centers_new[i] = nom/denom;
        }
        cell_to_vertex_interpolation_new();
        set_boundary_conditions_U();
        max_diff = 0.0;
        for(i=0; i<center_num; i++){
            max_diff = max(fabs(U_centers_old[i] - U_centers_new[i]), max_diff);
            U_centers_old[i] = U_centers_new[i];
        }
        printf("max diff: %.12f\n", max_diff);
    }while(max_diff>eps);
}
double analytical(int node_ind){
    double n_coor[3];
    for(i=0; i<3; i++){
        n_coor[i] = nodes[node_ind][i];
    }
    return cos(3*n_coor[0] + n_coor[1] - 2*n_coor[2]) + exp(n_coor[0]-n_coor[2]) + 1;
}
int main(){
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    int i=0;
    nodes_coor();
    cout<<"1"<<"\n";
    cell_to_node();
    cout<<"2"<<"\n";
    center_coordinates_generator();
    cout<<"3"<<"\n";
    face_links();
    cout<<"4"<<"\n";
    cout<<"Bound "<<boundary_face.size()<<endl;
    vector<double> v1, v2;
    double A, V;
    cout<<face_num;
    for(i=0; i<face_num; i++){
        v1 = normal(i);
        v2 = face_center(i);
        f_centers.push_back(v2);
        A = area(i);
        face_areas.push_back(A);
        sn.push_back(v1);
        v1.clear();
        v2.clear();
        U_faces.push_back(0);
    }
    for(i=0; i<center_num; i++){
        V = volume(i);
        f_volume[i] = V;
        U_centers_new[i]=0;
        U_centers_old[i]=0;
    }
    for(i=0; i<node_num; i++){
        U_nodes[i] = 0;
    }
    set_boundary_conditions_U();
    cout<<"5"<<"\n";
    fvm();
    cout<<"6"<<"\n";
    ofstream MFile("3d_sol.dat");
    MFile << "VARIABLES= \"X\", \"Y\", \"Z\", \"U\"" <<endl;
    max_diff = 0.0;
    for(int i=0; i<node_num; i++){
        MFile << nodes[i][0] << " " << nodes[i][1] << " " << nodes[i][2]<<" "<<U_centers_new[i]<<endl;
        max_diff = max(max_diff, fabs(U_centers_new[i] - analytical(i)));
    }
    printf("max error: %.12f\n", max_diff);
    MFile.close();
    return 0;
}