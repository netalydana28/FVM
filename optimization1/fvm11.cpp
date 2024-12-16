//Navier Stokes on unstructured mesh optimization 
#include <bits/stdc++.h>  
using namespace std;
const int node_num=4042, center_num=7882;
double nodes[node_num][2], centers[center_num][2], 
    P_nodes[node_num], U_nodes[node_num], 
    V_nodes[node_num], C_nodes[node_num],
    P_centers_old[center_num], P_centers_new[center_num], 
    U_centers_old[center_num], U_centers_new[center_num], 
    V_centers_old[center_num], V_centers_new[center_num], 
    C_centers_old[center_num], C_centers_new[center_num], f[center_num], 
    eps=0.0001, nyu = 1.0/400.0, dt = 1.0e-006, ro = 1.0, d = 1.0/500.0, 
    max_diff, f_volume[center_num];
int face_num, l_cell_to_face[center_num][3], l_cell_to_node[center_num][3], i=0, j=0;
vector<vector<int> > l_face_to_cell, l_face_to_node;
vector<int> node_to_cell[node_num], boundary_face, source_node;
vector<vector<double> > sn, f_centers, speed_over_time;
vector<double> face_areas, U_faces, V_faces;
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
            if(my_str[0]== 't' or my_str[0]== 'r' or my_str[0]== 'i'){
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
void get_source(){
    string myText;
    ifstream SourceFile("concentration.txt");
    while (getline(SourceFile, myText)) {
        source_node = split_int(myText, ' ');
    }
    cout<<source_node[0]<<endl;
    SourceFile.close();
}
void get_speed(){
    string myText;
    ifstream SourceFile("uspeed.txt");
    while(getline(SourceFile, myText)) {
        speed_over_time.push_back(split_double(myText, ' '));
    }
    cout<<speed_over_time[1][0]<<" "<<speed_over_time[1][1]<<" "<<speed_over_time[1][2]<<endl;
    SourceFile.close();
}
void nodes_coor(){
    string myText;
    ifstream NodesFile("mesh_coor.txt");
    vector<double> v; 
    int cnt = 0;
    while (getline(NodesFile, myText)) {
        v = split_double(myText, ' ');
        v.erase(v.begin());
        v.erase(v.end()-1);
        nodes[cnt][0] = v[0];
        nodes[cnt][1] = v[1];
        v.clear();
        cnt++;
    }
    NodesFile.close();
    // cout<<nodes[807][0]<<" "<<nodes[807][1]<<endl;
}
void cell_to_node(){

    string myText;
    ifstream NodesFile("mesh_links.txt");
    vector<int> v; 
    while (getline(NodesFile, myText)) {
        v = split_int(myText, ' ');
        int cell = v[0];
        v.erase(v.begin());
        v.erase(v.begin());
        // cout<<"here1"<<endl;
        for(i=0; i<3; i++){
            // cout<<"here2"<<endl;
            l_cell_to_node[cell-1][i] = v[i]-1;
            node_to_cell[v[i]-1].push_back(cell-1);
        }
        v.clear();
    }
    NodesFile.close();
    cout<<"here"<<endl;
}
void center_coordinates_generator(){
    for(i = 0; i < center_num; i++){
        int node1 = l_cell_to_node[i][0], node2 = l_cell_to_node[i][1], node3 = l_cell_to_node[i][2];
        centers[i][0] = ((nodes[node1][0] + nodes[node2][0] + nodes[node3][0]) / 3);
        centers[i][1] = ((nodes[node1][1] + nodes[node2][1] + nodes[node3][1]) / 3);
    }
}
void face_links(){
    map< set<int>, vector<int> > nodes_to_cell;
    vector<int> key1, trg_nodes;
    int cnt = 1;
    set<int> nodes_ind;
    vector< vector<int> > c_t_f(center_num);
    for(i=0; i<center_num; i++){
        for(j=0; j<3; j++){
            register int trg_node = l_cell_to_node[i][j];
            trg_nodes.push_back(trg_node);
        }
        register int m;
        for(j=0; j<3; j++){
            
            if(j==2){
                m = 0;
            }else{
                m = j+1;
            }
            nodes_ind.insert(trg_nodes[j]); 
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
    for(i=0; i<center_num; i++){
        f[i] = 0;
    }
    // cout<<l_face_to_cell.size()<<endl;
    cout<<l_face_to_node.size()<<endl;
}
vector<double> face_center(int face_ind){
    double n1[2], n2[2];
    n1[0] = nodes[l_face_to_node[face_ind][0]][0];
    n1[1] = nodes[l_face_to_node[face_ind][0]][1];
    n2[0] = nodes[l_face_to_node[face_ind][1]][0];
    n2[1] = nodes[l_face_to_node[face_ind][1]][1];    
    vector<double> ans; 
    ans.push_back((n1[0] + n2[0])/2);
    ans.push_back((n1[1] + n2[1])/2);
    return ans;
}
vector<double> normal(int face_ind){
    vector <double> ans; 
    double n1[2], n2[2];
    n1[0] = nodes[l_face_to_node[face_ind][0]][0];
    n1[1] = nodes[l_face_to_node[face_ind][0]][1];
    n2[0] = nodes[l_face_to_node[face_ind][1]][0];
    n2[1] = nodes[l_face_to_node[face_ind][1]][1];
    ans.push_back(n1[1] - n2[1]); 
    ans.push_back(n2[0] - n1[0]);  
    double len = pow(pow(ans[0], 2) + pow(ans[1], 2), 0.5); 
    for(i=0; i<ans.size(); i++){
        ans[i] /= len;
    }
    return ans; 
}
double area(int face_ind){
    double n1[2], n2[2];
    n1[0] = nodes[l_face_to_node[face_ind][0]][0];
    n1[1] = nodes[l_face_to_node[face_ind][0]][1];
    n2[0] = nodes[l_face_to_node[face_ind][1]][0];
    n2[1] = nodes[l_face_to_node[face_ind][1]][1];
    return pow((pow((n1[0] - n2[0]), 2) + pow((n1[1] - n2[1]), 2)) , 0.5);
}
double volume(int cen_ind){
    int n1 = l_cell_to_node[cen_ind][0], n2 = l_cell_to_node[cen_ind][1], n3 = l_cell_to_node[cen_ind][2];
    double x1 = nodes[n1][0], x2 = nodes[n2][0], x3 = nodes[n3][0], 
        y1 = nodes[n1][1], y2 = nodes[n2][1], y3 = nodes[n3][1];
    return (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) / 2;
}
void cell_to_vertex_interpolation(){ 
    for(int node_ind=0; node_ind<node_num; node_ind++){
        double node_coor[2], center_coor[2], nom_U = 0, denom = 0, nom_V = 0, nom_P = 0, nom_C = 0, P, V, U, C, d;  
        node_coor[0]= nodes[node_ind][0]; 
        node_coor[1]= nodes[node_ind][1];  
        int cell_ind; 
        for(i=0; i<node_to_cell[node_ind].size(); i++){ 
            cell_ind = node_to_cell[node_ind][i]; 
            center_coor[0] = centers[cell_ind][0]; 
            center_coor[1] = centers[cell_ind][1]; 
            d = pow(pow((center_coor[0] - node_coor[0]), 2) + pow((center_coor[1] - node_coor[1]), 2), 0.5); 
            P = P_centers_old[cell_ind],  U = U_centers_old[cell_ind], 
                V = V_centers_old[cell_ind], C = C_centers_old[cell_ind]; 
            nom_P += P/d;
            nom_U += U/d; 
            nom_V += V/d; 
            nom_C += C/d; 
            denom += 1/d; 
        } 
        U_nodes[node_ind] = nom_U/denom;
        V_nodes[node_ind] = nom_V/denom;
        P_nodes[node_ind] = nom_P/denom;
        C_nodes[node_ind] = nom_C/denom;
    }
}
void write_file(){
    ofstream yFile("out5.dat");
    // yFile<< "VARIABLES= \"X\", \"Y\", \"P\", \"U\", \"V\", \"C\"" <<"\n";
    for(i = 0; i<node_num; i++){
        yFile << nodes[i][0] << " " << nodes[i][1] <<" "<<P_nodes[i]<<" "<<U_nodes[i]<<" "<<V_nodes[i]<<" "<<C_nodes[i]<<"\n";
    }
    yFile.close();
} 
void write_file_1st(){
    ofstream yFile("out1st_iteration.dat");
    // yFile<< "VARIABLES= \"X\", \"Y\", \"P\", \"U\", \"V\", \"C\"" <<"\n";
    for(i = 0; i<node_num; i++){
        yFile << nodes[i][0] << " " << nodes[i][1] <<" "<<P_nodes[i]<<" "<<U_nodes[i]<<" "<<V_nodes[i]<<" "<<C_nodes[i]<<"\n";
    }
    yFile.close();
}  
void set_boundary_conditions(double side[], double u0, double v0){
    double x_side[2], y_side[2];
    int b_f, n[2];
    for(j =0; j<source_node.size(); j++){
        for(i=0; i<node_to_cell[source_node[j]].size(); i++){
            C_nodes[l_cell_to_node[node_to_cell[source_node[j]][i]][0]] = 1;
            C_nodes[l_cell_to_node[node_to_cell[source_node[j]][i]][1]] = 1;
            C_nodes[l_cell_to_node[node_to_cell[source_node[j]][i]][2]] = 1;
        }
    }
    for(int l=0; l<boundary_face.size(); l++){
        b_f = boundary_face[l];
        n[0] = l_face_to_node[b_f][0];
        n[1] = l_face_to_node[b_f][1];
        x_side[0] = nodes[n[0]][0]; 
        x_side[1] = nodes[n[1]][0];
        y_side[0] = nodes[n[0]][1]; 
        y_side[1] = nodes[n[1]][1];
        if(x_side[0] == x_side[1] and x_side[0] == side[0]){
            U_nodes[n[0]] = u0;
            U_nodes[n[1]] = u0; 
            V_nodes[n[0]] = v0;
            V_nodes[n[1]] = v0; 
            // cout<<"Speed"<<endl; 
        }else if(y_side[0] == y_side[1] and y_side[0] == side[1]){
            U_nodes[n[0]] = u0;
            U_nodes[n[1]] = u0; 
            V_nodes[n[0]] = v0;
            V_nodes[n[1]] = v0;  
            // cout<<"Speed"<<endl; 
        }else{
            U_nodes[n[0]] = 0;
            U_nodes[n[1]] = 0; 
            V_nodes[n[0]] = 0;
            V_nodes[n[1]] = 0;
        }
        V_faces[b_f] = (V_nodes[n[0]] + V_nodes[n[1]])/2;
        U_faces[b_f] = (U_nodes[n[0]] + U_nodes[n[1]])/2;
    }
    for(int l=0; l<boundary_face.size(); l++){
        b_f = boundary_face[l];
        n[0] = l_face_to_node[b_f][0];
        n[1] = l_face_to_node[b_f][1];
        x_side[0] = nodes[n[0]][0]; 
        x_side[1] = nodes[n[1]][0];
        y_side[0] = nodes[n[0]][1]; 
        y_side[1] = nodes[n[1]][1];
        if(x_side[0] == x_side[1]){
            if(x_side[0] == side[0]){
                P_nodes[n[0]] = 1;
                P_nodes[n[1]] = 1;
                // cout<<"Pressure"<<endl;  
            }
        }else if(y_side[0] == y_side[1]){
            if(y_side[0] == side[1]){
                P_nodes[n[0]] = 1;
                P_nodes[n[1]] = 1;  
            }
        }
        // P_faces[b_f] = (P_nodes[n[0]] + P_nodes[n[1]])/2;
    }
}
bool is_boundary_face(int face_ind){
    if(l_face_to_cell[face_ind][0] == l_face_to_cell[face_ind][1]){
        return true;
    }
    return false;
}
void fvm(double degrees, double t_start, double t_end){
    double u0 = cos(degrees*(M_PI/180)); 
    double v0 = sin(degrees*(M_PI/180)); 
    u0 /= pow(pow(u0, 2)+pow(v0, 2), 0.5);
    v0 /= pow(pow(u0, 2)+pow(v0, 2), 0.5);
    cout<<"U0 "<<u0<<" V0 "<<v0<<endl;
    double side[2];
    if(u0>0){
        side[0] = 0;
    }else{
        side[0] = 1;
    }
    if(v0>0){
        side[1] = 0;
    }else{
        side[1] = 1;
    }
    cout<<"Side "<<side[0]<<" "<<side[1]<<endl;
    set_boundary_conditions(side, u0, v0);
    write_file_1st();
    cout<<"5"<<"\n";
    double au, av, ac, center_coor[2], deltx, delty, nx, ny, valid_denom, nom, denom, c2[2];
    register int out_center_ind, face_ind;
    vector<double> face_coor;
    for(int main=0; main < (t_end - t_start)/dt +1; main++){
        printf("%d \n", main);
        // write_file_1st();
        // double start_time =  clock();
        for(i=0; i<center_num; i++){
            au=0, av=0, ac = 0;
            center_coor[0] = centers[i][0];
            center_coor[1] = centers[i][1];
            for(j=0; j<3; j++){
                face_ind = l_cell_to_face[i][j];
                if(is_boundary_face(face_ind)){   
                    face_coor = f_centers[face_ind];
                    deltx = face_coor[0] - center_coor[0];
                    valid_denom=0;
                    if(deltx != 0){
                        nx = sn[face_ind][0];
                        if(((deltx>0) and (nx<0)) or ((deltx<0) and (nx>0))){
                            nx*=(-1);
                        }
                        valid_denom += (nx/deltx);
                    }
                    delty = face_coor[1] - center_coor[1];
                    if(delty != 0){
                        ny = sn[face_ind][1];
                        if(((delty>0) and (ny<0)) or ((delty<0) and (ny>0))){
                            ny*=(-1);
                        }
                        valid_denom += (ny/delty);
                    }
                    au += face_areas[face_ind]*((U_faces[face_ind] - U_centers_old[i])*(nyu*valid_denom) - U_faces[face_ind]*U_centers_old[i]*nx - U_faces[face_ind]*V_centers_old[i]*ny);
                    av += face_areas[face_ind]*((V_faces[face_ind] - V_centers_old[i])*(nyu*valid_denom) - V_faces[face_ind]*U_centers_old[i]*nx - V_faces[face_ind]*V_centers_old[i]*ny);
                }else{
                    out_center_ind = l_face_to_cell[face_ind][0]; 
                    if(out_center_ind == i){
                        out_center_ind = l_face_to_cell[face_ind][1]; 
                    }
                    c2[0] = centers[out_center_ind][0];
                    c2[1] = centers[out_center_ind][1];
                    deltx = c2[0] - center_coor[0];
                    valid_denom=0;
                    if(deltx != 0){
                        nx = sn[face_ind][0];
                        if(((deltx>0) and (nx<0)) or ((deltx<0) and (nx>0))){
                            nx*=(-1);
                        }
                        valid_denom += (nx/deltx);
                    }
                    delty = c2[1] - center_coor[1];
                    if(delty != 0){
                        ny = sn[face_ind][1];
                        if(((delty>0) and (ny<0)) or ((delty<0) and (ny>0))){
                            ny*=(-1);
                        }
                        valid_denom += (ny/delty);
                    }
                    au += face_areas[face_ind]*((U_centers_old[out_center_ind] - U_centers_old[i]) *(nyu*valid_denom) - U_centers_old[out_center_ind]*U_centers_old[i]*nx - U_centers_old[out_center_ind]*V_centers_old[i]*ny);
                    av += face_areas[face_ind]*((V_centers_old[out_center_ind] - V_centers_old[i]) * (nyu*valid_denom) - V_centers_old[out_center_ind]*U_centers_old[i]*nx - V_centers_old[out_center_ind]*V_centers_old[i]*ny);
                    ac += face_areas[face_ind]*((C_centers_old[out_center_ind] - C_centers_old[i]) * (d*valid_denom) - C_centers_old[out_center_ind]*U_centers_old[i]*nx - C_centers_old[out_center_ind]*V_centers_old[i]*ny);
                }
            }
            U_centers_new[i] = U_centers_old[i] + (dt/f_volume[i]) * au;
            V_centers_new[i] = V_centers_old[i] + (dt/f_volume[i]) * av; 
            C_centers_new[i] = C_centers_old[i] + (dt/f_volume[i]) * ac;
        }
        for(i=0; i<center_num; i++){
            au=0, av=0;
            center_coor[0] = centers[i][0];
            center_coor[1] = centers[i][1];
            for(j=0; j<3; j++){
                face_ind = l_cell_to_face[i][j];
                if(is_boundary_face(face_ind)){   
                    face_coor = f_centers[face_ind];
                    deltx = face_coor[0] - center_coor[0];
                    delty = face_coor[1] - center_coor[1];
                    nx = sn[face_ind][0];
                    ny = sn[face_ind][1];
                    if(((deltx>0) and (nx<0)) or ((deltx<0) and (nx>0))){
                        nx*=(-1);
                    }
                    if(((delty>0) and (ny<0)) or ((delty<0) and (ny>0))){
                        ny*=(-1);
                    }
                    au += face_areas[face_ind]*U_faces[face_ind]*nx;
                    av += face_areas[face_ind]*V_faces[face_ind]*ny;
                }else{
                    nx = sn[face_ind][0];
                    out_center_ind = l_face_to_cell[face_ind][0]; 
                    ny = sn[face_ind][1];
                    if(out_center_ind == i){
                        out_center_ind = l_face_to_cell[face_ind][1];
                    }
                    c2[0] = centers[out_center_ind][0];
                    c2[1] = centers[out_center_ind][1];
                    deltx = c2[0] - center_coor[0];
                    delty = c2[1] - center_coor[1];
                    if(((deltx>0) and (nx<0)) or ((deltx<0) and (nx>0))){
                        nx*=(-1);
                    }
                    if(((delty>0) and (ny<0)) or ((delty<0) and (ny>0))){
                        ny*=(-1);
                    }
                    au += U_centers_new[out_center_ind]*face_areas[face_ind]*nx;
                    av += V_centers_new[out_center_ind]*face_areas[face_ind]*ny;
                }
            }
            f[i] = (ro/dt) * (au + av);
        }
        do{
            for(i=0; i<center_num; i++){
                nom=0, denom=0;
                center_coor[0] = centers[i][0];
                center_coor[1] = centers[i][1];
                for(j=0; j<3; j++){
                    face_ind = l_cell_to_face[i][j];
                    if(is_boundary_face(face_ind)){   
                        face_coor = f_centers[face_ind];
                        deltx = face_coor[0] - center_coor[0];
                        valid_denom=0;
                        if(deltx != 0){
                            nx = sn[face_ind][0];
                            if(((deltx>0) and (nx<0)) or ((deltx<0) and (nx>0))){
                                nx*=(-1);
                            }
                            valid_denom += (nx/deltx);
                        }
                        delty = face_coor[1] - center_coor[1];
                        if(delty != 0){
                            ny = sn[face_ind][1];
                            if(((delty>0) and (ny<0)) or ((delty<0) and (ny>0))){
                                ny*=(-1);
                            }
                            valid_denom += (ny/delty);
                        }
                        if(face_coor[0] == side[0] or face_coor[1] == side[1]){
                            denom += face_areas[face_ind]*valid_denom;
                        }
                    }else{
                        out_center_ind = l_face_to_cell[face_ind][0]; 
                        if(out_center_ind == i){
                            out_center_ind = l_face_to_cell[face_ind][1]; 
                        }
                        c2[0] = centers[out_center_ind][0];
                        c2[1] = centers[out_center_ind][1];
                        deltx = c2[0] - center_coor[0];
                        valid_denom=0;
                        if(deltx != 0){
                            nx = sn[face_ind][0];
                            if(((deltx>0) and (nx<0)) or ((deltx<0) and (nx>0))){
                                nx*=(-1);
                            }
                            valid_denom += (nx/deltx);
                        }
                        delty = c2[1] - center_coor[1];
                        if(delty != 0){
                            ny = sn[face_ind][1];
                            if(((delty>0) and (ny<0)) or ((delty<0) and (ny>0))){
                                ny*=(-1);
                            }
                            valid_denom += (ny/delty);
                        }
                        nom += face_areas[face_ind]*valid_denom*P_centers_old[out_center_ind];
                        denom += face_areas[face_ind]*valid_denom;
                    }
                }
                nom -= f[i]*f_volume[i];
                P_centers_new[i] = nom/denom;
            }
            max_diff = 0.0; 
            for(i = 0; i< center_num; i++){
                if(max_diff<fabs(P_centers_new[i] - P_centers_old[i])){
                    max_diff = fabs(P_centers_new[i] - P_centers_old[i]);
                }
                P_centers_old[i] = P_centers_new[i];
            }
            printf("   %.6f\n", max_diff);
        }while(max_diff>eps);
        for(i=0; i<center_num; i++){
            au=0, av=0;
            center_coor[0] = centers[i][0];
            center_coor[1] = centers[i][1];
            for(j=0; j<3; j++){
                face_ind = l_cell_to_face[i][j];
                if(is_boundary_face(face_ind)){   
                    face_coor = f_centers[face_ind];
                    deltx = face_coor[0] - center_coor[0];
                    delty = face_coor[1] - center_coor[1];
                    nx = sn[face_ind][0];
                    ny = sn[face_ind][1];
                    if(((deltx>0) and (nx<0)) or ((deltx<0) and (nx>0))){
                        nx*=(-1);
                    }
                    if(((delty>0) and (ny<0)) or ((delty<0) and (ny>0))){
                        ny*=(-1);
                    }
                    if(face_coor[0] == side[0] or face_coor[1] == side[1]){
                        au += face_areas[face_ind]*nx;
                        av += face_areas[face_ind]*ny;
                    }
                }else{
                    nx = sn[face_ind][0];
                    out_center_ind = l_face_to_cell[face_ind][0]; 
                    ny = sn[face_ind][1];
                    if(out_center_ind == i){
                        out_center_ind = l_face_to_cell[face_ind][1];
                    }
                    c2[0] = centers[out_center_ind][0];
                    c2[1] = centers[out_center_ind][1];
                    deltx = c2[0] - center_coor[0];
                    delty = c2[1] - center_coor[1];
                    if(((deltx>0) and (nx<0)) or ((deltx<0) and (nx>0))){
                        nx*=(-1);
                    }
                    if(((delty>0) and (ny<0)) or ((delty<0) and (ny>0))){
                        ny*=(-1);
                    }
                    au += face_areas[face_ind]*P_centers_new[out_center_ind]*nx;
                    av += face_areas[face_ind]*P_centers_new[out_center_ind]*ny;
                }
            }
            U_centers_old[i] = U_centers_new[i] - (dt/(ro * f_volume[i])) * au;
            V_centers_old[i] = V_centers_new[i] - (dt/(ro * f_volume[i]))* av;
            C_centers_old[i] = C_centers_new[i];
            P_centers_old[i] = P_centers_new[i];
        }
        for(j =0; j<source_node.size(); j++){
            for(i=0; i<node_to_cell[source_node[j]].size(); i++){
                C_centers_old[node_to_cell[source_node[j]][i]] = 1;
            }
        }
        // if(main== 0){
        //     cell_to_vertex_interpolation();
        //     for(j=0; j<source_node.size(); j++){
                
        //         for(i=0; i<node_to_cell[source_node[j]].size(); i++){
        //             C_nodes[l_cell_to_node[node_to_cell[source_node[j]][i]][0]] = 1;
        //             C_nodes[l_cell_to_node[node_to_cell[source_node[j]][i]][1]] = 1;
        //             C_nodes[l_cell_to_node[node_to_cell[source_node[j]][i]][2]] = 1;
        //         }
        //     }
        //     set_boundary_conditions(side, u0, v0);
        //     write_file_1st();
        // }
        if((main)%5000 == 0){
            cell_to_vertex_interpolation();
            for(j=0; j<source_node.size(); j++){
                
                for(i=0; i<node_to_cell[source_node[j]].size(); i++){
                    C_nodes[l_cell_to_node[node_to_cell[source_node[j]][i]][0]] = 1;
                    C_nodes[l_cell_to_node[node_to_cell[source_node[j]][i]][1]] = 1;
                    C_nodes[l_cell_to_node[node_to_cell[source_node[j]][i]][2]] = 1;
                }
            }
            set_boundary_conditions(side, u0, v0);
            write_file();
        }
    }
}
int main(){
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    cout<<"1"<<"\n";
    get_speed();
    // cout<<axis<<" "<<side<<" "<<left_lim<<" "<<right<<" "<<u0<<" "<<v0<<endl;
    int i = 0;
    get_source();
    nodes_coor();
    cout<<"2"<<"\n";
    cell_to_node();
    center_coordinates_generator();
    cout<<"3"<<"\n";
    face_links();
    cout<<"4"<<"\n";
    cout<<"Bound "<<boundary_face.size()<<endl;
    vector<double> v1, v2;
    double A, V;
    // cout<<face_num;
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
        V_faces.push_back(0);
    }
    for(i=0; i<center_num; i++){
        V = volume(i);
        f_volume[i] = V;
        U_centers_new[i]=0;
        U_centers_old[i]=0;
        P_centers_new[i]=0;
        P_centers_old[i]=0;
        V_centers_new[i]=0;
        V_centers_old[i]=0;
        C_centers_new[i]=0;
        C_centers_old[i]=0;
    }
    for(i=0; i<node_num; i++){
        U_nodes[i] = 0;
        P_nodes[i] = 0;
        V_nodes[i] = 0;
        C_nodes[i] = 0;
    }
    for(i=0; i<speed_over_time.size(); i++){
        fvm(speed_over_time[i][0], speed_over_time[i][1], speed_over_time[i][2]);
    }
    return 0;
}